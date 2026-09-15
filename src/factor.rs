//! Stage 8: factoring polynomials over the rationals.
//!
//! Common factors first; then univariate polynomials are split with the
//! rational root theorem and Kronecker's method, and polynomials in two
//! variables whose terms all have the same total degree are factored through
//! their one-variable version. Anything else is reported as not fully factored.

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational as Q;
use num_traits::{One, Signed, ToPrimitive, Zero};

use crate::poly::{divide, from_uni, is_polynomial, leading, raw_term, uni_q};
use crate::solve::{SolveError, Step};
use crate::sym::{Atom, Mono, Sym};
use crate::upoly;

#[derive(Debug, Clone, PartialEq)]
pub struct Factorization {
    pub content: Q,
    pub factors: Vec<(Sym, usize)>,
    /// False when a search limit was hit or the case isn't supported, so a
    /// factor may still be reducible.
    pub complete: bool,
}

impl Factorization {
    pub fn to_sym(&self) -> Sym {
        self.factors.iter().fold(Sym::constant(self.content.clone()), |acc, (f, k)| {
            acc.mul(&f.powi(*k as i64).expect("positive power"))
        })
    }

    pub fn display(&self) -> String {
        let parts: Vec<String> = self
            .factors
            .iter()
            .map(|(f, k)| {
                let body = if f.terms.len() > 1 { format!("({f})") } else { f.to_string() };
                if *k > 1 {
                    format!("{body}^{k}")
                } else {
                    body
                }
            })
            .collect();
        let c = &self.content;
        if parts.is_empty() {
            return c.to_string();
        }
        if c.is_one() && self.factors.len() == 1 && self.factors[0].1 == 1 {
            return self.factors[0].0.to_string();
        }
        let prefix = if c.is_one() {
            String::new()
        } else if *c == -Q::one() {
            "-".to_string()
        } else if c.is_integer() {
            c.to_string()
        } else {
            format!("({c})")
        };
        format!("{prefix}{}", parts.join(""))
    }
}

fn unsupported() -> SolveError {
    SolveError::Unsupported("factor() works on polynomials with rational coefficients".into())
}

fn is_homogeneous(s: &Sym) -> Option<Q> {
    let mut degs = s.terms.keys().map(|m| m.values().cloned().sum::<Q>());
    let first = degs.next()?;
    degs.all(|d| d == first).then_some(first)
}

pub fn factor(s: &Sym) -> Result<(Factorization, Vec<Step>), SolveError> {
    if !is_polynomial(s) || s.terms.keys().any(|m| m.keys().any(|a| !matches!(a, Atom::Var(_)))) {
        return Err(unsupported());
    }
    let mut steps = Vec::new();
    if s.is_zero() {
        return Ok((Factorization { content: Q::zero(), factors: vec![], complete: true }, steps));
    }

    // Common numeric and monomial factor.
    let num_gcd = s.terms.values().fold(BigInt::zero(), |g, c| g.gcd(c.numer()));
    let den_lcm = s.terms.values().fold(BigInt::one(), |l, c| l.lcm(c.denom()));
    let sign = if leading(s).is_some_and(|(_, c)| c.is_negative()) { -1 } else { 1 };
    let mut content = Q::new(num_gcd * sign, den_lcm);
    let mut common: Mono = s.terms.keys().next().unwrap().clone();
    for m in s.terms.keys() {
        common = common.into_iter().filter_map(|(a, e)| m.get(&a).map(|e2| (a, e.min(e2.clone())))).collect();
    }
    let divisor = raw_term(common.clone(), content.clone());
    let rest = divide(s, &divisor).ok_or_else(unsupported)?;
    let mut factors: Vec<(Sym, usize)> = common
        .iter()
        .map(|(a, e)| (Sym::atom(a.clone()), e.to_integer().to_usize().unwrap_or(1)))
        .collect();
    if !divisor.is_one() && rest.terms.len() > 1 {
        steps.push(Step {
            description: "Take out the common factor".into(),
            snapshot: vec![format!("{s} = {divisor}·({rest})")],
        });
    }

    let vars: Vec<String> = rest.vars().into_iter().collect();
    let mut complete = true;
    let mut pieces: Vec<(Sym, usize)> = Vec::new();
    match vars.len() {
        0 => content *= rest.as_constant().unwrap_or_else(Q::one),
        1 => {
            let x = &vars[0];
            let f = upoly::factor(&uni_q(&rest, x).ok_or_else(unsupported)?);
            content *= &f.content;
            complete = f.complete;
            pieces = f.factors.iter().map(|(p, k)| (from_uni(p, x), *k)).collect();
            describe_univariate(&f, x, &mut steps);
        }
        2 if is_homogeneous(&rest).is_some() => {
            // f(x, y) = y^d · f(x/y, 1): factor the one-variable version.
            let (x, y) = (&vars[0], &vars[1]);
            let total = is_homogeneous(&rest).unwrap().to_integer().to_usize().unwrap_or(0);
            let uni = uni_q(&rest.subst(y, &Sym::int(1))?, x).ok_or_else(unsupported)?;
            let f = upoly::factor(&uni);
            content *= &f.content;
            complete = f.complete;
            steps.push(Step {
                description: format!("All terms have degree {total}: factor the version with {y} = 1"),
                snapshot: vec![from_uni(&uni, x).to_string()],
            });
            let mut used = 0;
            for (p, k) in &f.factors {
                let d = upoly::degree(p);
                used += d * k;
                let mut h = Sym::zero();
                for (i, c) in p.iter().enumerate() {
                    let m = Mono::from_iter(
                        [(Atom::Var(x.clone()), i), (Atom::Var(y.clone()), d - i)]
                            .into_iter()
                            .filter(|(_, e)| *e > 0)
                            .map(|(a, e)| (a, Q::from_integer(e.into()))),
                    );
                    h.absorb(raw_term(m, c.clone()));
                }
                pieces.push((h, *k));
            }
            if total > used {
                pieces.push((Sym::var(y), total - used));
            }
        }
        _ => {
            complete = false;
            pieces.push((rest.clone(), 1));
        }
    }
    factors.extend(pieces.into_iter().filter(|(f, _)| !f.is_one()));
    let result = Factorization { content, factors, complete };
    steps.push(Step { description: "Result".into(), snapshot: vec![result.display()] });
    Ok((result, steps))
}

fn describe_univariate(f: &upoly::Factored, x: &str, steps: &mut Vec<Step>) {
    let roots: Vec<String> = f
        .factors
        .iter()
        .filter(|(p, _)| upoly::degree(p) == 1)
        .map(|(p, _)| format!("{x} = {}", -&p[0] / &p[1]))
        .collect();
    if !roots.is_empty() {
        steps.push(Step {
            description: "Rational root theorem: test ±(factors of constant)/(factors of leading coefficient)".into(),
            snapshot: vec![format!("roots found: {}", roots.join(", "))],
        });
    }
    let higher: Vec<String> = f
        .factors
        .iter()
        .filter(|(p, _)| upoly::degree(p) >= 2)
        .map(|(p, _)| from_uni(p, x).to_string())
        .collect();
    if !higher.is_empty() {
        let note = if f.complete {
            "Remaining factors have no rational roots and don't split further over the rationals"
        } else {
            "Search limit reached: remaining factors may split further"
        };
        steps.push(Step { description: note.into(), snapshot: higher });
    }
}
