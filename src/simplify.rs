//! Simplification beyond the canonical form: cancel common factors in
//! fractions, e.g. `(x^2-1)/(x-1)` → `x + 1`.

use num_rational::BigRational as Q;
use num_traits::{One, Signed};

use crate::poly::{divide, from_uni, is_polynomial, raw_term, uni_q};
use crate::sym::{Atom, Mono, Sym};
use crate::upoly;

pub fn simplify(s: &Sym) -> Sym {
    let mut cur = s.clone();
    let mut tried: Vec<Sym> = Vec::new();
    for _ in 0..16 {
        let Some((d, k)) = next_denominator(&cur, &tried) else { break };
        tried.push(d.clone());
        if let Some(next) = cancel(&cur, &d, k) {
            cur = next;
            tried.clear();
        }
    }
    pythagorean(&cur)
}

/// Rewrite sin(u)^2 as 1 - cos(u)^2 when cos(u) also appears and doing so
/// makes the expression shorter (e.g. sin(x)^2 + cos(x)^2 = 1).
fn pythagorean(s: &Sym) -> Sym {
    use crate::sym::Atom;
    for (m, c) in &s.terms {
        for (a, e) in m {
            let Atom::Func(name, args) = a else { continue };
            if name != "sin" || !e.is_integer() || *e < crate::sym::q(2) {
                continue;
            }
            let cos_atom = Atom::Func("cos".into(), args.clone());
            if !s.terms.keys().any(|m2| m2.contains_key(&cos_atom)) {
                continue;
            }
            let mut rest = m.clone();
            let e2 = e - crate::sym::q(2);
            if e2 == crate::sym::q(0) {
                rest.remove(a);
            } else {
                rest.insert(a.clone(), e2);
            }
            let Ok(cos2) = Sym::atom(cos_atom).powi(2) else { continue };
            let replaced = crate::poly::raw_term(rest, c.clone()).mul(&Sym::int(1).sub(&cos2));
            let cand = s.sub(&crate::poly::raw_term(m.clone(), c.clone())).add(&replaced);
            if cand.terms.len() < s.terms.len() {
                return pythagorean(&cand);
            }
        }
    }
    s.clone()
}

/// A grouped denominator `(d)^(-k)` not tried yet.
fn next_denominator(s: &Sym, tried: &[Sym]) -> Option<(Sym, i64)> {
    for m in s.terms.keys() {
        for (a, e) in m {
            if let Atom::Group(b) = a {
                if e.is_integer() && e.is_negative() && !tried.contains(b) {
                    return Some(((**b).clone(), -e.to_integer().try_into().ok()?));
                }
            }
        }
    }
    None
}

fn group_pow(d: &Sym, k: i64) -> Sym {
    Sym::term(Q::one(), Mono::from([(Atom::Group(Box::new(d.clone())), Q::from_integer(k.into()))]))
}

fn cancel(s: &Sym, d: &Sym, k: i64) -> Option<Sym> {
    // Numerator N with s = N / d^k.
    let n = s.mul(&raw_term(
        Mono::from([(Atom::Group(Box::new(d.clone())), Q::from_integer(k.into()))]),
        Q::one(),
    ));
    if !is_polynomial(&n) {
        return None;
    }
    if let Some(quot) = divide(&n, d) {
        return Some(quot.mul(&group_pow(d, -(k - 1))));
    }
    // Univariate: cancel the polynomial gcd.
    let vars: Vec<String> = n.vars().union(&d.vars()).cloned().collect();
    if k != 1 || vars.len() != 1 {
        return None;
    }
    let x = &vars[0];
    let (np, dp) = (uni_q(&n, x)?, uni_q(d, x)?);
    let g = upoly::gcd(&np, &dp);
    if upoly::degree(&g) == 0 {
        return None;
    }
    let (n2, _) = upoly::divrem(&np, &g);
    let (d2, _) = upoly::divrem(&dp, &g);
    from_uni(&n2, x).div(&from_uni(&d2, x)).ok()
}
