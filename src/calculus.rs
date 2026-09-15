//! Stages 5–6: symbolic differentiation with a step trace, partial
//! derivatives, gradient, Jacobian and Hessian.

use num_rational::BigRational as Q;
use num_traits::One;

use crate::poly::{atom_has_var, raw_term};
use crate::simplify::simplify;
use crate::solve::{SolveError, Step};
use crate::sym::{func, half, Atom, Mono, Sym};

const MAX_STEPS: usize = 40;

pub struct Differ<'a> {
    pub x: &'a str,
    pub steps: Vec<Step>,
}

impl<'a> Differ<'a> {
    pub fn new(x: &'a str) -> Self {
        Differ { x, steps: Vec::new() }
    }

    fn note(&mut self, rule: &str, line: String) {
        if self.steps.len() < MAX_STEPS {
            self.steps.push(Step { description: rule.to_string(), snapshot: vec![line] });
        }
    }

    pub fn diff(&mut self, s: &Sym) -> Result<Sym, SolveError> {
        let mut out = Sym::zero();
        for (m, c) in &s.terms {
            out.absorb(self.term(c, m)?);
        }
        Ok(out)
    }

    fn term(&mut self, c: &Q, m: &Mono) -> Result<Sym, SolveError> {
        let x = self.x;
        let (dep, indep): (Vec<_>, Vec<_>) = m.iter().partition(|(a, _)| atom_has_var(a, x));
        if dep.is_empty() {
            return Ok(Sym::zero());
        }
        let k = raw_term(indep.into_iter().map(|(a, e)| (a.clone(), e.clone())).collect(), c.clone());
        let mut total = Sym::zero();
        for i in 0..dep.len() {
            let d = self.power(dep[i].0, dep[i].1)?;
            let others: Mono = dep
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != i)
                .map(|(_, (a, e))| ((*a).clone(), (*e).clone()))
                .collect();
            total.absorb(k.mul(&Sym::term(Q::one(), others)).mul(&d));
        }
        if dep.len() > 1 {
            let t = raw_term(m.clone(), c.clone());
            self.note("Product rule", format!("d/d{x} [{t}] = {total}"));
        }
        Ok(total)
    }

    fn power(&mut self, a: &Atom, e: &Q) -> Result<Sym, SolveError> {
        let inner = self.atom(a)?;
        if e.is_one() {
            return Ok(inner);
        }
        let base = match a {
            Atom::Group(b) => (**b).clone(),
            other => Sym::atom(other.clone()),
        };
        let r = base.pow_q(&(e - Q::one()))?.scale(e).mul(&inner);
        let label = if matches!(a, Atom::Var(_)) { "Power rule" } else { "Power rule with chain rule" };
        let x = self.x;
        self.note(label, format!("d/d{x} [{}] = {r}", base.pow_q(e)?));
        Ok(r)
    }

    fn atom(&mut self, a: &Atom) -> Result<Sym, SolveError> {
        let x = self.x;
        Ok(match a {
            Atom::Var(n) => Sym::int(if n == x { 1 } else { 0 }),
            Atom::Group(b) => self.diff(b)?,
            Atom::Pow(b, y) => {
                // d(b^y) = b^y · (y'·ln(b) + y·b'/b)
                let by = Sym::atom(a.clone());
                let (db, dy) = (self.diff(b)?, self.diff(y)?);
                let mut inner = Sym::zero();
                if !dy.is_zero() {
                    inner.absorb(dy.mul(&func("ln", vec![(**b).clone()])?));
                }
                if !db.is_zero() {
                    inner.absorb(y.mul(&db).div(b)?);
                }
                let r = by.mul(&inner);
                self.note("Exponential rule", format!("d/d{x} [{by}] = {r}"));
                r
            }
            Atom::Func(name, args) => {
                let du = self.diff(&args[0])?;
                if du.is_zero() {
                    return Ok(Sym::zero());
                }
                let outer = outer_derivative(name, &args[0])?;
                let r = outer.mul(&du);
                let f = Sym::atom(a.clone());
                if du.is_one() {
                    self.note(&format!("Derivative of {name}"), format!("d/d{x} [{f}] = {r}"));
                } else {
                    self.note("Chain rule", format!("d/d{x} [{f}] = ({outer})·({du}) = {r}"));
                }
                r
            }
            _ => Sym::zero(),
        })
    }
}

/// f'(u) for a function f applied to u.
fn outer_derivative(name: &str, u: &Sym) -> Result<Sym, SolveError> {
    let one = Sym::int(1);
    let f = |n: &str| func(n, vec![u.clone()]);
    let u2 = u.powi(2)?;
    Ok(match name {
        "sin" => f("cos")?,
        "cos" => f("sin")?.neg(),
        "tan" => one.add(&f("tan")?.powi(2)?),
        "asin" => one.sub(&u2).pow_q(&-half())?,
        "acos" => one.sub(&u2).pow_q(&-half())?.neg(),
        "atan" => one.add(&u2).recip()?,
        "sinh" => f("cosh")?,
        "cosh" => f("sinh")?,
        "tanh" => one.sub(&f("tanh")?.powi(2)?),
        "ln" => u.recip()?,
        "log" => u.mul(&func("ln", vec![Sym::int(10)])?).recip()?,
        "log2" => u.mul(&func("ln", vec![Sym::int(2)])?).recip()?,
        "abs" => u.div(&f("abs")?)?,
        _ => {
            return Err(SolveError::Unsupported(format!(
                "{name}() can't be differentiated symbolically"
            )))
        }
    })
}

/// n-th derivative with respect to x, simplified, with the rules used.
pub fn derivative(s: &Sym, x: &str, order: usize) -> Result<(Sym, Vec<Step>), SolveError> {
    let mut d = Differ::new(x);
    let mut cur = s.clone();
    for _ in 0..order {
        cur = simplify(&d.diff(&cur)?);
    }
    Ok((cur, d.steps))
}

pub fn partial(s: &Sym, x: &str) -> Result<Sym, SolveError> {
    Ok(derivative(s, x, 1)?.0)
}

pub fn gradient(f: &Sym, vars: &[String]) -> Result<Vec<Sym>, SolveError> {
    vars.iter().map(|v| partial(f, v)).collect()
}

pub fn jacobian(fs: &[Sym], vars: &[String]) -> Result<Vec<Vec<Sym>>, SolveError> {
    fs.iter().map(|f| gradient(f, vars)).collect()
}

pub fn hessian(f: &Sym, vars: &[String]) -> Result<Vec<Vec<Sym>>, SolveError> {
    gradient(f, vars)?.iter().map(|g| gradient(g, vars)).collect()
}
