//! Stage 7: symbolic integration.
//!
//! Not a general algorithm: linearity, a table of standard forms (with linear
//! inner arguments), rational functions with distinct linear or a single
//! quadratic denominator, e^(ax)·sin/cos(bx), u-substitution and integration
//! by parts, tried in that order. Every antiderivative is checked by
//! differentiating it numerically at sample points; if nothing applies or the
//! check fails, it says so instead of guessing.

use std::collections::HashMap;

use num_complex::Complex64;
use num_rational::BigRational as Q;
use num_traits::{One, Signed, Zero};

use crate::ast::Expr;
use crate::calculus::Differ;
use crate::commands::{var_arg, with_approx};
use crate::error::Error;
use crate::poly::{atom_has_var, coeffs_in, from_uni, raw_term, uni_q};
use crate::simplify::simplify;
use crate::solve::{SolveError, Step};
use crate::sym::{from_expr, func, half, q, rational_power, Atom, Konst, Mono, Sym};
use crate::upoly;
use crate::value::Value;
use crate::Outcome;

const MAX_DEPTH: usize = 6;

type R = Result<Sym, String>;

struct Integrator {
    x: String,
    steps: Vec<Step>,
}

fn e_sym() -> Sym {
    Sym::atom(Atom::Const(Konst::E))
}

fn ln(s: Sym) -> R {
    func("ln", vec![s]).map_err(|e| e.to_string())
}

fn ln_abs(s: Sym) -> R {
    ln(func("abs", vec![s]).map_err(|e| e.to_string())?)
}

fn f1(name: &str, s: Sym) -> R {
    func(name, vec![s]).map_err(|e| e.to_string())
}

fn ev<T>(r: Result<T, crate::error::EvalError>) -> Result<T, String> {
    r.map_err(|e| e.to_string())
}

fn derivative(s: &Sym, x: &str) -> R {
    Differ::new(x).diff(s).map(|d| simplify(&d)).map_err(|e| e.to_string())
}

impl Integrator {
    fn note(&mut self, description: String, line: String) {
        if self.steps.len() < 40 {
            self.steps.push(Step { description, snapshot: vec![line] });
        }
    }

    /// Slope a if u = a·x + b (b free of x).
    fn slope(&self, u: &Sym) -> Option<Q> {
        let c = coeffs_in(u, &self.x)?;
        if c.len() == 2 {
            c[1].as_constant()
        } else {
            None
        }
    }

    fn integrate(&mut self, f: &Sym, depth: usize) -> R {
        let x = self.x.clone();
        if depth > MAX_DEPTH {
            return Err("too many nested rules".into());
        }
        if !f.has_var(&x) {
            return Ok(f.mul(&Sym::var(&x)));
        }
        if f.terms.len() > 1 {
            if depth == 0 {
                self.note("Integrate term by term (linearity)".into(), format!("∫ ({f}) d{x}"));
            }
            let mut total = Sym::zero();
            for (m, c) in &f.terms {
                total.absorb(self.integrate(&raw_term(m.clone(), c.clone()), depth)?);
            }
            return Ok(total);
        }
        let (m, c) = f.single().unwrap();
        let (dep, indep): (Mono, Mono) = m.iter().map(|(a, e)| (a.clone(), e.clone())).partition(|(a, _)| atom_has_var(a, &x));
        let k = raw_term(indep, c.clone());
        let p = raw_term(dep, Q::one());
        Ok(k.mul(&self.product(&p, depth)?))
    }

    fn product(&mut self, p: &Sym, depth: usize) -> R {
        let (m, _) = p.single().unwrap();
        let factors: Vec<(Atom, Q)> = m.iter().map(|(a, e)| (a.clone(), e.clone())).collect();
        if factors.len() == 1 {
            if let Some(r) = self.arcsin(&factors[0].0, &factors[0].1)? {
                return Ok(r);
            }
            if let Some(r) = self.table(&factors[0].0, &factors[0].1)? {
                return Ok(r);
            }
        }
        if let Some(r) = self.rational(p)? {
            return Ok(r);
        }
        if let Some(r) = self.exp_trig(&factors)? {
            return Ok(r);
        }
        if let Some(r) = self.substitution(p, depth)? {
            return Ok(r);
        }
        if let Some(r) = self.by_parts(&factors, depth)? {
            return Ok(r);
        }
        Err(format!("no rule applies to {p}"))
    }

    /// ∫ (a − b·x²)^(−1/2) dx = asin(x·√(b/a)) / √b for a, b > 0.
    fn arcsin(&mut self, at: &Atom, e: &Q) -> Result<Option<Sym>, String> {
        let x = self.x.clone();
        let Atom::Group(g) = at else { return Ok(None) };
        if *e != -half() {
            return Ok(None);
        }
        let Some(c) = coeffs_in(g, &x) else { return Ok(None) };
        if c.len() != 3 || !c[1].is_zero() {
            return Ok(None);
        }
        let (Some(a), Some(b)) = (c[0].as_constant(), c[2].as_constant().map(|v| -v)) else { return Ok(None) };
        if !a.is_positive() || !b.is_positive() {
            return Ok(None);
        }
        let k = rational_power(&(&b / &a), &half());
        let r = ev(f1("asin", Sym::var(&x).mul(&k))?.div(&rational_power(&b, &half())))?;
        self.note("Standard integral ∫ 1/√(a − b·x²) dx = asin(x·√(b/a))/√b".into(), format!("= {r}"));
        Ok(Some(r))
    }

    /// Standard forms a^e with a linear inner argument.
    fn table(&mut self, a: &Atom, e: &Q) -> Result<Option<Sym>, String> {
        let x = self.x.clone();
        let one = Q::one();
        let (u, is_var) = match a {
            Atom::Var(_) => (Sym::var(&x), true),
            Atom::Group(b) => ((**b).clone(), false),
            Atom::Func(_, args) => (args[0].clone(), false),
            Atom::Pow(_, y) => ((**y).clone(), false),
            _ => return Ok(None),
        };
        let Some(slope) = self.slope(&u) else { return Ok(None) };
        let inv = slope.recip();
        let r = match a {
            Atom::Var(_) | Atom::Group(_) => {
                if *e == -one.clone() {
                    ln_abs(u.clone())?.scale(&inv)
                } else {
                    let e1 = e + &one;
                    ev(u.pow_q(&e1))?.scale(&(e1.recip() * &inv))
                }
            }
            Atom::Func(name, _) => match (name.as_str(), e) {
                (n, e) if e.is_one() => match n {
                    "sin" => f1("cos", u.clone())?.neg(),
                    "cos" => f1("sin", u.clone())?,
                    "tan" => ln_abs(f1("cos", u.clone())?)?.neg(),
                    "sinh" => f1("cosh", u.clone())?,
                    "cosh" => f1("sinh", u.clone())?,
                    "tanh" => ln(f1("cosh", u.clone())?)?,
                    "ln" => u.mul(&ln(u.clone())?).sub(&u),
                    "log" => u.mul(&ln(u.clone())?).sub(&u).div(&ln(Sym::int(10))?).map_err(|e| e.to_string())?,
                    "asin" => u.mul(&f1("asin", u.clone())?).add(&ev(Sym::int(1).sub(&ev(u.powi(2))?).pow_q(&half()))?),
                    "acos" => u.mul(&f1("acos", u.clone())?).sub(&ev(Sym::int(1).sub(&ev(u.powi(2))?).pow_q(&half()))?),
                    "atan" => u.mul(&f1("atan", u.clone())?).sub(&ln(Sym::int(1).add(&ev(u.powi(2))?))?.scale(&half())),
                    "abs" => u.mul(&f1("abs", u.clone())?).scale(&half()),
                    _ => return Ok(None),
                }
                .scale(&inv),
                ("sin", e) if *e == q(2) => {
                    u.scale(&half()).sub(&f1("sin", u.scale(&q(2)))?.scale(&Q::new(1.into(), 4.into()))).scale(&inv)
                }
                ("cos", e) if *e == q(2) => {
                    u.scale(&half()).add(&f1("sin", u.scale(&q(2)))?.scale(&Q::new(1.into(), 4.into()))).scale(&inv)
                }
                ("tan", e) if *e == q(2) => f1("tan", u.clone())?.sub(&u).scale(&inv),
                ("cos", e) if *e == q(-2) => f1("tan", u.clone())?.scale(&inv),
                ("sin", e) if *e == q(-2) => {
                    ev(f1("cos", u.clone())?.div(&f1("sin", u.clone())?))?.neg().scale(&inv)
                }
                _ => return Ok(None),
            },
            Atom::Pow(b, _) if e.is_one() && !b.has_var(&x) => {
                // ∫ b^u dx = b^u / (u'·ln b)
                let bu = Sym::atom(a.clone());
                ev(bu.div(&ln((**b).clone())?))?.scale(&inv)
            }
            _ => return Ok(None),
        };
        let what = match a {
            Atom::Var(_) if *e == -Q::one() => "∫ 1/x dx = ln|x|".to_string(),
            Atom::Var(_) => "Power rule: ∫ x^n dx = x^(n+1)/(n+1)".to_string(),
            _ if is_var || slope.is_one() => "Standard integral".to_string(),
            _ => format!("Standard integral with linear argument (divide by {slope})"),
        };
        let integrand = Sym::term(Q::one(), Mono::from([(a.clone(), e.clone())]));
        self.note(what, format!("∫ {integrand} d{x} = {r}"));
        Ok(Some(r))
    }

    /// Rational functions N(x)/D(x) whose denominator factors into distinct
    /// linear factors, or is a single irreducible quadratic.
    fn rational(&mut self, p: &Sym) -> Result<Option<Sym>, String> {
        let x = self.x.clone();
        let (m, _) = p.single().unwrap();
        let mut num = Sym::int(1);
        let mut den = Sym::int(1);
        for (a, e) in m {
            match a {
                Atom::Var(_) if e.is_integer() => {
                    let t = Sym::term(Q::one(), Mono::from([(a.clone(), e.clone())]));
                    if e.is_positive() { num = num.mul(&t) } else { den = den.mul(&ev(t.recip())?) }
                }
                Atom::Group(b) if e.is_integer() && e.is_negative() && uni_q(b, &x).is_some() => {
                    den = den.mul(&ev(b.powi(-e.to_integer().try_into().unwrap_or(1)))?);
                }
                _ => return Ok(None),
            }
        }
        let (Some(np), Some(dp)) = (uni_q(&num, &x), uni_q(&den, &x)) else { return Ok(None) };
        if upoly::degree(&dp) < 2 {
            return Ok(None); // table handles 1/(ax+b)^k
        }
        let (quot, rem) = upoly::divrem(&np, &dp);
        let mut result = self.integrate(&from_uni(&quot, &x), MAX_DEPTH)?;
        let fac = upoly::factor(&dp);
        let d_prime = upoly::derivative(&dp);
        if fac.factors.iter().all(|(f, k)| upoly::degree(f) == 1 && *k == 1) {
            // Σ R(r)/D'(r) · ln|x − r|
            let mut parts = Vec::new();
            for (f, _) in &fac.factors {
                let r = -&f[0] / &f[1];
                let coeff = upoly::eval_q(&rem, &r) / upoly::eval_q(&d_prime, &r);
                parts.push(format!("{coeff}/({x} - {r})"));
                result.absorb(ln_abs(Sym::var(&x).sub(&Sym::constant(r)))?.scale(&coeff));
            }
            self.note("Partial fractions".into(), parts.join(" + "));
            return Ok(Some(result));
        }
        if fac.factors.len() == 1 && fac.factors[0].1 == 1 && upoly::degree(&dp) == 2 && rem.len() <= 2 {
            let (c0, b, a) = (&dp[0], &dp[1], &dp[2]);
            let bb = rem.get(1).cloned().unwrap_or_else(Q::zero);
            let cc = rem.first().cloned().unwrap_or_else(Q::zero);
            let qx = from_uni(&dp, &x);
            // (B/2a)·ln|q| + (C − B·b/2a)·∫ 1/q
            let two_a = a * q(2);
            result.absorb(ln_abs(qx)?.scale(&(&bb / &two_a)));
            let k = &cc - &bb * b / &two_a;
            let delta = a * c0 * q(4) - b * b;
            let lin = Sym::var(&x).scale(&two_a).add(&Sym::constant(b.clone()));
            let inner = if delta.is_positive() {
                let s = rational_power(&delta, &half());
                ev(f1("atan", ev(lin.div(&s))?)?.scale(&q(2)).div(&s))?
            } else {
                let s = rational_power(&-&delta, &half());
                ev(ln_abs(ev(lin.sub(&s).div(&lin.add(&s)))?)?.div(&s))?
            };
            result.absorb(inner.scale(&k));
            self.note("Quadratic denominator: split into ln and arctan parts".into(), format!("∫ {p} d{x}"));
            return Ok(Some(result));
        }
        Ok(None)
    }

    /// e^(a·x+b) · sin/cos(c·x+d).
    fn exp_trig(&mut self, factors: &[(Atom, Q)]) -> Result<Option<Sym>, String> {
        if factors.len() != 2 || !factors.iter().all(|(_, e)| e.is_one()) {
            return Ok(None);
        }
        let (mut ex, mut tr) = (None, None);
        for (a, _) in factors {
            match a {
                Atom::Pow(b, u) if **b == e_sym() => ex = Some((**u).clone()),
                Atom::Func(n, args) if n == "sin" || n == "cos" => tr = Some((n.clone(), args[0].clone())),
                _ => {}
            }
        }
        let (Some(u), Some((name, v))) = (ex, tr) else { return Ok(None) };
        let (Some(a), Some(b)) = (self.slope(&u), self.slope(&v)) else { return Ok(None) };
        let (s, c) = (f1("sin", v.clone())?, f1("cos", v.clone())?);
        let comb = if name == "sin" {
            s.scale(&a).sub(&c.scale(&b))
        } else {
            c.scale(&a).add(&s.scale(&b))
        };
        let r = ev(e_sym().pow(&u))?.mul(&comb).scale(&(&a * &a + &b * &b).recip());
        self.note("Standard form ∫ e^(ax)·sin/cos(bx) (integrate by parts twice)".into(), format!("= {r}"));
        Ok(Some(r))
    }

    /// u-substitution: try u = each inner expression g with q = f/g' free of x.
    fn substitution(&mut self, p: &Sym, depth: usize) -> Result<Option<Sym>, String> {
        let x = self.x.clone();
        let (m, _) = p.single().unwrap();
        let mut candidates: Vec<Sym> = Vec::new();
        for a in m.keys() {
            match a {
                Atom::Func(_, args) => {
                    candidates.push(args[0].clone());
                    // The function itself: ∫ sin·cos (u = sin), ∫ ln(x)/x (u = ln x).
                    candidates.push(Sym::atom(a.clone()));
                }
                Atom::Group(b) => candidates.push((**b).clone()),
                Atom::Pow(b, y) => {
                    candidates.push((**y).clone());
                    candidates.push((**b).clone());
                }
                _ => {}
            }
        }
        // Also powers of x inside functions, e.g. x^2 in x·e^(x^2).
        let u_name = if p.has_var("u") { "u_1" } else { "u" };
        for g in candidates {
            if !g.has_var(&x) || g == Sym::var(&x) || self.slope(&g).is_some() {
                continue;
            }
            let Ok(gp) = derivative(&g, &x) else { continue };
            if gp.is_zero() {
                continue;
            }
            let Ok(ratio) = p.div(&gp) else { continue };
            let ratio = simplify(&ratio);
            let Some(qu) = replace(&ratio, &g, u_name, &x) else { continue };
            if qu.has_var(&x) {
                continue;
            }
            let mut sub = Integrator { x: u_name.to_string(), steps: vec![] };
            let Ok(iu) = sub.integrate(&qu, depth + 1) else { continue };
            let Ok(back) = iu.subst(u_name, &g) else { continue };
            self.note(
                format!("Substitute {u_name} = {g}, d{u_name} = ({gp}) d{x}"),
                format!("∫ {qu} d{u_name} = {iu}"),
            );
            self.steps.extend(sub.steps);
            return Ok(Some(back));
        }
        Ok(None)
    }

    /// Integration by parts, choosing u by LIATE.
    fn by_parts(&mut self, factors: &[(Atom, Q)], depth: usize) -> Result<Option<Sym>, String> {
        let x = self.x.clone();
        let priority = |a: &Atom, e: &Q| -> u8 {
            match a {
                Atom::Func(n, _) if matches!(n.as_str(), "ln" | "log" | "log2") && e.is_one() => 5,
                Atom::Func(n, _) if matches!(n.as_str(), "asin" | "acos" | "atan") && e.is_one() => 4,
                Atom::Var(_) if e.is_integer() && e.is_positive() => 3,
                _ => 0,
            }
        };
        let Some(i) = (0..factors.len()).max_by_key(|&i| priority(&factors[i].0, &factors[i].1)) else {
            return Ok(None);
        };
        if priority(&factors[i].0, &factors[i].1) == 0 || factors.len() < 2 {
            return Ok(None);
        }
        let u = Sym::term(Q::one(), Mono::from([factors[i].clone()]));
        let dv: Mono = factors.iter().enumerate().filter(|(j, _)| *j != i).map(|(_, f)| f.clone()).collect();
        let dv = Sym::term(Q::one(), dv);
        let Ok(v) = self.integrate(&dv, depth + 1) else { return Ok(None) };
        let du = derivative(&u, &x)?;
        let rest = simplify(&v.mul(&du));
        let Ok(w) = self.integrate(&rest, depth + 1) else { return Ok(None) };
        self.note(
            "Integration by parts: ∫ u dv = u·v − ∫ v du".into(),
            format!("u = {u}, dv = {dv} d{x}, v = {v}"),
        );
        Ok(Some(u.mul(&v).sub(&w)))
    }
}

/// Rewrite `s` in terms of u where g appears; None if it can't.
fn replace(s: &Sym, g: &Sym, u: &str, x: &str) -> Option<Sym> {
    // g = c·x^k: replace x^(j·k) by (u/c)^j.
    if let Some((gm, gc)) = g.single() {
        if gm.len() == 1 {
            if let Some((Atom::Var(n), k)) = gm.iter().next() {
                if n == x {
                    let mut out = Sym::zero();
                    for (m, c) in &s.terms {
                        let mut t = Sym::constant(c.clone());
                        for (a, e) in m {
                            let piece = match a {
                                Atom::Var(n) if n == x => {
                                    let j = e / k;
                                    Sym::var(u).scale(&gc.recip()).pow_q(&j).ok()?
                                }
                                _ => replace_atom(a, g, u, x)?.pow_q(e).ok()?,
                            };
                            t = t.mul(&piece);
                        }
                        out.absorb(t);
                    }
                    return Some(out);
                }
            }
        }
    }
    let mut out = Sym::zero();
    for (m, c) in &s.terms {
        let mut t = Sym::constant(c.clone());
        for (a, e) in m {
            t = t.mul(&replace_atom(a, g, u, x)?.pow_q(e).ok()?);
        }
        out.absorb(t);
    }
    Some(out)
}

fn replace_atom(a: &Atom, g: &Sym, u: &str, x: &str) -> Option<Sym> {
    let as_sym = match a {
        Atom::Group(b) => (**b).clone(),
        other => Sym::atom(other.clone()),
    };
    if as_sym == *g {
        return Some(Sym::var(u));
    }
    match a {
        Atom::Func(name, args) => {
            let args = args.iter().map(|s| replace(s, g, u, x)).collect::<Option<Vec<_>>>()?;
            func(name, args).ok()
        }
        Atom::Group(b) => replace(b, g, u, x),
        Atom::Pow(b, y) => replace(b, g, u, x)?.pow(&replace(y, g, u, x)?).ok(),
        _ => Some(as_sym),
    }
}

fn eval_at(s: &Sym, x: &str, t: f64, others: &HashMap<String, Value>) -> Option<Complex64> {
    let mut env = others.clone();
    env.insert(x.to_string(), Value::Real(t));
    let v = s.eval(&env).ok()?.to_c64();
    v.is_finite().then_some(v)
}

/// Check d/dx F = f at sample points. Err if they disagree anywhere.
fn verify(f: &Sym, big_f: &Sym, x: &str) -> Result<bool, String> {
    let d = derivative(big_f, x)?;
    let others: HashMap<String, Value> = f
        .vars()
        .into_iter()
        .filter(|v| v != x)
        .enumerate()
        .map(|(i, v)| (v, Value::Real(0.61 + 0.37 * i as f64)))
        .collect();
    let mut checked = 0;
    for t in [0.37, 0.71, 1.29, 2.13, -0.43, -1.66, 3.7] {
        if let (Some(a), Some(b)) = (eval_at(f, x, t, &others), eval_at(&d, x, t, &others)) {
            if (a - b).norm() > 1e-6 * (1.0 + a.norm()) {
                return Err(format!("check failed at {x} = {t}"));
            }
            checked += 1;
        }
    }
    Ok(checked > 0)
}

/// Adaptive Simpson on a real integrand; None if it's undefined somewhere.
fn simpson(f: &dyn Fn(f64) -> Option<f64>, a: f64, b: f64) -> Option<f64> {
    fn rec(f: &dyn Fn(f64) -> Option<f64>, a: f64, b: f64, fa: f64, fm: f64, fb: f64, whole: f64, eps: f64, depth: u32) -> Option<f64> {
        let m = (a + b) / 2.0;
        let (lm, rm) = ((a + m) / 2.0, (m + b) / 2.0);
        let (flm, frm) = (f(lm)?, f(rm)?);
        let left = (m - a) / 6.0 * (fa + 4.0 * flm + fm);
        let right = (b - m) / 6.0 * (fm + 4.0 * frm + fb);
        if depth == 0 || (left + right - whole).abs() <= 15.0 * eps {
            return Some(left + right + (left + right - whole) / 15.0);
        }
        Some(rec(f, a, m, fa, flm, fm, left, eps / 2.0, depth - 1)? + rec(f, m, b, fm, frm, fb, right, eps / 2.0, depth - 1)?)
    }
    let (fa, fb, fm) = (f(a)?, f(b)?, f((a + b) / 2.0)?);
    let whole = (b - a) / 6.0 * (fa + 4.0 * fm + fb);
    rec(f, a, b, fa, fm, fb, whole, 1e-10, 40)
}

fn unsupported(msg: String) -> Error {
    Error::Solve(SolveError::Unsupported(msg))
}

pub fn command(args: &[Expr]) -> Result<Outcome, Error> {
    if !matches!(args.len(), 1 | 2 | 4) {
        return Err(unsupported("usage: integrate(f), integrate(f, x), or integrate(f, x, a, b)".into()));
    }
    let f = simplify(&from_expr(&args[0])?);
    let x = var_arg(args.get(1), &f, "integrate")?;
    let mut it = Integrator { x: x.clone(), steps: vec![] };
    let found = it.integrate(&f, 0).map(|r| simplify(&r));
    let found = match found {
        Ok(big_f) => match verify(&f, &big_f, &x) {
            Ok(_) => Ok(big_f),
            Err(why) => Err(format!("the candidate antiderivative {big_f} didn't verify ({why})")),
        },
        Err(e) => Err(e),
    };
    let mut steps = it.steps;

    if args.len() <= 2 {
        let big_f = found.map_err(|why| {
            unsupported(format!(
                "couldn't integrate {f} symbolically: {why}. Integration has no general algorithm; \
                 this covers standard forms, substitution, by parts and simple rational functions"
            ))
        })?;
        steps.push(Step {
            description: "Check: differentiating the result gives back the integrand".into(),
            snapshot: vec![format!("d/d{x} [{big_f}] = {f}")],
        });
        return Ok(Outcome::Text { steps, answer: format!("∫ {f} d{x} = {big_f} + C") });
    }

    // Definite integral.
    let (a, b) = (simplify(&from_expr(&args[2])?), simplify(&from_expr(&args[3])?));
    let num = |s: &Sym| -> Result<f64, Error> {
        match s.eval(&HashMap::new())? {
            v => v.as_real_f64().ok_or_else(|| unsupported("integration bounds must be real numbers".into())),
        }
    };
    let (af, bf) = (num(&a)?, num(&b)?);
    if f.vars().iter().any(|v| *v != x) {
        return Err(unsupported("a definite integral can only involve the integration variable".into()));
    }
    let g = |t: f64| -> Option<f64> {
        let v = eval_at(&f, &x, t, &HashMap::new())?;
        (v.im.abs() < 1e-12).then_some(v.re)
    };
    let numeric = simpson(&g, af, bf);
    let exact = found.ok().and_then(|big_f| {
        let fb = simplify(&big_f.subst(&x, &b).ok()?);
        let fa = simplify(&big_f.subst(&x, &a).ok()?);
        Some((big_f, fb.sub(&fa)))
    });
    let heading = format!("∫ from {a} to {b} of {f} d{x}");
    match (exact, numeric) {
        (Some((big_f, val)), Some(n)) => {
            let v = val.eval(&HashMap::new())?.to_c64();
            if (v.re - n).abs() <= 1e-6 * (1.0 + n.abs()) && v.im.abs() < 1e-9 {
                steps.push(Step {
                    description: "Evaluate the antiderivative at the bounds: F(b) − F(a)".into(),
                    snapshot: vec![format!("F({x}) = {big_f}")],
                });
                Ok(Outcome::Text { steps, answer: format!("{heading} = {}", with_approx(&val)) })
            } else {
                Ok(Outcome::Text {
                    steps: vec![],
                    answer: format!(
                        "{heading} ≈ {}\n(numerical; the antiderivative doesn't apply across the whole interval)",
                        crate::format::fmt_f64(n)
                    ),
                })
            }
        }
        (None, Some(n)) => Ok(Outcome::Text {
            steps: vec![],
            answer: format!("{heading} ≈ {}\n(numerical; no symbolic antiderivative found)", crate::format::fmt_f64(n)),
        }),
        (_, None) => Err(unsupported(
            "the integrand is undefined or complex somewhere on the interval; improper integrals aren't supported".into(),
        )),
    }
}
