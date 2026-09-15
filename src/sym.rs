//! Symbolic expressions in canonical form: a sum of terms, each an exact
//! rational coefficient times a product of atoms raised to rational powers.
//!
//! Canonical means equal expressions usually compare equal: `2x+3x` and `5x`
//! are the same `Sym`. Radicals are normalized (`sqrt(8)` = `2sqrt(2)`),
//! `i^2` = `-1`, and integer powers of sums are expanded.

use std::collections::btree_map::Entry;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::f64::consts::{E, PI};

use num_bigint::BigInt;
use num_complex::Complex;
use num_integer::Integer;
use num_rational::BigRational as Q;
use num_traits::{One, Signed, ToPrimitive, Zero};

use crate::ast::Expr;
use crate::builtins::{is_command, FUNCTIONS};
use crate::error::{EvalError, NameKind};
use crate::suggest::suggest;
use crate::value::{rat_ipow, Value};

/// Integer powers of sums up to this are expanded; larger stay grouped.
const MAX_EXPAND_POWER: i64 = 64;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Konst {
    Pi,
    E,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Atom {
    Var(String),
    Const(Konst),
    I,
    /// Positive integer base; only appears with a fractional exponent (a radical).
    Num(BigInt),
    Func(String, Vec<Sym>),
    /// A sub-expression that can't be merged into the term, e.g. `(x+1)^(-1)`.
    Group(Box<Sym>),
    /// `base^exponent` with a non-constant exponent, e.g. `e^x`.
    Pow(Box<Sym>, Box<Sym>),
}

pub type Mono = BTreeMap<Atom, Q>;

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Default)]
pub struct Sym {
    /// Coefficients are never zero; exponents are never zero.
    pub terms: BTreeMap<Mono, Q>,
}

pub(crate) fn q(n: i64) -> Q {
    Q::from_integer(n.into())
}

pub(crate) fn half() -> Q {
    Q::new(1.into(), 2.into())
}

impl Sym {
    pub fn zero() -> Sym {
        Sym::default()
    }

    pub fn constant(c: Q) -> Sym {
        let mut s = Sym::zero();
        if !c.is_zero() {
            s.terms.insert(Mono::new(), c);
        }
        s
    }

    pub fn int(n: i64) -> Sym {
        Sym::constant(q(n))
    }

    pub fn var(name: &str) -> Sym {
        Sym::atom(Atom::Var(name.to_string()))
    }

    pub fn atom(a: Atom) -> Sym {
        normalize_term(Q::one(), Mono::from([(a, Q::one())]))
    }

    pub fn term(c: Q, m: Mono) -> Sym {
        normalize_term(c, m)
    }

    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn as_constant(&self) -> Option<Q> {
        match self.terms.len() {
            0 => Some(Q::zero()),
            1 => {
                let (m, c) = self.terms.iter().next().unwrap();
                m.is_empty().then(|| c.clone())
            }
            _ => None,
        }
    }

    pub fn is_one(&self) -> bool {
        self.as_constant().is_some_and(|c| c.is_one())
    }

    /// The only term, if there is exactly one.
    pub fn single(&self) -> Option<(&Mono, &Q)> {
        if self.terms.len() == 1 {
            self.terms.iter().next()
        } else {
            None
        }
    }

    fn add_term(&mut self, m: Mono, c: Q) {
        if c.is_zero() {
            return;
        }
        match self.terms.entry(m) {
            Entry::Vacant(v) => {
                v.insert(c);
            }
            Entry::Occupied(mut o) => {
                let s = o.get() + &c;
                if s.is_zero() {
                    o.remove();
                } else {
                    *o.get_mut() = s;
                }
            }
        }
    }

    pub fn absorb(&mut self, other: Sym) {
        for (m, c) in other.terms {
            self.add_term(m, c);
        }
    }

    pub fn add(&self, other: &Sym) -> Sym {
        let mut r = self.clone();
        for (m, c) in &other.terms {
            r.add_term(m.clone(), c.clone());
        }
        r
    }

    pub fn scale(&self, k: &Q) -> Sym {
        if k.is_zero() {
            return Sym::zero();
        }
        Sym { terms: self.terms.iter().map(|(m, c)| (m.clone(), c * k)).collect() }
    }

    pub fn neg(&self) -> Sym {
        self.scale(&q(-1))
    }

    pub fn sub(&self, other: &Sym) -> Sym {
        self.add(&other.neg())
    }

    pub fn mul(&self, other: &Sym) -> Sym {
        let mut r = Sym::zero();
        for (m1, c1) in &self.terms {
            for (m2, c2) in &other.terms {
                let mut m = m1.clone();
                for (a, e) in m2 {
                    let s = m.get(a).cloned().unwrap_or_else(Q::zero) + e;
                    if s.is_zero() {
                        m.remove(a);
                    } else {
                        m.insert(a.clone(), s);
                    }
                }
                r.absorb(normalize_term(c1 * c2, m));
            }
        }
        r
    }

    pub fn powi(&self, n: i64) -> Result<Sym, EvalError> {
        if n < 0 {
            return self.powi(-n)?.recip();
        }
        let mut result = Sym::int(1);
        let mut base = self.clone();
        let mut e = n as u64;
        while e > 0 {
            if e & 1 == 1 {
                result = result.mul(&base);
            }
            e >>= 1;
            if e > 0 {
                base = base.mul(&base);
            }
        }
        Ok(result)
    }

    pub fn recip(&self) -> Result<Sym, EvalError> {
        if self.is_zero() {
            return Err(EvalError::DivisionByZero);
        }
        if let Some((m, c)) = self.single() {
            let inv: Mono = m.iter().map(|(a, e)| (a.clone(), -e)).collect();
            return Ok(normalize_term(c.recip(), inv));
        }
        Ok(normalize_term(Q::one(), Mono::from([(Atom::Group(Box::new(self.clone())), q(-1))])))
    }

    pub fn div(&self, other: &Sym) -> Result<Sym, EvalError> {
        Ok(self.mul(&other.recip()?))
    }

    pub fn pow(&self, exp: &Sym) -> Result<Sym, EvalError> {
        if let Some(k) = exp.as_constant() {
            return self.pow_q(&k);
        }
        if self.is_zero() {
            return Err(EvalError::Domain("0 raised to a symbolic power is undefined".into()));
        }
        if self.is_one() {
            return Ok(Sym::int(1));
        }
        Ok(Sym::atom(Atom::Pow(Box::new(self.clone()), Box::new(exp.clone()))))
    }

    pub fn pow_q(&self, k: &Q) -> Result<Sym, EvalError> {
        if k.is_zero() {
            return Ok(Sym::int(1)); // convention: 0^0 = 1
        }
        if self.is_zero() {
            return if k.is_positive() { Ok(Sym::zero()) } else { Err(EvalError::DivisionByZero) };
        }
        if let Some(c) = self.as_constant() {
            return Ok(rational_power(&c, k));
        }
        if k.is_integer() {
            if let Some(n) = k.to_integer().to_i64().filter(|n| n.abs() <= MAX_EXPAND_POWER) {
                if n > 0 || self.single().is_some() {
                    return self.powi(n);
                }
            }
        }
        if let Some((m, c)) = self.single() {
            // (c·x)^k = c^k·x^k is only safe on the principal branch when c > 0
            // and every exponent is 1 (sqrt(x^2) is |x|, not x).
            if c.is_positive() && m.values().all(One::is_one) {
                let m2: Mono = m.iter().map(|(a, e)| (a.clone(), e * k)).collect();
                return Ok(rational_power(c, k).mul(&normalize_term(Q::one(), m2)));
            }
        }
        Ok(normalize_term(Q::one(), Mono::from([(Atom::Group(Box::new(self.clone())), k.clone())])))
    }

    /// All variable names, including inside functions and groups.
    pub fn vars(&self) -> BTreeSet<String> {
        let mut out = BTreeSet::new();
        self.collect_vars(&mut out);
        out
    }

    fn collect_vars(&self, out: &mut BTreeSet<String>) {
        for m in self.terms.keys() {
            for a in m.keys() {
                match a {
                    Atom::Var(n) => {
                        out.insert(n.clone());
                    }
                    Atom::Func(_, args) => args.iter().for_each(|s| s.collect_vars(out)),
                    Atom::Group(b) => b.collect_vars(out),
                    Atom::Pow(b, x) => {
                        b.collect_vars(out);
                        x.collect_vars(out);
                    }
                    _ => {}
                }
            }
        }
    }

    pub fn has_var(&self, x: &str) -> bool {
        self.vars().contains(x)
    }

    /// Replace variable `var` with `val` everywhere.
    pub fn subst(&self, var: &str, val: &Sym) -> Result<Sym, EvalError> {
        let mut out = Sym::zero();
        for (m, c) in &self.terms {
            let mut t = Sym::constant(c.clone());
            for (a, e) in m {
                t = t.mul(&subst_atom(a, var, val)?.pow_q(e)?);
            }
            out.absorb(t);
        }
        Ok(out)
    }

    pub fn eval(&self, env: &HashMap<String, Value>) -> Result<Value, EvalError> {
        let mut total = Value::int(0);
        for (m, c) in &self.terms {
            let mut t = Value::Rational(c.clone());
            for (a, e) in m {
                t = t.mul(&atom_value(a, env)?.pow(&Value::Rational(e.clone()))?);
            }
            total = total.add(&t);
        }
        Ok(total)
    }
}

fn subst_atom(a: &Atom, var: &str, val: &Sym) -> Result<Sym, EvalError> {
    Ok(match a {
        Atom::Var(n) if n == var => val.clone(),
        Atom::Func(name, args) => {
            func(name, args.iter().map(|s| s.subst(var, val)).collect::<Result<_, _>>()?)?
        }
        Atom::Group(b) => b.subst(var, val)?,
        Atom::Pow(b, x) => b.subst(var, val)?.pow(&x.subst(var, val)?)?,
        other => Sym::atom(other.clone()),
    })
}

fn atom_value(a: &Atom, env: &HashMap<String, Value>) -> Result<Value, EvalError> {
    Ok(match a {
        Atom::Var(n) => env.get(n).cloned().ok_or_else(|| EvalError::UnknownName {
            name: n.clone(),
            kind: NameKind::Variable,
            suggestion: None,
        })?,
        Atom::Const(Konst::Pi) => Value::Real(PI),
        Atom::Const(Konst::E) => Value::Real(E),
        Atom::I => Value::ExactComplex(Complex::new(Q::zero(), Q::one())),
        Atom::Num(n) => Value::Rational(Q::from_integer(n.clone())),
        Atom::Func(name, args) => {
            let vals = args.iter().map(|s| s.eval(env)).collect::<Result<Vec<_>, _>>()?;
            if name == "factorial" {
                vals[0].factorial()?
            } else {
                crate::eval::call(name, &vals)?
            }
        }
        Atom::Group(b) => b.eval(env)?,
        Atom::Pow(b, x) => b.eval(env)?.pow(&x.eval(env)?)?,
    })
}

/// c^k exactly, as a rational times radicals (and `i` for negative bases).
pub fn rational_power(c: &Q, k: &Q) -> Sym {
    if c.is_zero() {
        return Sym::zero();
    }
    if k.is_integer() {
        if let Some(n) = k.to_integer().to_i64().filter(|n| n.abs() <= 100_000) {
            return Sym::constant(rat_ipow(c, n));
        }
    }
    let sign = if c.is_negative() {
        let (p, d) = (k.numer(), k.denom());
        if d.is_odd() {
            // Real odd root of a negative number, matching the calculator.
            Sym::int(if p.is_odd() { -1 } else { 1 })
        } else if *d == BigInt::from(2) {
            normalize_term(Q::one(), Mono::from([(Atom::I, Q::from_integer(p.clone()))]))
        } else {
            normalize_term(Q::one(), Mono::from([(Atom::Group(Box::new(Sym::int(-1))), k.clone())]))
        }
    } else {
        Sym::int(1)
    };
    let a = c.abs();
    let num = normalize_term(Q::one(), Mono::from([(Atom::Num(a.numer().clone()), k.clone())]));
    let den = normalize_term(Q::one(), Mono::from([(Atom::Num(a.denom().clone()), -k)]));
    sign.mul(&num).mul(&den)
}

fn small_factor(n: &BigInt) -> Vec<(BigInt, u32)> {
    let mut n = n.clone();
    let mut out = Vec::new();
    let mut p = BigInt::from(2);
    let limit = BigInt::from(1_000_000);
    while &p * &p <= n && p <= limit {
        let mut k = 0;
        while (&n % &p).is_zero() {
            n /= &p;
            k += 1;
        }
        if k > 0 {
            out.push((p.clone(), k));
        }
        p += 1u32;
    }
    if n > BigInt::one() {
        out.push((n, 1));
    }
    out
}

/// n^frac = s · t^e with t as small as possible (0 < frac < 1).
fn extract_radical(n: &BigInt, frac: &Q) -> (BigInt, BigInt, Q) {
    let (Some(a), Some(b)) = (frac.numer().to_u32(), frac.denom().to_u32()) else {
        return (BigInt::one(), n.clone(), frac.clone());
    };
    let factors = small_factor(n);
    let mut s = BigInt::one();
    let mut rems = Vec::new();
    for (p, k) in &factors {
        let total = k * a;
        s *= num_traits::pow(p.clone(), (total / b) as usize);
        rems.push((p, total % b));
    }
    let g = rems.iter().fold(b, |g, (_, r)| g.gcd(r));
    let mut t = BigInt::one();
    for (p, r) in rems {
        t *= num_traits::pow(p.clone(), (r / g) as usize);
    }
    (s, t, Q::new(1.into(), (b / g).into()))
}

/// Move integer parts of radical exponents into the coefficient and merge
/// radicals that become equal, until stable.
fn normalize_nums(c: &mut Q, nums: BTreeMap<BigInt, Q>) -> BTreeMap<BigInt, Q> {
    let mut cur = nums;
    loop {
        let mut next: BTreeMap<BigInt, Q> = BTreeMap::new();
        for (n, e) in &cur {
            if n.is_one() || e.is_zero() {
                continue;
            }
            let whole = e.floor();
            let Some(w) = whole.to_integer().to_i64().filter(|w| w.abs() <= 100_000) else {
                next.insert(n.clone(), e.clone());
                continue;
            };
            *c *= rat_ipow(&Q::from_integer(n.clone()), w);
            let frac = e - &whole;
            if !frac.is_zero() {
                let (s, t, ex) = extract_radical(n, &frac);
                *c *= Q::from_integer(s);
                if !t.is_one() {
                    *next.entry(t).or_insert_with(Q::zero) += ex;
                }
            }
        }
        if next == cur {
            return next;
        }
        cur = next;
    }
}

fn normalize_term(c: Q, m: Mono) -> Sym {
    if c.is_zero() {
        return Sym::zero();
    }
    if m.keys().all(|a| matches!(a, Atom::Var(_) | Atom::Const(_) | Atom::Func(..))) {
        let mut s = Sym::zero();
        s.terms.insert(m.into_iter().filter(|(_, e)| !e.is_zero()).collect(), c);
        return s;
    }
    let mut c = c;
    let mut out = Mono::new();
    let mut nums = BTreeMap::new();
    let mut pows: BTreeMap<Sym, Sym> = BTreeMap::new();
    let mut extra: Vec<Sym> = Vec::new();
    for (a, e) in m {
        if e.is_zero() {
            continue;
        }
        match a {
            Atom::Num(n) => *nums.entry(n).or_insert_with(Q::zero) += e,
            Atom::I if e.is_integer() => {
                let k = e.to_integer().mod_floor(&BigInt::from(4)).to_u32().unwrap_or(0);
                if k >= 2 {
                    c = -c;
                }
                if k % 2 == 1 {
                    out.insert(Atom::I, Q::one());
                }
            }
            Atom::Group(base) => {
                let n = if e.is_integer() { e.to_integer().to_i64() } else { None }
                    .filter(|n| n.abs() <= MAX_EXPAND_POWER);
                match n {
                    Some(n) if n > 0 || base.single().is_some() => {
                        extra.push(base.powi(n).expect("group bases are nonzero"))
                    }
                    _ => {
                        out.insert(Atom::Group(base), e);
                    }
                }
            }
            Atom::Pow(b, x) => {
                let entry = pows.entry(*b).or_insert_with(Sym::zero);
                *entry = entry.add(&x.scale(&e));
            }
            other => {
                out.insert(other, e);
            }
        }
    }
    for (n, e) in normalize_nums(&mut c, nums) {
        out.insert(Atom::Num(n), e);
    }
    for (b, x) in pows {
        match x.as_constant() {
            Some(k) if k.is_zero() => {}
            Some(k) => extra.push(b.pow_q(&k).expect("power bases are nonzero")),
            None => {
                out.insert(Atom::Pow(Box::new(b), Box::new(x)), Q::one());
            }
        }
    }
    let mut s = Sym::zero();
    s.terms.insert(out, c);
    for x in extra {
        s = s.mul(&x);
    }
    s
}

fn expected_arity(name: &str) -> &'static [usize] {
    match name {
        "log" => &[1, 2],
        "root" => &[2],
        _ => &[1],
    }
}

/// Build a function application, simplifying exact special values.
pub fn func(name: &str, args: Vec<Sym>) -> Result<Sym, EvalError> {
    if is_command(name) {
        return Err(EvalError::Domain(format!("{name}() must be used on its own")));
    }
    if name != "factorial" && !FUNCTIONS.contains(&name) {
        return Err(EvalError::UnknownName {
            name: name.to_string(),
            kind: NameKind::Function,
            suggestion: suggest(name, FUNCTIONS),
        });
    }
    let allowed = expected_arity(name);
    if !allowed.contains(&args.len()) {
        let expected = allowed.iter().map(ToString::to_string).collect::<Vec<_>>().join(" or ");
        return Err(EvalError::Arity { name: name.to_string(), expected, got: args.len() });
    }
    match name {
        "sqrt" => return args[0].pow_q(&half()),
        "cbrt" => return args[0].pow_q(&Q::new(1.into(), 3.into())),
        "root" => {
            let n = args[1].as_constant().filter(|n| !n.is_zero()).ok_or_else(|| {
                EvalError::Domain("root index must be a nonzero number".into())
            })?;
            return args[0].pow_q(&n.recip());
        }
        "exp" => return Sym::atom(Atom::Const(Konst::E)).pow(&args[0]),
        "ln" => {
            if let Some(r) = ln_simplify(&args[0]) {
                return Ok(r);
            }
        }
        "log" if args.len() == 2 => {
            if let (Some(x), Some(b)) = (args[0].as_constant(), args[1].as_constant()) {
                if let Value::Rational(r) =
                    crate::eval::call("log", &[Value::Rational(x), Value::Rational(b)])?
                {
                    return Ok(Sym::constant(r));
                }
            }
            return func("ln", vec![args[0].clone()])?.div(&func("ln", vec![args[1].clone()])?);
        }
        _ => {}
    }
    // Exact numeric evaluation when every argument is a rational constant.
    if let Some(vals) = args.iter().map(|a| a.as_constant().map(Value::Rational)).collect::<Option<Vec<_>>>() {
        let v = if name == "factorial" { vals[0].factorial()? } else { crate::eval::call(name, &vals)? };
        match v {
            Value::Rational(r) => return Ok(Sym::constant(r)),
            Value::ExactComplex(c) => {
                return Ok(Sym::constant(c.re).add(&Sym::constant(c.im).mul(&Sym::atom(Atom::I))))
            }
            _ => {}
        }
    }
    if let Some(v) = trig_special(name, &args[0]) {
        return v;
    }
    Ok(Sym::atom(Atom::Func(name.to_string(), args)))
}

fn ln_simplify(x: &Sym) -> Option<Sym> {
    if x.is_one() {
        return Some(Sym::zero());
    }
    let (m, c) = x.single()?;
    if !c.is_one() || m.len() != 1 {
        return None;
    }
    let (a, e) = m.iter().next()?;
    match a {
        Atom::Const(Konst::E) => Some(Sym::constant(e.clone())),
        Atom::Pow(b, y) if e.is_one() && **b == Sym::atom(Atom::Const(Konst::E)) => Some((**y).clone()),
        _ => None,
    }
}

/// Exact sin/cos/tan at multiples of 30° and 45° (arguments k·pi).
fn trig_special(name: &str, arg: &Sym) -> Option<Result<Sym, EvalError>> {
    if !matches!(name, "sin" | "cos" | "tan") {
        return None;
    }
    let (m, k) = arg.single()?;
    if m.len() != 1 || m.get(&Atom::Const(Konst::Pi)) != Some(&Q::one()) {
        return None;
    }
    let deg = k * q(180);
    if !deg.is_integer() {
        return None;
    }
    let d = deg.to_integer().mod_floor(&BigInt::from(360)).to_i64()?;
    let (s, c) = (sin_deg(d)?, sin_deg(d + 90)?);
    Some(match name {
        "sin" => Ok(s),
        "cos" => Ok(c),
        _ if c.is_zero() => Err(EvalError::Domain("tan is undefined at this angle".into())),
        _ => s.div(&c),
    })
}

fn sin_deg(d: i64) -> Option<Sym> {
    let d = d.rem_euclid(360);
    if d >= 180 {
        return sin_deg(d - 180).map(|s| s.neg());
    }
    Some(match d {
        0 => Sym::zero(),
        30 | 150 => Sym::constant(half()),
        45 | 135 => rational_power(&q(2), &half()).scale(&half()),
        60 | 120 => rational_power(&q(3), &half()).scale(&half()),
        90 => Sym::int(1),
        _ => return None,
    })
}

pub fn from_expr(e: &Expr) -> Result<Sym, EvalError> {
    Ok(match e {
        Expr::Num(Value::Rational(r)) => Sym::constant(r.clone()),
        Expr::Num(Value::ExactComplex(c)) => {
            Sym::constant(c.re.clone()).add(&Sym::constant(c.im.clone()).mul(&Sym::atom(Atom::I)))
        }
        Expr::Num(_) => {
            return Err(EvalError::Domain("approximate numbers can't be used symbolically".into()))
        }
        Expr::Var(n) => match n.as_str() {
            "pi" => Sym::atom(Atom::Const(Konst::Pi)),
            "e" => Sym::atom(Atom::Const(Konst::E)),
            "i" => Sym::atom(Atom::I),
            "tau" => Sym::atom(Atom::Const(Konst::Pi)).scale(&q(2)),
            "phi" => Sym::int(1).add(&rational_power(&q(5), &half())).scale(&half()),
            _ => Sym::var(n),
        },
        Expr::Add(a, b) => from_expr(a)?.add(&from_expr(b)?),
        Expr::Sub(a, b) => from_expr(a)?.sub(&from_expr(b)?),
        Expr::Mul(a, b) => from_expr(a)?.mul(&from_expr(b)?),
        Expr::Div(a, b) => from_expr(a)?.div(&from_expr(b)?)?,
        Expr::Neg(a) => from_expr(a)?.neg(),
        Expr::Pow(a, b) => from_expr(a)?.pow(&from_expr(b)?)?,
        Expr::Factorial(a) => func("factorial", vec![from_expr(a)?])?,
        Expr::Call(name, args) => func(name, args.iter().map(from_expr).collect::<Result<_, _>>()?)?,
    })
}
