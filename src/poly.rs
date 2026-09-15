//! Polynomial views of `Sym`: monomial ordering, exact division, and
//! conversion to dense univariate coefficient vectors.

use std::cmp::Ordering;
use std::collections::BTreeSet;

use num_rational::BigRational as Q;
use num_traits::{Signed, ToPrimitive, Zero};

use crate::sym::{Atom, Mono, Sym};

/// Lexicographic monomial order over atoms (variables alphabetically first).
pub fn lex_cmp(a: &Mono, b: &Mono) -> Ordering {
    let keys: BTreeSet<&Atom> = a.keys().chain(b.keys()).collect();
    for k in keys {
        let ea = a.get(k).cloned().unwrap_or_else(Q::zero);
        let eb = b.get(k).cloned().unwrap_or_else(Q::zero);
        match ea.cmp(&eb) {
            Ordering::Equal => continue,
            o => return o,
        }
    }
    Ordering::Equal
}

pub fn leading(s: &Sym) -> Option<(Mono, Q)> {
    s.terms.iter().max_by(|a, b| lex_cmp(a.0, b.0)).map(|(m, c)| (m.clone(), c.clone()))
}

/// a / b if every exponent stays non-negative.
pub fn mono_div(a: &Mono, b: &Mono) -> Option<Mono> {
    let mut out = a.clone();
    for (k, e) in b {
        let r = out.get(k).cloned().unwrap_or_else(Q::zero) - e;
        if r.is_negative() {
            return None;
        }
        if r.is_zero() {
            out.remove(k);
        } else {
            out.insert(k.clone(), r);
        }
    }
    Some(out)
}

pub fn mono_lcm(a: &Mono, b: &Mono) -> Mono {
    let mut out = a.clone();
    for (k, e) in b {
        let cur = out.entry(k.clone()).or_insert_with(Q::zero);
        if *e > *cur {
            *cur = e.clone();
        }
    }
    out
}

pub fn raw_term(m: Mono, c: Q) -> Sym {
    let mut s = Sym::zero();
    if !c.is_zero() {
        s.terms.insert(m, c);
    }
    s
}

/// Every exponent is a non-negative integer.
pub fn is_polynomial(s: &Sym) -> bool {
    s.terms.keys().all(|m| m.values().all(|e| e.is_integer() && !e.is_negative()))
}

/// Exact quotient p / d, if d divides p as polynomials (atoms act as variables).
pub fn divide(p: &Sym, d: &Sym) -> Option<Sym> {
    if d.is_zero() || !is_polynomial(p) || !is_polynomial(d) {
        return None;
    }
    let (dm, dc) = leading(d)?;
    let mut rem = p.clone();
    let mut quot = Sym::zero();
    for _ in 0..10_000 {
        let Some((m, c)) = leading(&rem) else { return Some(quot) };
        let t = raw_term(mono_div(&m, &dm)?, &c / &dc);
        rem = rem.sub(&t.mul(d));
        quot.absorb(t);
    }
    None
}

pub fn atom_has_var(a: &Atom, x: &str) -> bool {
    match a {
        Atom::Var(n) => n == x,
        Atom::Func(_, args) => args.iter().any(|s| s.has_var(x)),
        Atom::Group(b) => b.has_var(x),
        Atom::Pow(b, e) => b.has_var(x) || e.has_var(x),
        _ => false,
    }
}

/// Coefficients of x^0, x^1, …, if `s` is a polynomial in `x` (other atoms
/// may appear in the coefficients as long as they don't involve x).
pub fn coeffs_in(s: &Sym, x: &str) -> Option<Vec<Sym>> {
    let xa = Atom::Var(x.to_string());
    let mut out: Vec<Sym> = Vec::new();
    for (m, c) in &s.terms {
        let e = match m.get(&xa) {
            None => 0,
            Some(e) if e.is_integer() && !e.is_negative() => e.to_integer().to_usize()?,
            Some(_) => return None,
        };
        let mut rest = m.clone();
        rest.remove(&xa);
        if rest.keys().any(|a| atom_has_var(a, x)) {
            return None;
        }
        if out.len() <= e {
            out.resize(e + 1, Sym::zero());
        }
        out[e].absorb(raw_term(rest, c.clone()));
    }
    while out.last().is_some_and(Sym::is_zero) {
        out.pop();
    }
    Some(out)
}

pub fn degree_in(s: &Sym, x: &str) -> Option<usize> {
    coeffs_in(s, x).map(|c| c.len().saturating_sub(1))
}

/// Dense rational coefficients (low → high), if `s` is a polynomial in x
/// with rational coefficients.
pub fn uni_q(s: &Sym, x: &str) -> Option<Vec<Q>> {
    coeffs_in(s, x)?.iter().map(Sym::as_constant).collect()
}

pub fn from_uni(coeffs: &[Q], x: &str) -> Sym {
    let mut out = Sym::zero();
    for (k, c) in coeffs.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let m = if k == 0 {
            Mono::new()
        } else {
            Mono::from([(Atom::Var(x.to_string()), Q::from_integer((k as i64).into()))])
        };
        out.absorb(raw_term(m, c.clone()));
    }
    out
}
