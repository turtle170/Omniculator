//! Dense univariate polynomials over the rationals (coefficients low → high):
//! arithmetic, rational roots, factoring over Q, and numeric roots.

use num_bigint::BigInt;
use num_complex::Complex64;
use num_integer::Integer;
use num_rational::BigRational as Q;
use num_traits::{One, Signed, ToPrimitive, Zero};

use crate::sym::{half, rational_power, Sym};

pub type UPoly = Vec<Q>;

/// Kronecker search gives up beyond this many candidate factors.
const KRONECKER_BUDGET: usize = 200_000;
const ROOT_CANDIDATE_BUDGET: usize = 50_000;

pub fn trim(mut p: UPoly) -> UPoly {
    while p.last().is_some_and(Zero::is_zero) {
        p.pop();
    }
    p
}

pub fn degree(p: &[Q]) -> usize {
    p.len().saturating_sub(1)
}

pub fn eval_q(p: &[Q], x: &Q) -> Q {
    p.iter().rev().fold(Q::zero(), |acc, c| acc * x + c)
}

pub fn eval_c(p: &[Complex64], z: Complex64) -> Complex64 {
    p.iter().rev().fold(Complex64::new(0.0, 0.0), |acc, c| acc * z + c)
}

pub fn derivative(p: &[Q]) -> UPoly {
    p.iter().enumerate().skip(1).map(|(k, c)| c * Q::from_integer(k.into())).collect()
}

pub fn mul(a: &[Q], b: &[Q]) -> UPoly {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut out = vec![Q::zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    trim(out)
}

/// (quotient, remainder); `b` must be nonzero.
pub fn divrem(a: &[Q], b: &[Q]) -> (UPoly, UPoly) {
    let b = trim(b.to_vec());
    let mut r = trim(a.to_vec());
    if r.len() < b.len() {
        return (vec![], r);
    }
    let lead = b.last().expect("nonzero divisor").clone();
    let mut quot = vec![Q::zero(); r.len() - b.len() + 1];
    while r.len() >= b.len() && !r.is_empty() {
        let shift = r.len() - b.len();
        let f = r.last().unwrap() / &lead;
        for (k, c) in b.iter().enumerate() {
            r[shift + k] -= &f * c;
        }
        quot[shift] = f;
        r.pop();
        r = trim(r);
    }
    (trim(quot), r)
}

/// Monic gcd.
pub fn gcd(a: &[Q], b: &[Q]) -> UPoly {
    let (mut a, mut b) = (trim(a.to_vec()), trim(b.to_vec()));
    while !b.is_empty() {
        let (_, r) = divrem(&a, &b);
        a = b;
        b = r;
    }
    monic(&a)
}

pub fn monic(p: &[Q]) -> UPoly {
    match p.last() {
        Some(l) if !l.is_zero() => p.iter().map(|c| c / l).collect(),
        _ => p.to_vec(),
    }
}

/// Scale to integer coefficients with gcd 1 and positive leading coefficient.
pub fn primitive(p: &[Q]) -> UPoly {
    let lcm = p.iter().fold(BigInt::one(), |l, c| l.lcm(c.denom()));
    let ints: Vec<BigInt> = p.iter().map(|c| (c * Q::from_integer(lcm.clone())).to_integer()).collect();
    let g = ints.iter().fold(BigInt::zero(), |g, c| g.gcd(c));
    if g.is_zero() {
        return p.to_vec();
    }
    let sign = if ints.last().is_some_and(Signed::is_negative) { -BigInt::one() } else { BigInt::one() };
    ints.into_iter().map(|c| Q::from_integer(c / &g * &sign)).collect()
}

/// Positive divisors, or None if |n| can't be fully factored cheaply.
fn divisors(n: &BigInt, limit: usize) -> Option<Vec<BigInt>> {
    let mut n = n.abs();
    if n.is_zero() {
        return None;
    }
    let mut factors: Vec<(BigInt, u32)> = Vec::new();
    let mut p = BigInt::from(2);
    while &p * &p <= n {
        if p > BigInt::from(2_000_000) {
            return None; // leftover might be composite; can't enumerate honestly
        }
        let mut k = 0;
        while (&n % &p).is_zero() {
            n /= &p;
            k += 1;
        }
        if k > 0 {
            factors.push((p.clone(), k));
        }
        p += 1u32;
    }
    if n > BigInt::one() {
        factors.push((n, 1));
    }
    let mut divs = vec![BigInt::one()];
    for (p, k) in factors {
        let mut next = Vec::new();
        for d in &divs {
            let mut pk = BigInt::one();
            for _ in 0..=k {
                next.push(d * &pk);
                pk *= &p;
            }
        }
        divs = next;
        if divs.len() > limit {
            return None;
        }
    }
    Some(divs)
}

/// All rational roots (distinct), or None if the candidates couldn't be enumerated.
pub fn rational_roots(p: &[Q]) -> Option<Vec<Q>> {
    let p = primitive(&trim(p.to_vec()));
    if p.len() <= 1 {
        return Some(vec![]);
    }
    let mut roots = Vec::new();
    let low = p.iter().position(|c| !c.is_zero()).unwrap_or(0);
    if low > 0 {
        roots.push(Q::zero());
    }
    let (a0, an) = (p[low].to_integer(), p.last().unwrap().to_integer());
    let (num, den) = (divisors(&a0, ROOT_CANDIDATE_BUDGET)?, divisors(&an, ROOT_CANDIDATE_BUDGET)?);
    if num.len() * den.len() > ROOT_CANDIDATE_BUDGET {
        return None;
    }
    for n in &num {
        for d in &den {
            for cand in [Q::new(n.clone(), d.clone()), Q::new(-n, d.clone())] {
                if !roots.contains(&cand) && eval_q(&p, &cand).is_zero() {
                    roots.push(cand);
                }
            }
        }
    }
    roots.sort();
    Some(roots)
}

/// Factorization over Q: `content · Π factor^multiplicity`, factors primitive
/// with positive leading coefficient. `complete` is false when a search limit
/// was hit, so some factor may still be reducible.
#[derive(Debug, Clone, PartialEq)]
pub struct Factored {
    pub content: Q,
    pub factors: Vec<(UPoly, usize)>,
    pub complete: bool,
}

fn push_factor(factors: &mut Vec<(UPoly, usize)>, f: UPoly) {
    match factors.iter_mut().find(|(g, _)| *g == f) {
        Some((_, k)) => *k += 1,
        None => factors.push((f, 1)),
    }
}

pub fn factor(p: &[Q]) -> Factored {
    let p = trim(p.to_vec());
    let mut out = Factored { content: Q::zero(), factors: vec![], complete: true };
    if p.is_empty() {
        return out;
    }
    let mut rest = p.clone();
    // Rational roots → linear factors (b·x − a), with multiplicity.
    match rational_roots(&rest) {
        Some(roots) => {
            for r in roots {
                let lin = primitive(&[-r.clone(), Q::one()]);
                loop {
                    let (quot, rem) = divrem(&rest, &lin);
                    if !rem.is_empty() {
                        break;
                    }
                    push_factor(&mut out.factors, lin.clone());
                    rest = quot;
                }
            }
        }
        None => out.complete = false,
    }
    // Higher-degree pieces via Kronecker's method.
    let mut stack = vec![rest];
    let mut leftover_const = Q::one();
    while let Some(f) = stack.pop() {
        if degree(&f) == 0 {
            leftover_const *= f.first().cloned().unwrap_or_else(Q::one);
            continue;
        }
        let prim = primitive(&f);
        leftover_const *= &f.last().unwrap().clone() / prim.last().unwrap();
        if degree(&prim) <= 3 {
            // No rational roots left, so degree 2–3 pieces are irreducible.
            push_factor(&mut out.factors, prim);
            continue;
        }
        let mut split = None;
        for d in 2..=degree(&prim) / 2 {
            match kronecker(&prim, d) {
                Ok(Some(g)) => {
                    split = Some(g);
                    break;
                }
                Ok(None) => {}
                Err(()) => {
                    out.complete = false;
                    break;
                }
            }
        }
        match split {
            Some(g) => {
                let (quot, _) = divrem(&prim, &g);
                stack.push(g);
                stack.push(quot);
            }
            None => push_factor(&mut out.factors, prim),
        }
    }
    // Whatever scalar is left over is the content.
    let prod = out.factors.iter().fold(vec![Q::one()], |acc, (f, k)| {
        (0..*k).fold(acc, |a, _| mul(&a, f))
    });
    out.content = p.last().unwrap() / prod.last().unwrap();
    let _ = leftover_const;
    out.factors.sort_by(|a, b| degree(&a.0).cmp(&degree(&b.0)).then_with(|| a.0.cmp(&b.0)));
    out
}

/// Find a factor of exact degree `d` of primitive `f` by Kronecker's method.
/// Err(()) if the search is too large.
fn kronecker(f: &[Q], d: usize) -> Result<Option<UPoly>, ()> {
    let mut points: Vec<Q> = Vec::new();
    let mut k: i64 = 0;
    while points.len() <= d {
        let a = Q::from_integer(k.into());
        if !eval_q(f, &a).is_zero() {
            points.push(a);
        }
        k = if k > 0 { -k } else { -k + 1 };
    }
    let mut choices: Vec<Vec<Q>> = Vec::new();
    let mut total: usize = 1;
    for (i, a) in points.iter().enumerate() {
        let v = eval_q(f, a).to_integer();
        let divs = divisors(&v, KRONECKER_BUDGET).ok_or(())?;
        let mut c: Vec<Q> = divs.iter().map(|x| Q::from_integer(x.clone())).collect();
        if i > 0 {
            c.extend(divs.iter().map(|x| Q::from_integer(-x)));
        }
        total = total.saturating_mul(c.len());
        if total > KRONECKER_BUDGET {
            return Err(());
        }
        choices.push(c);
    }
    // Lagrange basis polynomials for the chosen points.
    let basis: Vec<UPoly> = (0..points.len())
        .map(|i| {
            let mut l = vec![Q::one()];
            for (j, aj) in points.iter().enumerate() {
                if i != j {
                    let denom = &points[i] - aj;
                    l = mul(&l, &[-aj / &denom, Q::one() / &denom]);
                }
            }
            l
        })
        .collect();
    let mut idx = vec![0usize; points.len()];
    loop {
        let mut g = vec![Q::zero(); d + 1];
        for (i, &j) in idx.iter().enumerate() {
            for (k, c) in basis[i].iter().enumerate() {
                g[k] += &choices[i][j] * c;
            }
        }
        let g = trim(g);
        if degree(&g) == d && g.iter().all(Q::is_integer) {
            let (_, rem) = divrem(f, &g);
            if rem.is_empty() {
                return Ok(Some(primitive(&g)));
            }
        }
        // Advance the odometer.
        let mut pos = 0;
        loop {
            if pos == idx.len() {
                return Ok(None);
            }
            idx[pos] += 1;
            if idx[pos] < choices[pos].len() {
                break;
            }
            idx[pos] = 0;
            pos += 1;
        }
    }
}

/// All complex roots numerically (Aberth–Ehrlich), polished with Newton steps.
pub fn roots_numeric(p: &[Complex64]) -> Vec<Complex64> {
    let mut p = p.to_vec();
    while p.last().is_some_and(|c| c.norm() == 0.0) {
        p.pop();
    }
    let n = p.len().saturating_sub(1);
    if n == 0 {
        return vec![];
    }
    let lead = p[n];
    let p: Vec<Complex64> = p.iter().map(|c| c / lead).collect();
    let dp: Vec<Complex64> = p.iter().enumerate().skip(1).map(|(k, c)| c * k as f64).collect();
    let radius = 1.0 + p[..n].iter().map(|c| c.norm()).fold(0.0, f64::max);
    let mut z: Vec<Complex64> = (0..n)
        .map(|k| Complex64::from_polar(radius * 0.5, 2.0 * std::f64::consts::PI * k as f64 / n as f64 + 0.4))
        .collect();
    for _ in 0..1000 {
        let mut max_step: f64 = 0.0;
        for k in 0..n {
            let ratio = eval_c(&p, z[k]) / eval_c(&dp, z[k]);
            let sum: Complex64 = (0..n).filter(|&j| j != k).map(|j| 1.0 / (z[k] - z[j])).sum();
            let w = ratio / (1.0 - ratio * sum);
            if w.is_finite() {
                z[k] -= w;
                max_step = max_step.max(w.norm() / (1.0 + z[k].norm()));
            }
        }
        if max_step < 1e-15 {
            break;
        }
    }
    for r in z.iter_mut() {
        for _ in 0..3 {
            let step = eval_c(&p, *r) / eval_c(&dp, *r);
            if step.is_finite() {
                *r -= step;
            }
        }
    }
    z
}

pub fn to_c64(p: &[Q]) -> Vec<Complex64> {
    p.iter().map(|c| Complex64::new(c.to_f64().unwrap_or(f64::NAN), 0.0)).collect()
}

/// Exact roots of a·x² + b·x + c: (discriminant, (−b + √D)/2a, (−b − √D)/2a).
pub fn quadratic_roots(a: &Q, b: &Q, c: &Q) -> (Q, Sym, Sym) {
    let disc = b * b - Q::from_integer(4.into()) * a * c;
    let sqrt_d = rational_power(&disc, &half());
    let two_a = a * Q::from_integer(2.into());
    let minus_b = Sym::constant(-b);
    let r1 = minus_b.add(&sqrt_d).scale(&two_a.recip());
    let r2 = minus_b.sub(&sqrt_d).scale(&two_a.recip());
    (disc, r1, r2)
}
