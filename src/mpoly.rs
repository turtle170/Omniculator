//! Fast sparse multivariate polynomials over Q for Gröbner bases.
//!
//! Terms are (exponent vector, coefficient) sorted in descending lex order,
//! so the leading term is `terms[0]` and subtraction is a linear merge.

use std::cmp::Ordering;
use std::time::Instant;

use num_rational::BigRational as Q;
use num_traits::{One, Signed, ToPrimitive, Zero};

use crate::poly::raw_term;
use crate::sym::{Atom, Mono, Sym};

pub type Exps = Vec<u32>;

/// Give up (and let the caller fall back) past this many terms or basis size.
const MAX_TERMS: usize = 60_000;
const MAX_BASIS: usize = 600;

#[derive(Clone, Debug, PartialEq)]
pub struct MPoly {
    pub terms: Vec<(Exps, Q)>,
}

pub fn divides(a: &[u32], b: &[u32]) -> bool {
    a.iter().zip(b).all(|(x, y)| x <= y)
}

fn lcm(a: &[u32], b: &[u32]) -> Exps {
    a.iter().zip(b).map(|(x, y)| *x.max(y)).collect()
}

fn diff(b: &[u32], a: &[u32]) -> Exps {
    b.iter().zip(a).map(|(x, y)| x - y).collect()
}

fn add(a: &[u32], b: &[u32]) -> Exps {
    a.iter().zip(b).map(|(x, y)| x + y).collect()
}

fn tdeg(a: &[u32]) -> u32 {
    a.iter().sum()
}

fn coprime(a: &[u32], b: &[u32]) -> bool {
    a.iter().zip(b).all(|(x, y)| *x == 0 || *y == 0)
}

impl MPoly {
    /// None unless `s` is a polynomial in exactly these variables.
    pub fn from_sym(s: &Sym, vars: &[String]) -> Option<MPoly> {
        let mut terms = Vec::with_capacity(s.terms.len());
        for (m, c) in &s.terms {
            let mut e = vec![0u32; vars.len()];
            for (a, x) in m {
                let Atom::Var(n) = a else { return None };
                let i = vars.iter().position(|v| v == n)?;
                if !x.is_integer() || x.is_negative() {
                    return None;
                }
                e[i] = x.to_integer().to_u32()?;
            }
            terms.push((e, c.clone()));
        }
        terms.sort_by(|a, b| b.0.cmp(&a.0));
        Some(MPoly { terms })
    }

    pub fn to_sym(&self, vars: &[String]) -> Sym {
        let mut out = Sym::zero();
        for (e, c) in &self.terms {
            let m: Mono = e
                .iter()
                .enumerate()
                .filter(|(_, k)| **k > 0)
                .map(|(i, k)| (Atom::Var(vars[i].clone()), Q::from_integer((*k as i64).into())))
                .collect();
            out.absorb(raw_term(m, c.clone()));
        }
        out
    }

    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn lm(&self) -> &[u32] {
        &self.terms[0].0
    }

    pub fn total_degree(&self) -> u32 {
        self.terms.iter().map(|(e, _)| tdeg(e)).max().unwrap_or(0)
    }

    fn monic(mut self) -> MPoly {
        let l = self.terms[0].1.clone();
        if !l.is_one() {
            for t in &mut self.terms {
                t.1 = &t.1 / &l;
            }
        }
        self
    }

    fn shifted(&self, s: &[u32]) -> MPoly {
        MPoly { terms: self.terms.iter().map(|(e, c)| (add(e, s), c.clone())).collect() }
    }

    /// self − c·x^s·g
    fn sub_mul(&self, c: &Q, s: &[u32], g: &MPoly) -> MPoly {
        let mut out = Vec::with_capacity(self.terms.len() + g.terms.len());
        let mut a = self.terms.iter().peekable();
        let mut b = g.terms.iter().map(|(e, k)| (add(e, s), -(c * k))).peekable();
        loop {
            let ord = match (a.peek(), b.peek()) {
                (Some(x), Some(y)) => x.0.cmp(&y.0),
                (Some(_), None) => Ordering::Greater,
                (None, Some(_)) => Ordering::Less,
                (None, None) => break,
            };
            match ord {
                Ordering::Greater => out.push(a.next().unwrap().clone()),
                Ordering::Less => out.push(b.next().unwrap()),
                Ordering::Equal => {
                    let x = a.next().unwrap();
                    let (e, y) = b.next().unwrap();
                    let sum = &x.1 + y;
                    if !sum.is_zero() {
                        out.push((e, sum));
                    }
                }
            }
        }
        MPoly { terms: out }
    }

    pub fn eval_q(&self, x: &[Q]) -> Q {
        self.terms
            .iter()
            .map(|(e, c)| {
                e.iter().enumerate().fold(c.clone(), |acc, (i, &k)| acc * num_traits::pow(x[i].clone(), k as usize))
            })
            .fold(Q::zero(), |a, b| a + b)
    }
}

/// Full reduction of p modulo g. Err(()) past the deadline or size limit.
fn reduce(p: &MPoly, g: &[MPoly], deadline: Instant) -> Result<MPoly, ()> {
    let mut p = p.clone();
    let mut i = 0;
    let mut steps = 0u32;
    while i < p.terms.len() {
        let e = &p.terms[i].0;
        match g.iter().find(|gk| divides(gk.lm(), e)) {
            Some(gk) => {
                let shift = diff(e, gk.lm());
                let coef = &p.terms[i].1 / &gk.terms[0].1;
                p = p.sub_mul(&coef, &shift, gk);
                steps += 1;
                if steps % 32 == 0 && Instant::now() > deadline {
                    return Err(());
                }
                if p.terms.len() > MAX_TERMS {
                    return Err(());
                }
            }
            None => i += 1,
        }
    }
    Ok(p)
}

struct Pair {
    i: usize,
    j: usize,
    lcm: Exps,
    sugar: u32,
}

fn update(g: &mut Vec<MPoly>, sugar: &mut Vec<u32>, pairs: &mut Vec<Pair>, h: MPoly, s: u32) {
    let n = g.len();
    let hl = h.lm().to_vec();
    // Gebauer–Möller: drop pairs made redundant by h.
    pairs.retain(|p| {
        !(divides(&hl, &p.lcm) && lcm(g[p.i].lm(), &hl) != p.lcm && lcm(g[p.j].lm(), &hl) != p.lcm)
    });
    for k in 0..n {
        let gl = g[k].lm();
        if coprime(gl, &hl) {
            continue; // Buchberger's first criterion
        }
        let l = lcm(gl, &hl);
        let ps = (sugar[k] + tdeg(&diff(&l, gl))).max(s + tdeg(&diff(&l, &hl)));
        pairs.push(Pair { i: k, j: n, lcm: l, sugar: ps });
    }
    g.push(h);
    sugar.push(s);
}

/// Reduced lex Gröbner basis (Buchberger with the sugar strategy).
/// Err(()) if it doesn't finish before `deadline` or grows too large.
pub fn groebner(input: &[MPoly], deadline: Instant) -> Result<Vec<MPoly>, ()> {
    let mut g: Vec<MPoly> = Vec::new();
    let mut sugar: Vec<u32> = Vec::new();
    let mut pairs: Vec<Pair> = Vec::new();
    for p in input.iter().filter(|p| !p.is_zero()) {
        let r = reduce(p, &g, deadline)?;
        if !r.is_zero() {
            let s = r.total_degree();
            update(&mut g, &mut sugar, &mut pairs, r.monic(), s);
        }
    }
    while !pairs.is_empty() {
        if Instant::now() > deadline || g.len() > MAX_BASIS {
            return Err(());
        }
        let k = (0..pairs.len()).min_by_key(|&k| (pairs[k].sugar, tdeg(&pairs[k].lcm))).unwrap();
        let pr = pairs.swap_remove(k);
        let (fi, fj) = (&g[pr.i], &g[pr.j]);
        let sp = fi.shifted(&diff(&pr.lcm, fi.lm())).sub_mul(&Q::one(), &diff(&pr.lcm, fj.lm()), fj);
        let r = reduce(&sp, &g, deadline)?;
        if !r.is_zero() {
            update(&mut g, &mut sugar, &mut pairs, r.monic(), pr.sugar);
        }
    }
    // Minimal basis, then reduce every tail by the others.
    let keep: Vec<MPoly> = g
        .iter()
        .enumerate()
        .filter(|(i, p)| {
            !g.iter().enumerate().any(|(j, q)| j != *i && divides(q.lm(), p.lm()) && (q.lm() != p.lm() || j < *i))
        })
        .map(|(_, p)| p.clone())
        .collect();
    let mut out = Vec::with_capacity(keep.len());
    for i in 0..keep.len() {
        let others: Vec<MPoly> = keep.iter().enumerate().filter(|(j, _)| *j != i).map(|(_, p)| p.clone()).collect();
        let tail = MPoly { terms: keep[i].terms[1..].to_vec() };
        let mut terms = vec![keep[i].terms[0].clone()];
        terms.extend(reduce(&tail, &others, deadline)?.terms);
        out.push(MPoly { terms });
    }
    out.sort_by(|a, b| b.lm().cmp(a.lm()));
    Ok(out)
}
