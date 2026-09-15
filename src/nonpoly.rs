//! Numerical fallback for square systems that aren't polynomial (sin, ln,
//! fractional powers, ...). Damped Newton from many deterministic starting
//! points, run in parallel; only real solutions are reported. Unlike
//! homotopy continuation this can't promise it found every solution.

use std::collections::HashMap;

use rayon::prelude::*;

use crate::sym::Sym;
use crate::value::Value;

const STARTS: usize = 4096;
const RANGE: f64 = 10.0;

fn eval(polys: &[Sym], vars: &[String], x: &[f64]) -> Option<Vec<f64>> {
    let env: HashMap<String, Value> = vars.iter().cloned().zip(x.iter().map(|&v| Value::Real(v))).collect();
    polys
        .iter()
        .map(|p| {
            let c = p.eval(&env).ok()?.to_c64();
            (c.is_finite() && c.im.abs() <= 1e-9 * (1.0 + c.re.abs())).then_some(c.re)
        })
        .collect()
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|a| a * a).sum::<f64>().sqrt()
}

fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        let p = (col..n).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
        if a[p][col].abs() < 1e-300 {
            return None;
        }
        a.swap(col, p);
        b.swap(col, p);
        for r in col + 1..n {
            let f = a[r][col] / a[col][col];
            for k in col..n {
                a[r][k] -= f * a[col][k];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        let s: f64 = (r + 1..n).map(|k| a[r][k] * x[k]).sum();
        x[r] = (b[r] - s) / a[r][r];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

fn newton(polys: &[Sym], vars: &[String], mut x: Vec<f64>) -> Option<Vec<f64>> {
    let n = x.len();
    let mut f = eval(polys, vars, &x)?;
    for _ in 0..100 {
        let mut jac = vec![vec![0.0; n]; n];
        for j in 0..n {
            let h = 1e-7 * (1.0 + x[j].abs());
            let (mut xp, mut xm) = (x.clone(), x.clone());
            xp[j] += h;
            xm[j] -= h;
            let (fp, fm) = (eval(polys, vars, &xp)?, eval(polys, vars, &xm)?);
            for i in 0..n {
                jac[i][j] = (fp[i] - fm[i]) / (2.0 * h);
            }
        }
        let dx = solve(jac, f.iter().map(|v| -v).collect())?;
        // Backtrack until the residual shrinks (or the domain is respected).
        let mut step = 1.0;
        let current = norm(&f);
        let (nx, nf) = loop {
            let cand: Vec<f64> = x.iter().zip(&dx).map(|(a, d)| a + d * step).collect();
            if let Some(cf) = eval(polys, vars, &cand) {
                if norm(&cf) < current || step < 1e-3 {
                    break (cand, cf);
                }
            }
            step *= 0.5;
            if step < 1e-4 {
                return None;
            }
        };
        let moved = norm(&dx) * step;
        x = nx;
        f = nf;
        if norm(&x) > 1e6 {
            return None;
        }
        if moved <= 1e-14 * (1.0 + norm(&x)) {
            break;
        }
    }
    (norm(&f) < 1e-9).then_some(x)
}

/// Deterministic pseudo-random starting point number `k`.
fn start(k: usize, n: usize) -> Vec<f64> {
    let mut s = (k as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ 0xD1B5_4A32_D192_ED03;
    (0..n)
        .map(|_| {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            let u = (s >> 11) as f64 / (1u64 << 53) as f64;
            // Concentrate starts near the origin but still reach ±RANGE.
            let v = 2.0 * u - 1.0;
            v * v.abs() * RANGE
        })
        .collect()
}

/// Real solutions found by multi-start Newton, sorted and deduplicated.
pub fn solve_real(polys: &[Sym], vars: &[String]) -> Vec<Vec<f64>> {
    let n = vars.len();
    let found: Vec<Vec<f64>> = (0..STARTS).into_par_iter().filter_map(|k| newton(polys, vars, start(k, n))).collect();
    let mut out: Vec<Vec<f64>> = Vec::new();
    // Newton can jump far away on periodic functions; only keep the window
    // that was actually searched, where the list is reasonably complete.
    for s in found.into_iter().filter(|s| s.iter().all(|v| v.abs() <= RANGE)) {
        let s: Vec<f64> = s.iter().map(|&v| if v.abs() < 1e-12 { 0.0 } else { v }).collect();
        if !out.iter().any(|o| o.iter().zip(&s).all(|(a, b)| (a - b).abs() < 1e-7 * (1.0 + a.abs()))) {
            out.push(s);
        }
    }
    out.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    out
}
