//! Numerical polynomial homotopy continuation for square systems.
//!
//! Deforms the start system x_i^(d_i) = 1 (whose roots are known) into the
//! target system along H(x, t) = (1−t)·γ·g(x) + t·f(x), tracking every one of
//! the Π d_i paths with an RK4 predictor and Newton corrector. Paths run in
//! parallel. By the Bézout bound this finds every isolated solution; paths
//! that run off to infinity are dropped. γ is a fixed generic constant, so
//! results are deterministic.

use num_complex::Complex64 as C;
use num_traits::ToPrimitive;
use rayon::prelude::*;

use crate::mpoly::MPoly;

const GAMMA: C = C::new(0.832_328_4, 0.554_279_1);
const MAX_PATHS: usize = 1 << 18;

struct Term {
    c: C,
    e: Vec<u32>,
}

struct Sys {
    n: usize,
    f: Vec<Vec<Term>>,
    df: Vec<Vec<Vec<Term>>>,
    deg: Vec<u32>,
    maxdeg: usize,
}

impl Sys {
    fn new(polys: &[MPoly], n: usize) -> Sys {
        let f: Vec<Vec<Term>> = polys
            .iter()
            .map(|p| p.terms.iter().map(|(e, c)| Term { c: C::new(c.to_f64().unwrap_or(f64::NAN), 0.0), e: e.clone() }).collect())
            .collect();
        let df = f
            .iter()
            .map(|p| {
                (0..n)
                    .map(|j| {
                        p.iter()
                            .filter(|t| t.e[j] > 0)
                            .map(|t| {
                                let mut e = t.e.clone();
                                e[j] -= 1;
                                Term { c: t.c * t.e[j] as f64, e }
                            })
                            .collect()
                    })
                    .collect()
            })
            .collect();
        let deg = polys.iter().map(MPoly::total_degree).collect();
        let maxdeg = f.iter().flatten().flat_map(|t| t.e.iter()).copied().max().unwrap_or(0) as usize;
        Sys { n, f, df, deg, maxdeg }
    }

    fn powers(&self, x: &[C]) -> Vec<Vec<C>> {
        x.iter()
            .map(|&xi| {
                let mut v = Vec::with_capacity(self.maxdeg + 1);
                let mut p = C::new(1.0, 0.0);
                for _ in 0..=self.maxdeg {
                    v.push(p);
                    p *= xi;
                }
                v
            })
            .collect()
    }

    fn eval_terms(terms: &[Term], pw: &[Vec<C>]) -> C {
        terms.iter().map(|t| t.e.iter().enumerate().fold(t.c, |acc, (j, &k)| acc * pw[j][k as usize])).sum()
    }

    /// f(x) and its Jacobian.
    fn eval(&self, x: &[C]) -> (Vec<C>, Vec<Vec<C>>) {
        let pw = self.powers(x);
        let f = self.f.iter().map(|p| Self::eval_terms(p, &pw)).collect();
        let j = self.df.iter().map(|row| row.iter().map(|p| Self::eval_terms(p, &pw)).collect()).collect();
        (f, j)
    }

    /// H, ∂H/∂x, ∂H/∂t at (x, t).
    fn homotopy(&self, x: &[C], t: f64) -> (Vec<C>, Vec<Vec<C>>, Vec<C>) {
        let (f, mut hx) = self.eval(x);
        let s = 1.0 - t;
        let mut h = Vec::with_capacity(self.n);
        let mut ht = Vec::with_capacity(self.n);
        for i in 0..self.n {
            let d = self.deg[i] as i32;
            let g = x[i].powi(d) - 1.0;
            h.push(GAMMA * s * g + f[i] * t);
            ht.push(f[i] - GAMMA * g);
            for v in hx[i].iter_mut() {
                *v *= t;
            }
            hx[i][i] += GAMMA * s * x[i].powi(d - 1) * d as f64;
        }
        (h, hx, ht)
    }
}

/// Solve a·x = b by Gaussian elimination with partial pivoting.
fn solve(mut a: Vec<Vec<C>>, mut b: Vec<C>) -> Option<Vec<C>> {
    let n = b.len();
    for col in 0..n {
        let p = (col..n).max_by(|&i, &j| a[i][col].norm().total_cmp(&a[j][col].norm()))?;
        if a[p][col].norm() < 1e-300 {
            return None;
        }
        a.swap(col, p);
        b.swap(col, p);
        for r in col + 1..n {
            let f = a[r][col] / a[col][col];
            if f.norm() == 0.0 {
                continue;
            }
            for k in col..n {
                let v = a[col][k];
                a[r][k] -= f * v;
            }
            let v = b[col];
            b[r] -= f * v;
        }
    }
    let mut x = vec![C::new(0.0, 0.0); n];
    for r in (0..n).rev() {
        let s: C = (r + 1..n).map(|k| a[r][k] * x[k]).sum();
        x[r] = (b[r] - s) / a[r][r];
    }
    x.iter().all(|z| z.is_finite()).then_some(x)
}

fn norm(v: &[C]) -> f64 {
    v.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt()
}

fn axpy(x: &[C], k: &[C], h: f64) -> Vec<C> {
    x.iter().zip(k).map(|(a, b)| a + b * h).collect()
}

fn velocity(sys: &Sys, x: &[C], t: f64) -> Option<Vec<C>> {
    let (_, hx, ht) = sys.homotopy(x, t);
    solve(hx, ht.iter().map(|z| -z).collect())
}

fn newton_h(sys: &Sys, mut x: Vec<C>, t: f64) -> Option<Vec<C>> {
    for _ in 0..4 {
        let (h, hx, _) = sys.homotopy(&x, t);
        let dx = solve(hx, h.iter().map(|z| -z).collect())?;
        for (a, d) in x.iter_mut().zip(&dx) {
            *a += d;
        }
        if norm(&dx) <= 1e-9 * (1.0 + norm(&x)) {
            return Some(x);
        }
    }
    None
}

fn track(sys: &Sys, start: Vec<C>) -> Option<Vec<C>> {
    let mut x = start;
    let mut t = 0.0f64;
    let mut h = 0.01f64;
    for _ in 0..200_000 {
        if t >= 1.0 {
            return polish(sys, x);
        }
        let dt = h.min(1.0 - t);
        let predicted = (|| {
            let k1 = velocity(sys, &x, t)?;
            let k2 = velocity(sys, &axpy(&x, &k1, dt / 2.0), t + dt / 2.0)?;
            let k3 = velocity(sys, &axpy(&x, &k2, dt / 2.0), t + dt / 2.0)?;
            let k4 = velocity(sys, &axpy(&x, &k3, dt), t + dt)?;
            Some((0..x.len()).map(|i| x[i] + (k1[i] + k2[i] * 2.0 + k3[i] * 2.0 + k4[i]) * (dt / 6.0)).collect::<Vec<_>>())
        })();
        match predicted.and_then(|xp| newton_h(sys, xp, t + dt)) {
            Some(xc) => {
                x = xc;
                t += dt;
                h = (h * 1.6).min(0.05);
                if norm(&x) > 1e8 {
                    return None; // diverging: a solution at infinity
                }
            }
            None => {
                h *= 0.5;
                if h < 1e-14 {
                    return None;
                }
            }
        }
    }
    None
}

/// Newton on the target system, then a residual check.
fn polish(sys: &Sys, mut x: Vec<C>) -> Option<Vec<C>> {
    for _ in 0..80 {
        let (f, j) = sys.eval(&x);
        let Some(dx) = solve(j, f.iter().map(|z| -z).collect()) else { break };
        for (a, d) in x.iter_mut().zip(&dx) {
            *a += d;
        }
        if norm(&dx) <= 1e-15 * (1.0 + norm(&x)) {
            break;
        }
    }
    let pw = sys.powers(&x);
    for p in &sys.f {
        let value = Sys::eval_terms(p, &pw).norm();
        let scale: f64 = p.iter().map(|t| t.e.iter().enumerate().fold(t.c.norm(), |acc, (j, &k)| acc * pw[j][k as usize].norm())).sum();
        if value > 1e-8 * (1.0 + scale) {
            return None;
        }
    }
    Some(x)
}

/// All isolated complex solutions of a square system, and the number of
/// paths tracked.
pub fn solve_system(polys: &[MPoly], n: usize) -> Result<(Vec<Vec<C>>, usize), String> {
    let sys = Sys::new(polys, n);
    if polys.len() != n || sys.deg.iter().any(|&d| d == 0) {
        return Err("numerical solving needs as many equations as unknowns".into());
    }
    let paths = sys.deg.iter().try_fold(1usize, |acc, &d| acc.checked_mul(d as usize)).filter(|&p| p <= MAX_PATHS);
    let Some(paths) = paths else {
        return Err("this system has too many potential solutions to track".into());
    };
    let start = |mut k: usize| -> Vec<C> {
        sys.deg
            .iter()
            .map(|&d| {
                let j = k % d as usize;
                k /= d as usize;
                C::from_polar(1.0, 2.0 * std::f64::consts::PI * j as f64 / d as f64)
            })
            .collect()
    };
    let found: Vec<Vec<C>> = (0..paths).into_par_iter().filter_map(|k| track(&sys, start(k))).collect();
    let mut out: Vec<Vec<C>> = Vec::new();
    for s in found {
        let dup = out.iter().any(|o| norm(&o.iter().zip(&s).map(|(a, b)| a - b).collect::<Vec<_>>()) < 1e-6 * (1.0 + norm(&s)));
        if !dup {
            out.push(s);
        }
    }
    Ok((out, paths))
}
