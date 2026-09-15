//! Stage 3 and non-linear systems: polynomial equations.
//!
//! Linear systems go to `linear`. A single polynomial equation is solved with
//! the quadratic formula (degree 2) or by factoring over Q (higher degree),
//! falling back to numeric roots only for irreducible pieces of degree ≥ 3.
//! Systems of polynomial equations are triangularized with a lexicographic
//! Gröbner basis and solved by back-substitution.

use std::collections::{BTreeSet, HashMap};

use num_complex::Complex64;
use num_rational::BigRational as Q;
use num_traits::{One, Signed, Zero};

use std::time::{Duration, Instant};

use crate::commands::{best_form, with_approx};
use crate::homotopy;
use crate::mpoly::{self, MPoly};
use crate::sym::{half, q};

/// Exact Gröbner bases get this long before square systems switch to
/// numerical homotopy continuation.
const GROEBNER_BUDGET: Duration = Duration::from_millis(400);
const PARAMETRIC_BUDGET: Duration = Duration::from_secs(5);
use crate::error::Error;
use crate::factor::Factorization;
use crate::format::fmt_c64;
use crate::linear::LinearSystem;
use crate::parser::Equation;
use crate::poly::{coeffs_in, from_uni, is_polynomial, leading, mono_div, mono_lcm, raw_term, uni_q};
use crate::simplify::simplify;
use crate::solve::{Solve, SolveError, Step};
use crate::sym::{from_expr, Atom, Mono, Sym};
use crate::upoly::{self, quadratic_roots, roots_numeric};
use crate::value::Value;
use crate::Outcome;

const MAX_BASIS: usize = 200;
const MAX_PAIRS: usize = 20_000;
const TOL: f64 = 1e-7;

#[derive(Debug, Clone)]
pub struct Root {
    pub exact: Option<Sym>,
    pub approx: Complex64,
}

impl Root {
    fn exact(s: Sym) -> Root {
        let approx = s.eval(&HashMap::new()).map(|v| v.to_c64()).unwrap_or(Complex64::new(f64::NAN, 0.0));
        Root { exact: Some(s), approx }
    }

    fn numeric(z: Complex64) -> Root {
        Root { exact: None, approx: z }
    }

    fn line(&self, var: &str) -> String {
        match &self.exact {
            Some(s) => format!("{var} = {}", with_approx(s)),
            None => format!("{var} ≈ {}", fmt_c64(self.approx)),
        }
    }

    fn value(&self) -> Value {
        if self.approx.im == 0.0 {
            Value::Real(self.approx.re)
        } else {
            Value::Complex(self.approx)
        }
    }
}

fn is_real(z: Complex64) -> bool {
    z.im.abs() <= 1e-9 * (1.0 + z.re.abs())
}

pub fn solve_equations(equations: Vec<Equation>) -> Result<Outcome, Error> {
    let mut polys = Vec::new();
    let mut denominators = Vec::new();
    for eq in &equations {
        let p = simplify(&from_expr(&eq.lhs)?.sub(&from_expr(&eq.rhs)?));
        polys.push(clear_denominators(p, &mut denominators)?);
    }
    let all: Vec<String> = polys.iter().flat_map(Sym::vars).collect::<BTreeSet<_>>().into_iter().collect();
    let (unknowns, params) = split_unknowns(&all, equations.len());
    if params.is_empty() {
        match (LinearSystem { equations: equations.clone() }).solve() {
            Ok((solution, steps)) => return Ok(Outcome::Linear(solution, steps)),
            Err(SolveError::Unsupported(_)) => {}
            Err(e) => return Err(e.into()),
        }
    }
    let mut steps = vec![Step {
        description: "Move everything to one side".into(),
        snapshot: polys.iter().map(|p| format!("{p} = 0")).collect(),
    }];
    for p in &polys {
        if !is_polynomial(p) || p.terms.keys().any(|m| m.keys().any(|a| !matches!(a, Atom::Var(_)))) {
            if params.is_empty() && polys.len() == unknowns.len() {
                return Ok(non_polynomial(&polys, &unknowns, steps));
            }
            return Err(SolveError::Unsupported(format!(
                "'{p} = 0' isn't a polynomial equation; only polynomial equations can be solved exactly, \
                 and non-polynomial ones numerically only when there are as many equations as unknowns"
            ))
            .into());
        }
    }
    if !params.is_empty() {
        return parametric(&polys, &unknowns, &params, steps);
    }
    let vars = unknowns;
    let answer = match vars.len() {
        0 => {
            if polys.iter().all(Sym::is_zero) {
                "The equation is always true.".to_string()
            } else {
                "No solution: the equations reduce to a false statement.".to_string()
            }
        }
        1 if polys.len() == 1 => {
            let x = &vars[0];
            let coeffs = uni_q(&polys[0], x).expect("checked polynomial");
            let roots = uni_roots(&coeffs, x, &mut steps, true);
            let (roots, excluded) = filter_roots(roots.into_iter().map(|(r, _)| vec![(x.clone(), r)]).collect(), &denominators);
            render_solutions(&roots, &excluded)
        }
        _ => {
            let mp: Vec<MPoly> =
                polys.iter().map(|p| MPoly::from_sym(p, &vars).expect("checked polynomial")).collect();
            let basis = match mpoly::groebner(&mp, Instant::now() + GROEBNER_BUDGET) {
                Ok(b) => b.iter().map(|p| p.to_sym(&vars)).collect::<Vec<Sym>>(),
                Err(()) if polys.len() == vars.len() => return numeric(&mp, &vars, steps, &denominators),
                Err(()) => {
                    return Err(SolveError::Unsupported("this system is too large to solve exactly".into()).into())
                }
            };
            steps.push(Step {
                description: format!(
                    "Compute a Gröbner basis (lex order {}); it's triangular, so the last polynomial has one variable",
                    vars.join(" > ")
                ),
                snapshot: basis.iter().map(|p| format!("{p} = 0")).collect(),
            });
            if basis.iter().any(|p| p.as_constant().is_some_and(|c| !c.is_zero())) {
                "No solution: the equations are inconsistent (the basis contains a nonzero constant).".to_string()
            } else {
                match back_substitute(&basis, &vars, &mut steps) {
                    Ok(sols) => {
                        let sols = verify(sols, &polys);
                        let (sols, excluded) = filter_roots(sols, &denominators);
                        render_solutions(&sols, &excluded)
                    }
                    Err(free) => format!(
                        "Infinitely many solutions: {free} is not determined by the equations \
                         (parametric solutions of non-linear systems aren't supported)."
                    ),
                }
            }
        }
    };
    Ok(Outcome::Text { steps, answer })
}

/// Multiply through by denominators (x^-k, (x+1)^-k) so the equation becomes
/// polynomial, remembering them to reject roots that make them zero.
fn clear_denominators(mut p: Sym, dens: &mut Vec<Sym>) -> Result<Sym, SolveError> {
    for _ in 0..16 {
        let worst = p
            .terms
            .keys()
            .flat_map(|m| m.iter())
            .filter(|(a, e)| e.is_negative() && matches!(a, Atom::Var(_) | Atom::Group(_)))
            .min_by(|a, b| a.1.cmp(b.1))
            .map(|(a, e)| (a.clone(), e.clone()));
        let Some((a, e)) = worst else { return Ok(p) };
        let base = match &a {
            Atom::Group(b) => (**b).clone(),
            other => Sym::atom(other.clone()),
        };
        dens.push(base);
        p = simplify(&p.mul(&raw_term(Mono::from([(a, -e)]), Q::one())));
    }
    Ok(p)
}

// ---------- univariate ----------

fn uni_roots(p: &[Q], x: &str, steps: &mut Vec<Step>, verbose: bool) -> Vec<(Root, usize)> {
    let p = upoly::trim(p.to_vec());
    let deg = upoly::degree(&p);
    if p.len() <= 1 {
        return vec![];
    }
    if deg == 1 {
        let r = -&p[0] / &p[1];
        if verbose {
            steps.push(Step {
                description: format!("Solve the linear equation for {x}"),
                snapshot: vec![format!("{x} = {r}")],
            });
        }
        return vec![(Root::exact(Sym::constant(r)), 1)];
    }
    if deg == 2 {
        let (a, b, c) = (&p[2], &p[1], &p[0]);
        let (disc, r1, r2) = quadratic_roots(a, b, c);
        if verbose {
            let kind = if disc.is_zero() {
                "D = 0: one repeated root"
            } else if disc.is_negative() {
                "D < 0: two complex roots"
            } else {
                "D > 0: two real roots"
            };
            steps.push(Step {
                description: "Identify the coefficients of a·x^2 + b·x + c = 0".into(),
                snapshot: vec![format!("a = {a}, b = {b}, c = {c}")],
            });
            steps.push(Step {
                description: format!("Discriminant D = b^2 - 4ac ({kind})"),
                snapshot: vec![format!("D = ({b})^2 - 4·({a})·({c}) = {disc}")],
            });
            steps.push(Step {
                description: "Quadratic formula: x = (-b ± sqrt(D)) / (2a)".into(),
                snapshot: if disc.is_zero() {
                    vec![format!("{x} = {r1}")]
                } else {
                    vec![format!("{x} = {r1}"), format!("{x} = {r2}")]
                },
            });
        }
        return if disc.is_zero() {
            vec![(Root::exact(r1), 2)]
        } else {
            vec![(Root::exact(r1), 1), (Root::exact(r2), 1)]
        };
    }
    let f = upoly::factor(&p);
    if verbose {
        let fact = Factorization {
            content: f.content.clone(),
            factors: f.factors.iter().map(|(g, k)| (from_uni(g, x), *k)).collect(),
            complete: f.complete,
        };
        steps.push(Step {
            description: "Factor over the rationals (rational root theorem + Kronecker's method)".into(),
            snapshot: vec![format!("{} = 0", fact.display())],
        });
    }
    let mut out = Vec::new();
    for (g, k) in &f.factors {
        match upoly::degree(g) {
            1 | 2 => out.extend(uni_roots(g, x, steps, false).into_iter().map(|(r, m)| (r, m * k))),
            n if g[1..n].iter().all(Zero::is_zero) => {
                // x^n = c: roots |c|^(1/n)·(cos θ + i·sin θ), θ = (π·[c<0] + 2πj)/n.
                let c = -&g[0] / &g[n];
                if verbose {
                    steps.push(Step {
                        description: format!("{x}^{n} = {c}: take n-th roots around the circle"),
                        snapshot: vec![format!("{x} = |{c}|^(1/{n})·(cos θ + i·sin θ)")],
                    });
                }
                let r = crate::sym::rational_power(&c.abs(), &Q::new(1.into(), (n as i64).into()));
                let pi = Sym::atom(Atom::Const(crate::sym::Konst::Pi));
                let start = if c.is_negative() { Q::new(1.into(), (n as i64).into()) } else { Q::zero() };
                for j in 0..n {
                    let t = &start + Q::new((2 * j as i64).into(), (n as i64).into());
                    let angle = pi.scale(&t);
                    let rot = crate::sym::func("cos", vec![angle.clone()])
                        .and_then(|cs| Ok(cs.add(&crate::sym::func("sin", vec![angle])?.mul(&Sym::atom(Atom::I)))));
                    match rot {
                        Ok(rot) => out.push((Root::exact(r.mul(&rot)), *k)),
                        Err(_) => {}
                    }
                }
            }
            _ => {
                if verbose {
                    steps.push(Step {
                        description: format!(
                            "No exact formula is used for degree {}: find its roots numerically",
                            upoly::degree(g)
                        ),
                        snapshot: vec![format!("{} = 0", from_uni(g, x))],
                    });
                }
                out.extend(roots_numeric(&upoly::to_c64(g)).into_iter().map(|z| (Root::numeric(z), *k)));
            }
        }
    }
    out
}

// ---------- Gröbner basis ----------

fn monic(p: &Sym) -> Sym {
    match leading(p) {
        Some((_, c)) => p.scale(&c.recip()),
        None => p.clone(),
    }
}

fn spoly(f: &Sym, g: &Sym) -> Sym {
    let (mf, cf) = leading(f).unwrap();
    let (mg, cg) = leading(g).unwrap();
    let l = mono_lcm(&mf, &mg);
    let a = raw_term(mono_div(&l, &mf).unwrap(), cf.recip());
    let b = raw_term(mono_div(&l, &mg).unwrap(), cg.recip());
    f.mul(&a).sub(&g.mul(&b))
}

fn reduce(p: &Sym, basis: &[Sym]) -> Sym {
    let leads: Vec<(Mono, Q)> = basis.iter().map(|g| leading(g).unwrap()).collect();
    let mut p = p.clone();
    let mut rem = Sym::zero();
    while let Some((m, c)) = leading(&p) {
        match leads.iter().position(|(gm, _)| mono_div(&m, gm).is_some()) {
            Some(i) => {
                let t = raw_term(mono_div(&m, &leads[i].0).unwrap(), &c / &leads[i].1);
                p = p.sub(&t.mul(&basis[i]));
            }
            None => {
                p.terms.remove(&m);
                rem.absorb(raw_term(m, c));
            }
        }
    }
    rem
}

/// Reduced Gröbner basis in lex order (Buchberger's algorithm).
pub fn groebner(polys: &[Sym]) -> Result<Vec<Sym>, SolveError> {
    let too_big = || SolveError::Unsupported("this system is too large to solve exactly".into());
    let mut g: Vec<Sym> = polys.iter().filter(|p| !p.is_zero()).map(monic).collect();
    let mut pairs: Vec<(usize, usize)> = (0..g.len()).flat_map(|j| (0..j).map(move |i| (i, j))).collect();
    let mut processed = 0;
    while let Some((i, j)) = pairs.pop() {
        processed += 1;
        if processed > MAX_PAIRS || g.len() > MAX_BASIS {
            return Err(too_big());
        }
        let (mi, _) = leading(&g[i]).unwrap();
        let (mj, _) = leading(&g[j]).unwrap();
        if mi.keys().all(|k| !mj.contains_key(k)) {
            continue; // coprime leading monomials: S-polynomial reduces to 0
        }
        let r = reduce(&spoly(&g[i], &g[j]), &g);
        if !r.is_zero() {
            let n = g.len();
            pairs.extend((0..n).map(|k| (k, n)));
            g.push(monic(&r));
        }
    }
    // Minimal basis, then fully reduce each element by the others.
    let mut minimal: Vec<Sym> = Vec::new();
    for (i, p) in g.iter().enumerate() {
        let (m, _) = leading(p).unwrap();
        let redundant = g.iter().enumerate().any(|(j, q)| {
            let (qm, _) = leading(q).unwrap();
            j != i && mono_div(&m, &qm).is_some() && (qm != m || j < i)
        });
        if !redundant {
            minimal.push(p.clone());
        }
    }
    let reduced: Vec<Sym> = (0..minimal.len())
        .map(|i| {
            let others: Vec<Sym> = minimal.iter().enumerate().filter(|(j, _)| *j != i).map(|(_, p)| p.clone()).collect();
            let (m, c) = leading(&minimal[i]).unwrap();
            let tail = minimal[i].sub(&raw_term(m.clone(), c.clone()));
            monic(&raw_term(m, c).add(&reduce(&tail, &others)))
        })
        .collect();
    let mut reduced = reduced;
    reduced.sort_by(|a, b| crate::poly::lex_cmp(&leading(b).unwrap().0, &leading(a).unwrap().0));
    Ok(reduced)
}

// ---------- back-substitution ----------

type Solution = Vec<(String, Root)>;

/// Solve the triangular basis from the last variable up. Err(var) if `var`
/// is left free (infinitely many solutions).
fn back_substitute(basis: &[Sym], vars: &[String], steps: &mut Vec<Step>) -> Result<Vec<Solution>, String> {
    let mut partial: Vec<Solution> = vec![vec![]];
    for idx in (0..vars.len()).rev() {
        let v = &vars[idx];
        let later: BTreeSet<&String> = vars[idx..].iter().collect();
        let relevant: Vec<&Sym> = basis
            .iter()
            .filter(|p| p.has_var(v) && p.vars().iter().all(|w| later.contains(w)))
            .collect();
        let mut next = Vec::new();
        for sol in &partial {
            let roots = solve_for(v, &relevant, sol, steps)?;
            for r in roots {
                let mut s = sol.clone();
                s.push((v.clone(), r));
                next.push(s);
            }
        }
        partial = next;
    }
    for s in partial.iter_mut() {
        s.reverse();
    }
    Ok(partial)
}

fn solve_for(v: &str, polys: &[&Sym], known: &Solution, steps: &mut Vec<Step>) -> Result<Vec<Root>, String> {
    if polys.is_empty() {
        return Err(v.to_string());
    }
    let log = steps.len() < 30;
    // Exact substitution when every known value is exact.
    if known.iter().all(|(_, r)| r.exact.is_some()) {
        let mut subbed = Vec::new();
        for p in polys {
            let mut s = (*p).clone();
            for (name, r) in known {
                s = s.subst(name, r.exact.as_ref().unwrap()).map_err(|_| v.to_string())?;
            }
            subbed.push(simplify(&s));
        }
        let nonzero: Vec<&Sym> = subbed.iter().filter(|s| !s.is_zero()).collect();
        if nonzero.is_empty() {
            return Err(v.to_string());
        }
        if nonzero.iter().any(|s| s.as_constant().is_some()) {
            return Ok(vec![]);
        }
        // Rational coefficients: exact univariate gcd and root finding.
        if let Some(uni) = nonzero.iter().map(|s| uni_q(s, v)).collect::<Option<Vec<_>>>() {
            let g = uni.iter().skip(1).fold(uni[0].clone(), |acc, p| upoly::gcd(&acc, p));
            if log {
                steps.push(Step {
                    description: format!("Solve for {v}{}", given(known)),
                    snapshot: vec![format!("{} = 0", from_uni(&g, v))],
                });
            }
            let mut tmp = Vec::new();
            return Ok(uni_roots(&g, v, &mut tmp, false).into_iter().map(|(r, _)| r).collect());
        }
        // Linear in v with exact (radical) coefficients.
        if let Some(c) = nonzero.iter().filter_map(|s| coeffs_in(s, v)).find(|c| c.len() == 2) {
            if let Ok(r) = c[0].neg().div(&c[1]) {
                let r = simplify(&r);
                if log {
                    steps.push(Step {
                        description: format!("Solve for {v}{}", given(known)),
                        snapshot: vec![format!("{v} = {r}")],
                    });
                }
                let root = Root::exact(r);
                let ok = nonzero.iter().all(|s| {
                    s.subst(v, root.exact.as_ref().unwrap()).map(|z| simplify(&z)).is_ok_and(|z| {
                        z.is_zero() || z.eval(&HashMap::new()).is_ok_and(|val| val.to_c64().norm() < TOL)
                    })
                });
                return Ok(if ok { vec![root] } else { vec![] });
            }
        }
    }
    // Numeric path.
    let env: HashMap<String, Value> = known.iter().map(|(n, r)| (n.clone(), r.value())).collect();
    let mut cpolys: Vec<Vec<Complex64>> = Vec::new();
    for p in polys {
        let coeffs = coeffs_in(p, v).ok_or_else(|| v.to_string())?;
        let c: Vec<Complex64> = coeffs
            .iter()
            .map(|s| s.eval(&env).map(|x| x.to_c64()).unwrap_or(Complex64::new(f64::NAN, 0.0)))
            .collect();
        cpolys.push(trim_c(c));
    }
    let nonconst: Vec<&Vec<Complex64>> = cpolys.iter().filter(|c| c.len() > 1).collect();
    if nonconst.is_empty() {
        return if cpolys.iter().all(|c| c.is_empty()) { Err(v.to_string()) } else { Ok(vec![]) };
    }
    let base = nonconst.iter().min_by_key(|c| c.len()).unwrap();
    let mut roots: Vec<Complex64> = Vec::new();
    for z in roots_numeric(base) {
        let fits = cpolys.iter().all(|c| upoly::eval_c(c, z).norm() < TOL * (1.0 + z.norm()).powi(c.len() as i32));
        if fits && !roots.iter().any(|r| (r - z).norm() < 1e-8) {
            roots.push(z);
        }
    }
    if log {
        steps.push(Step {
            description: format!("Solve for {v} numerically{}", given(known)),
            snapshot: roots.iter().map(|z| format!("{v} ≈ {}", fmt_c64(*z))).collect(),
        });
    }
    Ok(roots.into_iter().map(Root::numeric).collect())
}

fn trim_c(mut c: Vec<Complex64>) -> Vec<Complex64> {
    let scale = c.iter().map(|z| z.norm()).fold(0.0, f64::max);
    while c.last().is_some_and(|z| z.norm() <= 1e-12 * scale.max(1.0)) {
        c.pop();
    }
    c
}

fn given(known: &Solution) -> String {
    if known.is_empty() {
        return String::new();
    }
    let parts: Vec<String> = known
        .iter()
        .map(|(n, r)| match &r.exact {
            Some(s) => format!("{n} = {s}"),
            None => format!("{n} ≈ {}", fmt_c64(r.approx)),
        })
        .collect();
    format!(" with {}", parts.join(", "))
}

/// Keep solutions that satisfy the original equations; drop duplicates.
fn verify(sols: Vec<Solution>, polys: &[Sym]) -> Vec<Solution> {
    let mut out: Vec<Solution> = Vec::new();
    for s in sols {
        let env: HashMap<String, Value> = s.iter().map(|(n, r)| (n.clone(), r.value())).collect();
        let ok = polys.iter().all(|p| p.eval(&env).is_ok_and(|v| v.to_c64().norm() < 1e-6));
        let dup = out.iter().any(|o| o.iter().zip(&s).all(|((_, a), (_, b))| (a.approx - b.approx).norm() < 1e-8));
        if ok && !dup {
            out.push(s);
        }
    }
    out
}

fn filter_roots(sols: Vec<Solution>, dens: &[Sym]) -> (Vec<Solution>, Vec<Solution>) {
    sols.into_iter().partition(|s| {
        let env: HashMap<String, Value> = s.iter().map(|(n, r)| (n.clone(), r.value())).collect();
        dens.iter().all(|d| d.eval(&env).map_or(true, |v| v.to_c64().norm() > 1e-12))
    })
}

fn render_solutions(sols: &[Solution], excluded: &[Solution]) -> String {
    let mut sols: Vec<&Solution> = sols.iter().collect();
    sols.sort_by(|a, b| {
        let key = |s: &Solution| (s.iter().any(|(_, r)| !is_real(r.approx)), s.iter().map(|(_, r)| r.approx.re).collect::<Vec<_>>());
        let (ka, kb) = (key(a), key(b));
        ka.0.cmp(&kb.0).then(ka.1.partial_cmp(&kb.1).unwrap_or(std::cmp::Ordering::Equal))
    });
    let mut lines = Vec::new();
    if sols.is_empty() {
        lines.push("No solution.".to_string());
    }
    for s in &sols {
        let complex = s.iter().any(|(_, r)| !is_real(r.approx));
        let body = s.iter().map(|(n, r)| r.line(n)).collect::<Vec<_>>().join(", ");
        lines.push(if complex { format!("{body}   (complex)") } else { body });
    }
    for s in excluded {
        let body = s.iter().map(|(n, r)| r.line(n)).collect::<Vec<_>>().join(", ");
        lines.push(format!("excluded: {body} (makes a denominator zero)"));
    }
    lines.join("\n")
}

// ---------- parameters, numerics ----------

/// With more letters than equations, solve for the usual unknown names
/// (x, y, z, w, …) and treat the rest as parameters.
fn split_unknowns(vars: &[String], n_eq: usize) -> (Vec<String>, Vec<String>) {
    if vars.len() <= n_eq {
        return (vars.to_vec(), vec![]);
    }
    const PREFERRED: [&str; 8] = ["x", "y", "z", "w", "u", "v", "t", "s"];
    let mut unknowns: Vec<String> =
        PREFERRED.iter().filter(|p| vars.iter().any(|v| v == *p)).take(n_eq).map(|s| s.to_string()).collect();
    for v in vars.iter().rev() {
        if unknowns.len() >= n_eq {
            break;
        }
        if !unknowns.contains(v) {
            unknowns.push(v.clone());
        }
    }
    unknowns.sort();
    let params = vars.iter().filter(|v| !unknowns.contains(v)).cloned().collect();
    (unknowns, params)
}

/// Solve for the unknowns in terms of the parameters via a lex Gröbner basis
/// with unknowns ordered before parameters.
fn parametric(polys: &[Sym], unknowns: &[String], params: &[String], mut steps: Vec<Step>) -> Result<Outcome, Error> {
    let order: Vec<String> = unknowns.iter().chain(params).cloned().collect();
    let mp: Vec<MPoly> = polys.iter().map(|p| MPoly::from_sym(p, &order).expect("checked polynomial")).collect();
    let basis = mpoly::groebner(&mp, Instant::now() + PARAMETRIC_BUDGET).map_err(|_| {
        SolveError::Unsupported("this system with parameters is too large to solve symbolically".into())
    })?;
    let basis: Vec<Sym> = basis.iter().map(|p| p.to_sym(&order)).collect();
    steps.push(Step {
        description: format!(
            "Treat {} as parameters; compute a Gröbner basis (lex order {})",
            params.join(", "),
            order.join(" > ")
        ),
        snapshot: basis.iter().map(|p| format!("{p} = 0")).collect(),
    });
    if basis.iter().any(|p| p.as_constant().is_some_and(|c| !c.is_zero())) {
        return Ok(Outcome::Text { steps, answer: "No solution: the equations are inconsistent.".into() });
    }
    let mut lines = vec![format!("Solving for {} in terms of {}:", unknowns.join(", "), params.join(", "))];
    for p in &basis {
        if unknowns.iter().all(|u| !p.has_var(u)) {
            lines.push(format!("  requires {p} = 0"));
        }
    }
    for (k, u) in unknowns.iter().enumerate().rev() {
        let cands: Vec<(&Sym, Vec<Sym>)> = basis
            .iter()
            .filter(|p| p.has_var(u) && unknowns[..k].iter().all(|w| !p.has_var(w)))
            .filter_map(|p| coeffs_in(p, u).map(|c| (p, c)))
            .collect();
        let Some((p, c)) = cands.into_iter().min_by_key(|(p, c)| (c.len(), p.terms.len())) else {
            lines.push(format!("  {u} can be anything"));
            continue;
        };
        match c.len() - 1 {
            1 => {
                let r = simplify(&c[0].neg().div(&c[1])?);
                lines.push(format!("  {u} = {}", best_form(&r)));
            }
            2 => {
                let (a, b, cc) = (&c[2], &c[1], &c[0]);
                let disc = simplify(&b.powi(2)?.sub(&a.mul(cc).scale(&q(4))));
                let root = disc.pow_q(&half())?;
                let two_a = a.scale(&q(2));
                let r1 = simplify(&b.neg().add(&root).div(&two_a)?);
                let r2 = simplify(&b.neg().sub(&root).div(&two_a)?);
                lines.push(format!("  {u} = {r1}"));
                lines.push(format!("  or {u} = {r2}"));
            }
            d => lines.push(format!("  {u} is a root of the degree-{d} equation {p} = 0")),
        }
    }
    Ok(Outcome::Text { steps, answer: lines.join("\n") })
}

/// Continued-fraction rational approximation with a bounded denominator.
fn rational_approx(v: f64, max_den: i64) -> Option<Q> {
    if !v.is_finite() || v.abs() > 1e12 {
        return None;
    }
    let (mut h0, mut h1, mut k0, mut k1) = (0i64, 1i64, 1i64, 0i64);
    let mut x = v;
    for _ in 0..40 {
        let a = x.floor();
        let ai = a as i64;
        let (h2, k2) = (ai.checked_mul(h1)?.checked_add(h0)?, ai.checked_mul(k1)?.checked_add(k0)?);
        if k2 > max_den {
            break;
        }
        (h0, h1, k0, k1) = (h1, h2, k1, k2);
        if (h1 as f64 / k1 as f64 - v).abs() < 1e-9 * (1.0 + v.abs()) {
            return Some(Q::new(h1.into(), k1.into()));
        }
        let frac = x - a;
        if frac.abs() < 1e-15 {
            break;
        }
        x = 1.0 / frac;
    }
    None
}

/// Numeric solution; recognized as exact if it's rational and checks out exactly.
fn to_solution(x: &[Complex64], vars: &[String], mp: &[MPoly]) -> Solution {
    // Round away float noise such as 1e-32 next to much larger values.
    let scale = x.iter().map(|z| z.norm()).fold(1.0, f64::max);
    let clean = |v: f64| if v.abs() < 1e-12 * scale { 0.0 } else { v };
    let x: Vec<Complex64> = x.iter().map(|z| Complex64::new(clean(z.re), clean(z.im))).collect();
    let x = &x[..];
    let rational: Option<Vec<Q>> = x
        .iter()
        .map(|z| if z.im.abs() < 1e-8 * (1.0 + z.re.abs()) { rational_approx(z.re, 10_000) } else { None })
        .collect();
    if let Some(qs) = rational {
        if mp.iter().all(|p| p.eval_q(&qs).is_zero()) {
            return vars.iter().cloned().zip(qs.into_iter().map(|v| Root::exact(Sym::constant(v)))).collect();
        }
    }
    vars.iter().cloned().zip(x.iter().map(|z| Root::numeric(*z))).collect()
}

fn non_polynomial(polys: &[Sym], vars: &[String], mut steps: Vec<Step>) -> Outcome {
    let sols = crate::nonpoly::solve_real(polys, vars);
    steps.push(Step {
        description: "Not polynomial, so no exact method applies; search for real solutions with -10 ≤ each \
                      unknown ≤ 10 using Newton's method from many starting points in parallel (this may miss solutions)"
            .into(),
        snapshot: vec![format!("{} real solutions found", sols.len())],
    });
    let answer = if sols.is_empty() {
        "No real solution found numerically (there may still be one outside the searched region).".to_string()
    } else {
        sols.iter()
            .map(|s| {
                let body: Vec<String> =
                    vars.iter().zip(s).map(|(n, &v)| format!("{n} ≈ {}", fmt_c64(Complex64::new(v, 0.0)))).collect();
                body.join(", ")
            })
            .collect::<Vec<_>>()
            .join("\n")
    };
    Outcome::Text { steps, answer }
}

fn numeric(mp: &[MPoly], vars: &[String], mut steps: Vec<Step>, dens: &[Sym]) -> Result<Outcome, Error> {
    let (sols, paths) = homotopy::solve_system(mp, vars.len()).map_err(SolveError::Unsupported)?;
    steps.push(Step {
        description: format!(
            "An exact Gröbner basis is too expensive here; use numerical homotopy continuation \
             (deform a start system with {paths} known roots into this one, tracking all paths in parallel)"
        ),
        snapshot: vec![format!("{} finite solutions found", sols.len())],
    });
    let solutions: Vec<Solution> = sols.iter().map(|x| to_solution(x, vars, mp)).collect();
    let (s, ex) = filter_roots(solutions, dens);
    Ok(Outcome::Text { steps, answer: render_solutions(&s, &ex) })
}
