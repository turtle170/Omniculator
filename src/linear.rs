//! Stage 2: systems of linear equations, solved by Gauss-Jordan elimination
//! over exact values, with a step trace.

use std::collections::HashMap;
use std::fmt;

use num_traits::Signed;

use crate::ast::Expr;
use crate::builtins::is_constant;
use crate::eval::Evaluate;
use crate::format::{format_result, format_value};
use crate::parser::Equation;
use crate::solve::{Solve, SolveError, Step};
use crate::value::Value;

pub struct LinearSystem {
    pub equations: Vec<Equation>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Solution {
    /// One value per variable, in order of first appearance.
    Unique(Vec<(String, Value)>),
    /// Each dependent variable expressed in terms of the free ones.
    Infinite { dependent: Vec<(String, LinearExpr)>, free: Vec<String> },
    None { reason: String },
}

/// `constant + c₁·v₁ + c₂·v₂ + …`
#[derive(Debug, Clone, PartialEq)]
pub struct LinearExpr {
    pub terms: Vec<(String, Value)>,
    pub constant: Value,
}

impl fmt::Display for LinearExpr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str(&fmt_combo(self.terms.iter().map(|(v, c)| (v.as_str(), c)), &self.constant))
    }
}

/// Coefficients of a linear expression while it's being extracted.
struct Lin {
    coeffs: HashMap<String, Value>,
    constant: Value,
}

impl Lin {
    fn constant(v: Value) -> Lin {
        Lin { coeffs: HashMap::new(), constant: v }
    }

    fn is_constant(&self) -> bool {
        self.coeffs.values().all(Value::is_zero)
    }

    fn scale(mut self, k: &Value) -> Lin {
        for v in self.coeffs.values_mut() {
            *v = v.mul(k);
        }
        self.constant = self.constant.mul(k);
        self
    }

    fn add(mut self, other: Lin) -> Lin {
        for (name, v) in other.coeffs {
            let entry = self.coeffs.entry(name).or_insert(Value::int(0));
            *entry = entry.add(&v);
        }
        self.constant = self.constant.add(&other.constant);
        self
    }
}

fn has_vars(e: &Expr) -> bool {
    match e {
        Expr::Num(_) => false,
        Expr::Var(n) => !is_constant(n),
        Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) | Expr::Pow(a, b) => {
            has_vars(a) || has_vars(b)
        }
        Expr::Neg(a) | Expr::Factorial(a) => has_vars(a),
        Expr::Call(_, args) => args.iter().any(has_vars),
    }
}

fn nonlinear(e: &Expr) -> SolveError {
    SolveError::Unsupported(format!(
        "'{e}' is not linear; only linear equations can be solved so far"
    ))
}

/// Extract coefficients, recording variables in order of first appearance.
fn linearize(e: &Expr, order: &mut Vec<String>) -> Result<Lin, SolveError> {
    if !has_vars(e) {
        return Ok(Lin::constant(e.eval()?));
    }
    let minus_one = Value::int(-1);
    match e {
        Expr::Var(n) => {
            if !order.contains(n) {
                order.push(n.clone());
            }
            Ok(Lin { coeffs: HashMap::from([(n.clone(), Value::int(1))]), constant: Value::int(0) })
        }
        Expr::Add(a, b) => Ok(linearize(a, order)?.add(linearize(b, order)?)),
        Expr::Sub(a, b) => Ok(linearize(a, order)?.add(linearize(b, order)?.scale(&minus_one))),
        Expr::Neg(a) => Ok(linearize(a, order)?.scale(&minus_one)),
        Expr::Mul(a, b) => {
            let (la, lb) = (linearize(a, order)?, linearize(b, order)?);
            if la.is_constant() {
                Ok(lb.scale(&la.constant))
            } else if lb.is_constant() {
                Ok(la.scale(&lb.constant))
            } else {
                Err(nonlinear(e))
            }
        }
        Expr::Div(a, b) => {
            let lb = linearize(b, order)?;
            if !lb.is_constant() {
                return Err(nonlinear(e));
            }
            let inv = Value::int(1).div(&lb.constant)?;
            Ok(linearize(a, order)?.scale(&inv))
        }
        Expr::Pow(base, exp) if !has_vars(exp) => {
            let n = exp.eval()?;
            if n == Value::int(1) {
                linearize(base, order)
            } else if n.is_zero() {
                Ok(Lin::constant(Value::int(1)))
            } else {
                Err(nonlinear(e))
            }
        }
        _ => Err(nonlinear(e)),
    }
}

/// Zero test that tolerates float noise in approximate values.
fn negligible(v: &Value) -> bool {
    if v.is_exact() {
        v.is_zero()
    } else {
        v.to_c64().norm() < 1e-10
    }
}

/// First nonzero entry when everything is exact; largest magnitude otherwise
/// (partial pivoting, for numerical stability).
fn pick_pivot(rows: &[Vec<Value>], from: usize, col: usize) -> Option<usize> {
    let cands: Vec<usize> = (from..rows.len()).filter(|&i| !negligible(&rows[i][col])).collect();
    if cands.iter().all(|&i| rows[i][col].is_exact()) {
        cands.first().copied()
    } else {
        let norm = |i: usize| rows[i][col].to_c64().norm();
        cands.into_iter().max_by(|&a, &b| norm(a).total_cmp(&norm(b)))
    }
}

/// Sign and body of one term (`var` None means a bare constant).
fn term(c: &Value, var: Option<&str>) -> (bool, String) {
    let (negative, mag, simple) = match c {
        Value::Rational(r) => (r.is_negative(), Value::Rational(r.abs()), r.is_integer()),
        Value::Real(x) => (*x < 0.0, Value::Real(x.abs()), true),
        _ => (false, c.clone(), false),
    };
    let m = format_value(&mag);
    let body = match var {
        None if simple || matches!(c, Value::Rational(_)) => m,
        None => format!("({m})"),
        Some(v) if mag == Value::int(1) => v.to_string(),
        Some(v) if simple => format!("{m}{v}"),
        Some(v) => format!("({m}){v}"),
    };
    (negative, body)
}

fn push_term(out: &mut String, negative: bool, body: &str) {
    if out.is_empty() {
        if negative {
            out.push('-');
        }
    } else {
        out.push_str(if negative { " - " } else { " + " });
    }
    out.push_str(body);
}

fn fmt_combo<'a>(terms: impl Iterator<Item = (&'a str, &'a Value)>, constant: &Value) -> String {
    let mut out = String::new();
    if !constant.is_zero() {
        let (neg, body) = term(constant, None);
        push_term(&mut out, neg, &body);
    }
    for (var, c) in terms {
        if !c.is_zero() {
            let (neg, body) = term(c, Some(var));
            push_term(&mut out, neg, &body);
        }
    }
    if out.is_empty() {
        out.push('0');
    }
    out
}

fn snapshot(rows: &[Vec<Value>], vars: &[String]) -> Vec<String> {
    let n = vars.len();
    rows.iter()
        .map(|row| {
            let lhs = fmt_combo(vars.iter().map(String::as_str).zip(&row[..n]), &Value::int(0));
            format!("{lhs} = {}", format_value(&row[n]))
        })
        .collect()
}

fn paren(v: &Value) -> String {
    let s = format_value(v);
    if s.contains(['/', ' ', '-']) {
        format!("({s})")
    } else {
        s
    }
}

/// `R2 ← R2 - 3·R1`
fn elim_op(target: usize, pivot: usize, factor: &Value) -> String {
    let (negative, mag) = term(factor, None);
    let scaled = if mag == "1" { format!("R{}", pivot + 1) } else { format!("{mag}·R{}", pivot + 1) };
    let sign = if negative { '+' } else { '-' };
    format!("R{0} ← R{0} {sign} {scaled}", target + 1)
}

impl Solve for LinearSystem {
    type Output = Solution;

    fn solve(&self) -> Result<(Solution, Vec<Step>), SolveError> {
        let mut vars = Vec::new();
        let mut lins = Vec::new();
        for eq in &self.equations {
            // lhs = rhs  →  lhs - rhs = 0
            let lhs = linearize(&eq.lhs, &mut vars)?;
            let rhs = linearize(&eq.rhs, &mut vars)?;
            lins.push(lhs.add(rhs.scale(&Value::int(-1))));
        }
        let n = vars.len();
        // Augmented matrix: [coefficients | right-hand side].
        let mut rows: Vec<Vec<Value>> = lins
            .iter()
            .map(|l| {
                let mut row: Vec<Value> =
                    vars.iter().map(|v| l.coeffs.get(v).cloned().unwrap_or(Value::int(0))).collect();
                row.push(l.constant.neg());
                row
            })
            .collect();

        let mut steps = vec![Step {
            description: "Write each equation in standard form".into(),
            snapshot: snapshot(&rows, &vars),
        }];
        let mut pivots: Vec<usize> = Vec::new();
        for col in 0..n {
            let r = pivots.len();
            if r == rows.len() {
                break;
            }
            let Some(p) = pick_pivot(&rows, r, col) else { continue };
            let mut ops = Vec::new();
            if p != r {
                rows.swap(p, r);
                ops.push(format!("swap R{} and R{}", r + 1, p + 1));
            }
            let pv = rows[r][col].clone();
            if pv != Value::int(1) {
                let inv = Value::int(1).div(&pv)?;
                for x in rows[r].iter_mut() {
                    *x = x.mul(&inv);
                }
                ops.push(format!("R{0} ← R{0} ÷ {1}", r + 1, paren(&pv)));
            }
            rows[r][col] = Value::int(1);
            for i in 0..rows.len() {
                if i == r || negligible(&rows[i][col]) {
                    continue;
                }
                let factor = rows[i][col].clone();
                let pivot_row = rows[r].clone();
                for (x, pr) in rows[i].iter_mut().zip(&pivot_row) {
                    *x = x.sub(&factor.mul(pr));
                }
                rows[i][col] = Value::int(0);
                ops.push(elim_op(i, r, &factor));
            }
            pivots.push(col);
            if !ops.is_empty() {
                steps.push(Step {
                    description: format!("Isolate {}: {}", vars[col], ops.join("; ")),
                    snapshot: snapshot(&rows, &vars),
                });
            }
        }

        if let Some(bad) =
            rows.iter().find(|row| row[..n].iter().all(negligible) && !negligible(&row[n]))
        {
            let reason = format!(
                "the equations contradict each other (they reduce to 0 = {})",
                format_value(&bad[n])
            );
            return Ok((Solution::None { reason }, steps));
        }

        if pivots.len() == n {
            let values = vars.into_iter().enumerate().map(|(i, v)| (v, rows[i][n].clone())).collect();
            return Ok((Solution::Unique(values), steps));
        }

        let free_cols: Vec<usize> = (0..n).filter(|c| !pivots.contains(c)).collect();
        let dependent = pivots
            .iter()
            .enumerate()
            .map(|(i, &c)| {
                let terms =
                    free_cols.iter().map(|&j| (vars[j].clone(), rows[i][j].neg())).collect();
                (vars[c].clone(), LinearExpr { terms, constant: rows[i][n].clone() })
            })
            .collect();
        let free = free_cols.iter().map(|&j| vars[j].clone()).collect();
        Ok((Solution::Infinite { dependent, free }, steps))
    }
}

pub fn format_solution(solution: &Solution, steps: &[Step]) -> String {
    let mut out = String::from("Steps:\n");
    for (i, s) in steps.iter().enumerate() {
        out += &format!("  {}. {}\n", i + 1, s.description);
        for line in &s.snapshot {
            out += &format!("       {line}\n");
        }
    }
    match solution {
        Solution::Unique(values) if values.is_empty() => out += "The equation is always true.",
        Solution::Unique(values) => {
            let lines: Vec<String> = values
                .iter()
                .map(|(var, v)| {
                    let r = format_result(v);
                    match r.strip_prefix("≈ ") {
                        Some(approx) => format!("{var} ≈ {approx}"),
                        None => format!("{var} = {r}"),
                    }
                })
                .collect();
            out += &lines.join("\n");
        }
        Solution::Infinite { dependent, free } => {
            out += &format!("Infinitely many solutions ({} can be anything):", free.join(", "));
            for (var, e) in dependent {
                out += &format!("\n  {var} = {e}");
            }
        }
        Solution::None { reason } => out += &format!("No solution: {reason}."),
    }
    out
}
