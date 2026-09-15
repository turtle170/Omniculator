//! Top-level symbolic commands: simplify, expand, factor, diff, grad,
//! jacobian, hessian (integrate lives in `integrate`).

use crate::ast::Expr;
use crate::calculus::{derivative, gradient, hessian, jacobian};
use crate::error::Error;
use crate::eval::Evaluate;
use crate::factor::factor;
use crate::format::format_result;
use crate::simplify::simplify;
use crate::solve::{SolveError, Step};
use crate::sym::{from_expr, Sym};
use crate::value::Value;
use crate::Outcome;

fn usage(u: &str) -> Error {
    Error::Solve(SolveError::Unsupported(format!("usage: {u}")))
}

fn text(steps: Vec<Step>, answer: String) -> Outcome {
    Outcome::Text { steps, answer }
}

/// Exact constant results get a decimal approximation alongside.
pub fn with_approx(s: &Sym) -> String {
    if s.vars().is_empty() && s.as_constant().is_none() {
        if let Ok(v) = s.eval(&Default::default()) {
            let approx = format_result(&v);
            let approx = approx.strip_prefix("≈ ").unwrap_or(&approx).to_string();
            return format!("{s} ≈ {approx}");
        }
    }
    s.to_string()
}

/// The shorter of the expanded and factored forms.
pub fn best_form(s: &Sym) -> String {
    let expanded = with_approx(s);
    match factor(s) {
        Ok((f, _)) if f.complete && f.factors.iter().map(|(_, k)| k).sum::<usize>() > 1 => {
            let factored = f.display();
            if factored.len() < expanded.len() {
                factored
            } else {
                expanded
            }
        }
        _ => expanded,
    }
}

pub fn var_arg(e: Option<&Expr>, f: &Sym, cmd: &str) -> Result<String, Error> {
    match e {
        Some(Expr::Var(n)) => Ok(n.clone()),
        Some(other) => Err(usage(&format!("{cmd}(): expected a variable, found '{other}'"))),
        None => {
            let vars = f.vars();
            match vars.len() {
                0 => Ok("x".into()),
                1 => Ok(vars.into_iter().next().unwrap()),
                _ => Err(usage(&format!("{cmd}(f, x): say which variable, e.g. {cmd}(f, {})", vars.iter().next().unwrap()))),
            }
        }
    }
}

fn var_list(args: &[Expr], fs: &[Sym], cmd: &str) -> Result<Vec<String>, Error> {
    if args.is_empty() {
        let mut all = std::collections::BTreeSet::new();
        for f in fs {
            all.extend(f.vars());
        }
        return Ok(all.into_iter().collect());
    }
    args.iter().map(|a| var_arg(Some(a), &Sym::zero(), cmd)).collect()
}

fn matrix(rows: &[Vec<Sym>]) -> String {
    rows.iter()
        .map(|r| format!("[{}]", r.iter().map(ToString::to_string).collect::<Vec<_>>().join(", ")))
        .collect::<Vec<_>>()
        .join("\n")
}

fn order_arg(e: Option<&Expr>) -> Result<usize, Error> {
    let Some(e) = e else { return Ok(1) };
    match e.eval()? {
        Value::Rational(r) if r.is_integer() && (1..=20).contains(&r.to_integer().try_into().unwrap_or(0i64)) => {
            Ok(r.to_integer().try_into().unwrap_or(1i64) as usize)
        }
        _ => Err(usage("diff(f, x, n) with n a whole number from 1 to 20")),
    }
}

fn superscript(n: usize) -> String {
    if n == 1 {
        String::new()
    } else {
        format!("^{n}")
    }
}

pub fn run_command(name: &str, args: &[Expr]) -> Result<Outcome, Error> {
    let sym = |i: usize| -> Result<Sym, Error> { Ok(from_expr(&args[i])?) };
    match name {
        "simplify" => {
            if args.len() != 1 {
                return Err(usage("simplify(expr)"));
            }
            Ok(text(vec![], best_form(&simplify(&sym(0)?))))
        }
        "expand" => {
            if args.len() != 1 {
                return Err(usage("expand(expr)"));
            }
            Ok(text(vec![], with_approx(&simplify(&sym(0)?))))
        }
        "factor" => {
            if args.len() != 1 {
                return Err(usage("factor(polynomial)"));
            }
            let (f, steps) = factor(&simplify(&sym(0)?))?;
            let mut answer = f.display();
            if !f.complete {
                answer += "\n(not fully factored: this may split further, but it's beyond what \
                           the factoring search can prove)";
            }
            Ok(text(steps, answer))
        }
        "diff" | "derivative" => {
            if args.is_empty() || args.len() > 3 {
                return Err(usage("diff(f, x) or diff(f, x, n)"));
            }
            let f = sym(0)?;
            let x = var_arg(args.get(1), &f, "diff")?;
            let n = order_arg(args.get(2))?;
            let (d, steps) = derivative(&f, &x, n)?;
            let sup = superscript(n);
            Ok(text(steps, format!("d{sup}/d{x}{sup} [{f}] = {}", with_approx(&d))))
        }
        "grad" | "gradient" => {
            if args.is_empty() {
                return Err(usage("grad(f) or grad(f, x, y, ...)"));
            }
            let f = sym(0)?;
            let vars = var_list(&args[1..], std::slice::from_ref(&f), "grad")?;
            let g = gradient(&f, &vars)?;
            let lines: Vec<String> = vars.iter().zip(&g).map(|(v, d)| format!("∂f/∂{v} = {d}")).collect();
            let vec = g.iter().map(ToString::to_string).collect::<Vec<_>>().join(", ");
            Ok(text(vec![], format!("{}\n∇f = ({vec})", lines.join("\n"))))
        }
        "jacobian" => {
            if args.is_empty() {
                return Err(usage("jacobian(f1, f2, ...)"));
            }
            let fs = args.iter().map(from_expr).collect::<Result<Vec<_>, _>>()?;
            let vars = var_list(&[], &fs, "jacobian")?;
            let j = jacobian(&fs, &vars)?;
            Ok(text(vec![], format!("variables: {}\n{}", vars.join(", "), matrix(&j))))
        }
        "hessian" => {
            if args.is_empty() {
                return Err(usage("hessian(f) or hessian(f, x, y, ...)"));
            }
            let f = sym(0)?;
            let vars = var_list(&args[1..], std::slice::from_ref(&f), "hessian")?;
            let h = hessian(&f, &vars)?;
            Ok(text(vec![], format!("variables: {}\n{}", vars.join(", "), matrix(&h))))
        }
        "integrate" | "int" => crate::integrate::command(args),
        _ => unreachable!("not a command: {name}"),
    }
}
