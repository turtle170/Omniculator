//! Omniculator: a deterministic math engine. No AI in the compute path.

pub mod ast;
pub mod builtins;
pub mod calculus;
pub mod commands;
pub mod error;
pub mod eval;
pub mod factor;
pub mod format;
pub mod integrate;
pub mod linear;
pub mod parser;
pub mod poly;
pub mod polysolve;
pub mod simplify;
pub mod solve;
pub mod suggest;
pub mod sym;
pub mod sym_fmt;
pub mod token;
pub mod upoly;
pub mod value;

pub use ast::Expr;
pub use error::{Error, EvalError, ParseError};
pub use eval::Evaluate;
pub use value::Value;

use error::NameKind;
use linear::Solution;
use parser::Input;
use solve::Step;
use sym::Sym;

#[derive(Debug)]
pub enum Outcome {
    Value(Value),
    /// An irrational constant with its exact form, e.g. `2sqrt(2)` ≈ 2.83.
    Exact(Sym, Value),
    Linear(Solution, Vec<Step>),
    Text { steps: Vec<Step>, answer: String },
}

/// Parse and evaluate a plain expression (the fast path, no step tracking).
pub fn calculate(input: &str) -> Result<Value, Error> {
    let expr = parser::parse(input)?;
    Ok(expr.eval()?)
}

/// Evaluate an expression, run a command, or solve equations.
pub fn run(input: &str) -> Result<Outcome, Error> {
    match parser::parse_input(input)? {
        Input::Expr(Expr::Call(name, args)) if builtins::is_command(&name) => {
            commands::run_command(&name, &args)
        }
        Input::Expr(e) => match e.eval() {
            Ok(v) if v.is_exact() => Ok(Outcome::Value(v)),
            Ok(v) => Ok(match sym::from_expr(&e) {
                Ok(s) if s.vars().is_empty() => match s.as_constant() {
                    Some(c) => Outcome::Value(Value::Rational(c)), // e.g. sin(pi/6) = 1/2 exactly
                    None => Outcome::Exact(s, v),
                },
                _ => Outcome::Value(v),
            }),
            // Free variables with no typo suggestion: simplify symbolically.
            Err(EvalError::UnknownName { kind: NameKind::Variable, suggestion: None, .. }) => {
                let s = simplify::simplify(&sym::from_expr(&e)?);
                Ok(Outcome::Text { steps: vec![], answer: commands::best_form(&s) })
            }
            Err(err) => Err(err.into()),
        },
        Input::Equations(equations) => polysolve::solve_equations(equations),
    }
}
