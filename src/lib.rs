//! Omniculator: a deterministic math engine. No AI in the compute path.

pub mod ast;
pub mod builtins;
pub mod error;
pub mod eval;
pub mod format;
pub mod linear;
pub mod parser;
pub mod solve;
pub mod suggest;
pub mod token;
pub mod value;

pub use ast::Expr;
pub use error::{Error, EvalError, ParseError};
pub use eval::Evaluate;
pub use value::Value;

use linear::{LinearSystem, Solution};
use parser::Input;
use solve::{Solve, Step};

#[derive(Debug)]
pub enum Outcome {
    Value(Value),
    Linear(Solution, Vec<Step>),
}

/// Parse and evaluate a plain expression (the fast path, no step tracking).
pub fn calculate(input: &str) -> Result<Value, Error> {
    let expr = parser::parse(input)?;
    Ok(expr.eval()?)
}

/// Evaluate an expression, or solve a system of equations with steps.
pub fn run(input: &str) -> Result<Outcome, Error> {
    match parser::parse_input(input)? {
        Input::Expr(e) => Ok(Outcome::Value(e.eval()?)),
        Input::Equations(equations) => {
            let (solution, steps) = LinearSystem { equations }.solve()?;
            Ok(Outcome::Linear(solution, steps))
        }
    }
}
