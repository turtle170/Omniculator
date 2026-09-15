//! Slow path for "hard problems": solving with a step trace.

use std::fmt;

use crate::error::EvalError;

#[derive(Debug, Clone, PartialEq)]
pub struct Step {
    pub description: String,
    /// The state after this step, one line per equation/row.
    pub snapshot: Vec<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum SolveError {
    /// Outside what the solver can handle; reported honestly, never guessed.
    Unsupported(String),
    Eval(EvalError),
}

impl From<EvalError> for SolveError {
    fn from(e: EvalError) -> Self {
        SolveError::Eval(e)
    }
}

impl fmt::Display for SolveError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            SolveError::Unsupported(m) => f.write_str(m),
            SolveError::Eval(e) => e.fmt(f),
        }
    }
}

pub trait Solve {
    type Output;
    fn solve(&self) -> Result<(Self::Output, Vec<Step>), SolveError>;
}
