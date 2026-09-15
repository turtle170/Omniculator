use std::fmt;

use crate::solve::SolveError;

#[derive(Debug, Clone, PartialEq)]
pub struct ParseError {
    pub message: String,
    /// Byte span in the input, for pointing at the problem.
    pub start: usize,
    pub end: usize,
}

impl ParseError {
    pub fn new(message: impl Into<String>, start: usize, end: usize) -> Self {
        ParseError { message: message.into(), start, end }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NameKind {
    Function,
    Variable,
}

#[derive(Debug, Clone, PartialEq)]
pub enum EvalError {
    DivisionByZero,
    Domain(String),
    UnknownName { name: String, kind: NameKind, suggestion: Option<String> },
    Arity { name: String, expected: String, got: usize },
    TooLarge(String),
}

#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    Parse(ParseError),
    Eval(EvalError),
    Solve(SolveError),
}

impl From<SolveError> for Error {
    fn from(e: SolveError) -> Self {
        Error::Solve(e)
    }
}

impl From<ParseError> for Error {
    fn from(e: ParseError) -> Self {
        Error::Parse(e)
    }
}

impl From<EvalError> for Error {
    fn from(e: EvalError) -> Self {
        Error::Eval(e)
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl fmt::Display for EvalError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            EvalError::DivisionByZero => f.write_str("division by zero"),
            EvalError::Domain(m) | EvalError::TooLarge(m) => f.write_str(m),
            EvalError::UnknownName { name, kind, suggestion } => {
                let what = match kind {
                    NameKind::Function => "function",
                    NameKind::Variable => "name",
                };
                write!(f, "unknown {what} '{name}'")?;
                if let Some(s) = suggestion {
                    write!(f, " (did you mean '{s}'?)")?;
                }
                Ok(())
            }
            EvalError::Arity { name, expected, got } => {
                write!(f, "{name}() takes {expected} argument(s), got {got}")
            }
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Error::Parse(e) => e.fmt(f),
            Error::Eval(e) => e.fmt(f),
            Error::Solve(e) => e.fmt(f),
        }
    }
}

impl std::error::Error for Error {}
