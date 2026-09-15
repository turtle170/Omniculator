use std::fmt;

use num_traits::Signed;

use crate::format::format_value;
use crate::value::Value;

#[derive(Debug, Clone, PartialEq)]
pub enum Expr {
    Num(Value),
    Var(String),
    Add(Box<Expr>, Box<Expr>),
    Sub(Box<Expr>, Box<Expr>),
    Mul(Box<Expr>, Box<Expr>),
    Div(Box<Expr>, Box<Expr>),
    Pow(Box<Expr>, Box<Expr>),
    Neg(Box<Expr>),
    Factorial(Box<Expr>),
    Call(String, Vec<Expr>),
}

/// Display precedence: higher binds tighter.
fn prec(e: &Expr) -> u8 {
    match e {
        Expr::Add(..) | Expr::Sub(..) => 1,
        Expr::Mul(..) | Expr::Div(..) => 2,
        Expr::Neg(_) => 3,
        Expr::Pow(..) => 4,
        Expr::Factorial(_) => 5,
        Expr::Num(Value::Rational(r)) if r.is_negative() => 3,
        Expr::Num(Value::Rational(r)) if !r.is_integer() => 2,
        Expr::Num(Value::Real(x)) if *x < 0.0 => 3,
        Expr::Num(Value::ExactComplex(_) | Value::Complex(_)) => 1,
        Expr::Num(_) | Expr::Var(_) | Expr::Call(..) => 6,
    }
}

fn child(f: &mut fmt::Formatter, e: &Expr, min: u8) -> fmt::Result {
    if prec(e) < min {
        write!(f, "({e})")
    } else {
        write!(f, "{e}")
    }
}

impl fmt::Display for Expr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let binary = |f: &mut fmt::Formatter, a, op, b, l, r| {
            child(f, a, l)?;
            f.write_str(op)?;
            child(f, b, r)
        };
        match self {
            Expr::Num(v) => f.write_str(&format_value(v)),
            Expr::Var(n) => f.write_str(n),
            Expr::Add(a, b) => binary(f, a, " + ", b, 1, 2),
            Expr::Sub(a, b) => binary(f, a, " - ", b, 1, 2),
            Expr::Mul(a, b) => binary(f, a, " * ", b, 2, 3),
            Expr::Div(a, b) => binary(f, a, " / ", b, 2, 3),
            Expr::Pow(a, b) => binary(f, a, "^", b, 5, 4),
            Expr::Neg(a) => {
                f.write_str("-")?;
                child(f, a, 4)
            }
            Expr::Factorial(a) => {
                child(f, a, 6)?;
                f.write_str("!")
            }
            Expr::Call(name, args) => {
                write!(f, "{name}(")?;
                for (i, a) in args.iter().enumerate() {
                    if i > 0 {
                        f.write_str(", ")?;
                    }
                    write!(f, "{a}")?;
                }
                f.write_str(")")
            }
        }
    }
}
