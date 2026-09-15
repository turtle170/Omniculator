//! Fast-path evaluation: no step tracking.

use std::f64::consts::{E, PI, TAU};

use num_complex::{Complex, Complex64 as C64};
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};

use crate::ast::Expr;
use crate::builtins::{is_command, CONSTANTS, FUNCTIONS};
use crate::error::{EvalError, NameKind};
use crate::suggest::suggest;
use crate::value::{rat_ipow, ratio_to_f64, Value};

pub trait Evaluate {
    fn eval(&self) -> Result<Value, EvalError>;
}

impl Evaluate for Expr {
    fn eval(&self) -> Result<Value, EvalError> {
        let v = match self {
            Expr::Num(v) => v.clone(),
            Expr::Var(name) => constant(name).ok_or_else(|| {
                let all: Vec<&str> = CONSTANTS.iter().chain(FUNCTIONS).copied().collect();
                EvalError::UnknownName {
                    name: name.clone(),
                    kind: NameKind::Variable,
                    suggestion: suggest(name, &all),
                }
            })?,
            Expr::Add(a, b) => a.eval()?.add(&b.eval()?),
            Expr::Sub(a, b) => a.eval()?.sub(&b.eval()?),
            Expr::Mul(a, b) => a.eval()?.mul(&b.eval()?),
            Expr::Div(a, b) => a.eval()?.div(&b.eval()?)?,
            Expr::Pow(a, b) => a.eval()?.pow(&b.eval()?)?,
            Expr::Neg(a) => a.eval()?.neg(),
            Expr::Factorial(a) => a.eval()?.factorial()?,
            Expr::Call(name, args) => {
                let vals = args.iter().map(Evaluate::eval).collect::<Result<Vec<_>, _>>()?;
                call(name, &vals)?
            }
        };
        finite(v)
    }
}

fn finite(v: Value) -> Result<Value, EvalError> {
    let ok = match &v {
        Value::Real(x) => x.is_finite(),
        Value::Complex(c) => c.re.is_finite() && c.im.is_finite(),
        _ => true,
    };
    if ok {
        Ok(v)
    } else {
        Err(EvalError::Domain("result is not a finite number (overflow or undefined)".into()))
    }
}

fn constant(name: &str) -> Option<Value> {
    Some(match name {
        "pi" => Value::Real(PI),
        "e" => Value::Real(E),
        "tau" => Value::Real(TAU),
        "phi" => Value::Real((1.0 + 5f64.sqrt()) / 2.0),
        "i" => Value::ExactComplex(Complex::new(BigRational::zero(), BigRational::one())),
        _ => return None,
    })
}

fn arity(name: &str, args: &[Value], allowed: &[usize]) -> Result<(), EvalError> {
    if allowed.contains(&args.len()) {
        return Ok(());
    }
    let expected = allowed.iter().map(ToString::to_string).collect::<Vec<_>>().join(" or ");
    Err(EvalError::Arity { name: name.to_string(), expected, got: args.len() })
}

pub(crate) fn call(name: &str, args: &[Value]) -> Result<Value, EvalError> {
    match name {
        c if is_command(c) => Err(EvalError::Domain(format!(
            "{c}() must be used on its own, e.g. {c}(x^2, x)"
        ))),
        "sqrt" => {
            arity(name, args, &[1])?;
            args[0].pow(&Value::rat(1, 2))
        }
        "cbrt" => {
            arity(name, args, &[1])?;
            args[0].pow(&Value::rat(1, 3))
        }
        "root" => {
            arity(name, args, &[2])?;
            if args[1].is_zero() {
                return Err(EvalError::Domain("root index can't be 0".into()));
            }
            args[0].pow(&Value::int(1).div(&args[1])?)
        }
        "abs" => {
            arity(name, args, &[1])?;
            args[0].abs()
        }
        "sin" | "cos" | "tan" | "asin" | "acos" | "atan" | "sinh" | "cosh" | "tanh" | "exp" => {
            arity(name, args, &[1])?;
            Ok(transcendental(name, &args[0]))
        }
        "ln" => {
            arity(name, args, &[1])?;
            log_base(&args[0], None)
        }
        "log" => {
            arity(name, args, &[1, 2])?;
            let base = args.get(1).cloned().unwrap_or(Value::int(10));
            log_base(&args[0], Some(&base))
        }
        "log2" => {
            arity(name, args, &[1])?;
            log_base(&args[0], Some(&Value::int(2)))
        }
        "floor" | "ceil" | "round" => {
            arity(name, args, &[1])?;
            rounding(name, &args[0])
        }
        "re" | "im" | "conj" | "arg" => {
            arity(name, args, &[1])?;
            complex_part(name, &args[0])
        }
        _ => Err(EvalError::UnknownName {
            name: name.to_string(),
            kind: NameKind::Function,
            suggestion: suggest(name, FUNCTIONS),
        }),
    }
}

type RealFn = fn(f64) -> f64;
type ComplexFn = fn(C64) -> C64;
type Domain = fn(f64) -> bool;

fn transcendental(name: &str, v: &Value) -> Value {
    // Exact special values, so e.g. sin(0) stays exactly 0.
    if v.is_zero() {
        match name {
            "sin" | "tan" | "sinh" | "tanh" | "asin" | "atan" => return Value::int(0),
            "cos" | "cosh" | "exp" => return Value::int(1),
            _ => {}
        }
    }
    if name == "acos" && *v == Value::int(1) {
        return Value::int(0);
    }
    let all: Domain = |_| true;
    let unit: Domain = |x| x.abs() <= 1.0;
    let (real, complex, domain): (RealFn, ComplexFn, Domain) = match name {
        "sin" => (f64::sin, |c: C64| c.sin(), all),
        "cos" => (f64::cos, |c: C64| c.cos(), all),
        "tan" => (f64::tan, |c: C64| c.tan(), all),
        "asin" => (f64::asin, |c: C64| c.asin(), unit),
        "acos" => (f64::acos, |c: C64| c.acos(), unit),
        "atan" => (f64::atan, |c: C64| c.atan(), all),
        "sinh" => (f64::sinh, |c: C64| c.sinh(), all),
        "cosh" => (f64::cosh, |c: C64| c.cosh(), all),
        "tanh" => (f64::tanh, |c: C64| c.tanh(), all),
        "exp" => (f64::exp, |c: C64| c.exp(), all),
        _ => unreachable!("not a transcendental function: {name}"),
    };
    match v.as_real_f64() {
        Some(x) if domain(x) => Value::Real(real(x)),
        _ => Value::Complex(complex(v.to_c64())).normalize(),
    }
}

/// Natural log when `base` is None.
fn log_base(x: &Value, base: Option<&Value>) -> Result<Value, EvalError> {
    if x.is_zero() {
        return Err(EvalError::Domain("logarithm of 0 is undefined".into()));
    }
    match base {
        Some(b) => {
            if b.is_zero() || *b == Value::int(1) {
                return Err(EvalError::Domain("logarithm base can't be 0 or 1".into()));
            }
            if let (Value::Rational(xr), Value::Rational(br)) = (x, b) {
                if let Some(k) = exact_log(xr, br) {
                    return Ok(Value::int(k));
                }
            }
        }
        None if *x == Value::int(1) => return Ok(Value::int(0)),
        None => {}
    }
    let ln = |v: &Value| match v.as_real_f64() {
        Some(r) if r > 0.0 => Value::Real(r.ln()),
        _ => Value::Complex(v.to_c64().ln()).normalize(),
    };
    match base {
        None => Ok(ln(x)),
        Some(b) => ln(x).div(&ln(b)),
    }
}

/// Integer k with b^k = x exactly, if one exists (log(1000) = 3, log(1/8, 2) = -3).
fn exact_log(x: &BigRational, b: &BigRational) -> Option<i64> {
    if !x.is_positive() || !b.is_positive() {
        return None;
    }
    let est = ratio_to_f64(x).ln() / ratio_to_f64(b).ln();
    let k = est.round();
    if !est.is_finite() || (est - k).abs() > 1e-9 || k.abs() > 10_000.0 {
        return None;
    }
    let k = k as i64;
    (rat_ipow(b, k) == *x).then_some(k)
}

fn rounding(name: &str, v: &Value) -> Result<Value, EvalError> {
    Ok(match v {
        Value::Rational(r) => Value::Rational(match name {
            "floor" => r.floor(),
            "ceil" => r.ceil(),
            _ => r.round(),
        }),
        Value::Real(x) => Value::Real(match name {
            "floor" => x.floor(),
            "ceil" => x.ceil(),
            _ => x.round(),
        }),
        _ => return Err(EvalError::Domain(format!("{name}() needs a real number"))),
    })
}

fn complex_part(name: &str, v: &Value) -> Result<Value, EvalError> {
    Ok(match (name, v) {
        ("re", Value::ExactComplex(c)) => Value::Rational(c.re.clone()),
        ("re", Value::Complex(c)) => Value::Real(c.re),
        ("im", Value::ExactComplex(c)) => Value::Rational(c.im.clone()),
        ("im", Value::Complex(c)) => Value::Real(c.im),
        ("im", _) => Value::int(0),
        ("conj", Value::ExactComplex(c)) => Value::ExactComplex(c.conj()),
        ("conj", Value::Complex(c)) => Value::Complex(c.conj()),
        ("re" | "conj", _) => v.clone(),
        _ => {
            if v.is_zero() {
                return Err(EvalError::Domain("arg(0) is undefined".into()));
            }
            match v {
                Value::Rational(r) if r.is_positive() => Value::int(0),
                _ => Value::Real(v.to_c64().arg()),
            }
        }
    })
}
