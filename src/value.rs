//! Numeric values. Exact (rational) arithmetic is the default; floats appear
//! only when an operation genuinely leaves the rationals, and each fallback
//! is its own variant so it's always visible in the type.

use num_bigint::BigInt;
use num_complex::{Complex, Complex64};
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, ToPrimitive, Zero};

use crate::error::EvalError;

/// Largest integer exponent computed exactly; beyond this results are absurdly large.
const MAX_EXACT_EXPONENT: i64 = 100_000;
const MAX_FACTORIAL: u64 = 10_000;
const MAX_ROOT_INDEX: u32 = 64;

#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    /// Exact rational.
    Rational(BigRational),
    /// Approximate real.
    Real(f64),
    /// Exact complex with rational parts; imaginary part is never zero.
    ExactComplex(Complex<BigRational>),
    /// Approximate complex; imaginary part is never exactly zero.
    Complex(Complex64),
}

/// Two operands promoted to a common representation.
enum Pair {
    Rational(BigRational, BigRational),
    Real(f64, f64),
    ExactComplex(Complex<BigRational>, Complex<BigRational>),
    Complex(Complex64, Complex64),
}

pub(crate) fn ratio_to_f64(r: &BigRational) -> f64 {
    r.to_f64().unwrap_or(f64::NAN)
}

impl Value {
    pub fn int(n: i64) -> Value {
        Value::Rational(BigRational::from_integer(BigInt::from(n)))
    }

    pub fn rat(n: i64, d: i64) -> Value {
        Value::Rational(BigRational::new(BigInt::from(n), BigInt::from(d)))
    }

    pub fn is_exact(&self) -> bool {
        matches!(self, Value::Rational(_) | Value::ExactComplex(_))
    }

    pub fn is_zero(&self) -> bool {
        match self {
            Value::Rational(r) => r.is_zero(),
            Value::Real(x) => *x == 0.0,
            Value::ExactComplex(c) => c.re.is_zero() && c.im.is_zero(),
            Value::Complex(c) => c.re == 0.0 && c.im == 0.0,
        }
    }

    /// The value as an f64, if it is real.
    pub fn as_real_f64(&self) -> Option<f64> {
        match self {
            Value::Rational(r) => Some(ratio_to_f64(r)),
            Value::Real(x) => Some(*x),
            _ => None,
        }
    }

    pub fn to_c64(&self) -> Complex64 {
        match self {
            Value::Rational(r) => Complex64::new(ratio_to_f64(r), 0.0),
            Value::Real(x) => Complex64::new(*x, 0.0),
            Value::ExactComplex(c) => Complex64::new(ratio_to_f64(&c.re), ratio_to_f64(&c.im)),
            Value::Complex(c) => *c,
        }
    }

    fn to_exact_complex(&self) -> Option<Complex<BigRational>> {
        match self {
            Value::Rational(r) => Some(Complex::new(r.clone(), BigRational::zero())),
            Value::ExactComplex(c) => Some(c.clone()),
            _ => None,
        }
    }

    /// Collapse complex values with a zero imaginary part back to reals.
    pub fn normalize(self) -> Value {
        match self {
            Value::ExactComplex(c) if c.im.is_zero() => Value::Rational(c.re),
            Value::Complex(c) if c.im == 0.0 => Value::Real(c.re),
            v => v,
        }
    }

    fn promote(&self, other: &Value) -> Pair {
        match (self, other) {
            (Value::Rational(a), Value::Rational(b)) => Pair::Rational(a.clone(), b.clone()),
            (
                Value::Rational(_) | Value::ExactComplex(_),
                Value::Rational(_) | Value::ExactComplex(_),
            ) => Pair::ExactComplex(
                self.to_exact_complex().unwrap(),
                other.to_exact_complex().unwrap(),
            ),
            (Value::Rational(_) | Value::Real(_), Value::Rational(_) | Value::Real(_)) => {
                Pair::Real(self.as_real_f64().unwrap(), other.as_real_f64().unwrap())
            }
            _ => Pair::Complex(self.to_c64(), other.to_c64()),
        }
    }

    pub fn neg(&self) -> Value {
        match self {
            Value::Rational(r) => Value::Rational(-r),
            Value::Real(x) => Value::Real(-x),
            Value::ExactComplex(c) => Value::ExactComplex(-c),
            Value::Complex(c) => Value::Complex(-c),
        }
    }

    pub fn add(&self, other: &Value) -> Value {
        match self.promote(other) {
            Pair::Rational(a, b) => Value::Rational(a + b),
            Pair::Real(a, b) => Value::Real(a + b),
            Pair::ExactComplex(a, b) => Value::ExactComplex(a + b).normalize(),
            Pair::Complex(a, b) => Value::Complex(a + b).normalize(),
        }
    }

    pub fn sub(&self, other: &Value) -> Value {
        match self.promote(other) {
            Pair::Rational(a, b) => Value::Rational(a - b),
            Pair::Real(a, b) => Value::Real(a - b),
            Pair::ExactComplex(a, b) => Value::ExactComplex(a - b).normalize(),
            Pair::Complex(a, b) => Value::Complex(a - b).normalize(),
        }
    }

    pub fn mul(&self, other: &Value) -> Value {
        match self.promote(other) {
            Pair::Rational(a, b) => Value::Rational(a * b),
            Pair::Real(a, b) => Value::Real(a * b),
            Pair::ExactComplex(a, b) => Value::ExactComplex(a * b).normalize(),
            Pair::Complex(a, b) => Value::Complex(a * b).normalize(),
        }
    }

    pub fn div(&self, other: &Value) -> Result<Value, EvalError> {
        if other.is_zero() {
            return Err(EvalError::DivisionByZero);
        }
        Ok(match self.promote(other) {
            Pair::Rational(a, b) => Value::Rational(a / b),
            Pair::Real(a, b) => Value::Real(a / b),
            Pair::ExactComplex(a, b) => Value::ExactComplex(a / b).normalize(),
            Pair::Complex(a, b) => Value::Complex(a / b).normalize(),
        })
    }

    pub fn pow(&self, exp: &Value) -> Result<Value, EvalError> {
        if self.is_zero() {
            return zero_pow(exp);
        }
        if let Value::Rational(e) = exp {
            if e.is_integer() {
                return self.int_pow(&e.to_integer());
            }
            if let Value::Rational(b) = self {
                if let Some(v) = exact_rational_power(b, e) {
                    return Ok(v);
                }
            }
        }
        // No exact result: fall back to floating point.
        match (self.as_real_f64(), exp.as_real_f64()) {
            (Some(x), Some(y)) if x >= 0.0 || y.fract() == 0.0 => Ok(Value::Real(x.powf(y))),
            _ => Ok(Value::Complex(self.to_c64().powc(exp.to_c64())).normalize()),
        }
    }

    fn int_pow(&self, n: &BigInt) -> Result<Value, EvalError> {
        match self {
            Value::Real(x) => return Ok(Value::Real(x.powf(n.to_f64().unwrap_or(f64::INFINITY)))),
            Value::Complex(c) => {
                return Ok(Value::Complex(c.powf(n.to_f64().unwrap_or(f64::INFINITY))).normalize())
            }
            // ±1 stays cheap for any exponent.
            Value::Rational(b) if b.abs().is_one() => {
                let negative = b.is_negative() && n.is_odd();
                return Ok(Value::int(if negative { -1 } else { 1 }));
            }
            _ => {}
        }
        let n = n
            .to_i64()
            .filter(|n| n.abs() <= MAX_EXACT_EXPONENT)
            .ok_or_else(|| {
                EvalError::TooLarge(format!("exponent too large (max {MAX_EXACT_EXPONENT})"))
            })?;
        Ok(match self {
            Value::Rational(b) => Value::Rational(rat_ipow(b, n)),
            Value::ExactComplex(b) => Value::ExactComplex(cplx_ipow(b, n)).normalize(),
            _ => unreachable!(),
        })
    }

    pub fn factorial(&self) -> Result<Value, EvalError> {
        match self {
            Value::Rational(r) if r.is_integer() && !r.is_negative() => {
                let n = r.to_integer().to_u64().filter(|&n| n <= MAX_FACTORIAL).ok_or_else(|| {
                    EvalError::TooLarge(format!("factorial argument too large (max {MAX_FACTORIAL})"))
                })?;
                let mut acc = BigInt::one();
                for k in 2..=n {
                    acc *= k;
                }
                Ok(Value::Rational(BigRational::from_integer(acc)))
            }
            _ => Err(EvalError::Domain(
                "factorial is only defined for non-negative integers".into(),
            )),
        }
    }

    pub fn abs(&self) -> Result<Value, EvalError> {
        Ok(match self {
            Value::Rational(r) => Value::Rational(r.abs()),
            Value::Real(x) => Value::Real(x.abs()),
            Value::Complex(c) => Value::Real(c.norm()),
            Value::ExactComplex(c) => {
                // |a+bi| = sqrt(a²+b²), exact when that's a perfect square.
                let sq = &c.re * &c.re + &c.im * &c.im;
                return Value::Rational(sq).pow(&Value::rat(1, 2));
            }
        })
    }
}

fn zero_pow(exp: &Value) -> Result<Value, EvalError> {
    if exp.is_zero() {
        return Ok(Value::int(1)); // convention: 0^0 = 1
    }
    if exp.to_c64().re > 0.0 {
        Ok(Value::int(0))
    } else if exp.as_real_f64().is_some() {
        Err(EvalError::DivisionByZero)
    } else {
        Err(EvalError::Domain("0 raised to this power is undefined".into()))
    }
}

fn ipow<T: Clone + std::ops::Mul<Output = T>>(base: &T, one: T, mut e: u64) -> T {
    let mut result = one;
    let mut b = base.clone();
    while e > 0 {
        if e & 1 == 1 {
            result = result * b.clone();
        }
        e >>= 1;
        if e > 0 {
            b = b.clone() * b;
        }
    }
    result
}

/// `b^n`; `b` must be nonzero when `n` is negative.
pub(crate) fn rat_ipow(b: &BigRational, n: i64) -> BigRational {
    let r = ipow(b, BigRational::one(), n.unsigned_abs());
    if n < 0 {
        r.recip()
    } else {
        r
    }
}

fn cplx_ipow(b: &Complex<BigRational>, n: i64) -> Complex<BigRational> {
    let one = Complex::new(BigRational::one(), BigRational::zero());
    let r = ipow(b, one.clone(), n.unsigned_abs());
    if n < 0 {
        one / r
    } else {
        r
    }
}

fn exact_nth_root(x: &BigInt, q: u32) -> Option<BigInt> {
    let r = x.nth_root(q);
    (num_traits::pow(r.clone(), q as usize) == *x).then_some(r)
}

/// `b^(p/q)` when it has an exact rational or exact complex answer
/// (e.g. `8^(2/3) = 4`, `(-8)^(1/3) = -2`, `(-4)^(1/2) = 2i`).
fn exact_rational_power(b: &BigRational, e: &BigRational) -> Option<Value> {
    let q = e.denom().to_u32().filter(|&q| q <= MAX_ROOT_INDEX)?;
    let p = e.numer().to_i64().filter(|p| p.abs() <= MAX_EXACT_EXPONENT)?;
    let mag = b.abs();
    let root = BigRational::new(exact_nth_root(mag.numer(), q)?, exact_nth_root(mag.denom(), q)?);
    if !b.is_negative() {
        Some(Value::Rational(rat_ipow(&root, p)))
    } else if q % 2 == 1 {
        // Odd root of a negative number is real and negative.
        Some(Value::Rational(rat_ipow(&-root, p)))
    } else if q == 2 {
        // Principal square root of a negative number: i·√|b|.
        let base = Complex::new(BigRational::zero(), root);
        Some(Value::ExactComplex(cplx_ipow(&base, p)).normalize())
    } else {
        None
    }
}
