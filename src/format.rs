//! Output formatting: exact fraction first, decimal alongside.

use num_bigint::BigInt;
use num_complex::{Complex, Complex64};
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};

use crate::value::{ratio_to_f64, Value};

const DECIMALS: usize = 12;

/// Full result line: `1/3 ≈ 0.333333333333`, `1/4 = 0.25`, `≈ 1.41421356237`.
pub fn format_result(v: &Value) -> String {
    match v {
        Value::Rational(r) if r.is_integer() => r.to_string(),
        Value::Rational(r) => match decimal(r) {
            Some((s, true)) => format!("{r} = {s}"),
            Some((s, false)) => format!("{r} ≈ {s}"),
            None => format!("{r} ≈ {}", fmt_f64(ratio_to_f64(r))),
        },
        Value::Real(x) => format!("≈ {}", fmt_f64(*x)),
        Value::ExactComplex(c) => {
            let s = fmt_exact_complex(c);
            if c.re.is_integer() && c.im.is_integer() {
                s
            } else {
                format!("{s} ≈ {}", fmt_c64(v.to_c64()))
            }
        }
        Value::Complex(c) => format!("≈ {}", fmt_c64(*c)),
    }
}

/// Compact form without the decimal annotation (used inside expressions).
pub fn format_value(v: &Value) -> String {
    match v {
        Value::Rational(r) => r.to_string(),
        Value::Real(x) => fmt_f64(*x),
        Value::ExactComplex(c) => fmt_exact_complex(c),
        Value::Complex(c) => fmt_c64(*c),
    }
}

/// Decimal expansion rounded to DECIMALS places, and whether it's exact.
/// None when the value is too small to show at that precision.
fn decimal(r: &BigRational) -> Option<(String, bool)> {
    let scale = BigRational::from_integer(num_traits::pow(BigInt::from(10), DECIMALS));
    let scaled = r.abs() * scale;
    let exact = scaled.is_integer();
    let digits = scaled.round().to_integer();
    if digits.is_zero() {
        return None;
    }
    let mut s = format!("{:0>width$}", digits.to_string(), width = DECIMALS + 1);
    s.insert(s.len() - DECIMALS, '.');
    let s = s.trim_end_matches('0').trim_end_matches('.');
    Some((if r.is_negative() { format!("-{s}") } else { s.to_string() }, exact))
}

/// About 12 significant digits, trailing zeros trimmed.
pub fn fmt_f64(x: f64) -> String {
    if x == 0.0 {
        return "0".into();
    }
    let a = x.abs();
    if !(1e-6..1e15).contains(&a) {
        let s = format!("{x:.11e}");
        let (mantissa, exp) = s.split_once('e').expect("scientific format");
        let mantissa = mantissa.trim_end_matches('0').trim_end_matches('.');
        return format!("{mantissa}e{exp}");
    }
    let decimals = (11 - a.log10().floor() as i32).clamp(0, 17) as usize;
    let s = format!("{x:.decimals$}");
    let s = if s.contains('.') { s.trim_end_matches('0').trim_end_matches('.') } else { &s };
    if s == "-0" { "0".into() } else { s.to_string() }
}

fn join_complex(re: Option<String>, im_negative: bool, im: Option<String>) -> String {
    match (re, im) {
        (Some(r), Some(i)) => format!("{r} {} {i}", if im_negative { '-' } else { '+' }),
        (None, Some(i)) if im_negative => format!("-{i}"),
        (None, Some(i)) => i,
        (Some(r), None) => r,
        (None, None) => "0".into(),
    }
}

pub fn fmt_c64(c: Complex64) -> String {
    // Hide float noise like 1e-17 next to a much larger component.
    let scale = c.re.abs().max(c.im.abs());
    let clean = |v: f64| if v.abs() < scale * 1e-12 { 0.0 } else { v };
    let (re, im) = (clean(c.re), clean(c.im));
    let im_str = fmt_f64(im.abs());
    let im_term = if im_str == "1" { "i".to_string() } else { format!("{im_str}i") };
    join_complex((re != 0.0).then(|| fmt_f64(re)), im < 0.0, (im != 0.0).then_some(im_term))
}

fn fmt_exact_complex(c: &Complex<BigRational>) -> String {
    let mag = c.im.abs();
    let im_term = if mag.is_one() {
        "i".to_string()
    } else if mag.is_integer() {
        format!("{mag}i")
    } else {
        format!("({mag})i")
    };
    join_complex(
        (!c.re.is_zero()).then(|| c.re.to_string()),
        c.im.is_negative(),
        (!c.im.is_zero()).then_some(im_term),
    )
}
