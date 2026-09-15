//! Human-readable (and re-parseable) display of symbolic expressions:
//! `x^2 + 3x + 2`, `sqrt(2)/2`, `x/(2y)`, `e^x·y`.

use std::cmp::Ordering;
use std::fmt;

use num_rational::BigRational as Q;
use num_traits::{One, Signed, Zero};

use crate::poly::lex_cmp;
use crate::sym::{Atom, Konst, Mono, Sym};

fn var_degree(m: &Mono) -> Q {
    m.iter().filter(|(a, _)| matches!(a, Atom::Var(_))).map(|(_, e)| e.clone()).sum()
}

/// Highest total degree first, then lexicographic, constants last.
pub fn display_order(a: &(&Mono, &Q), b: &(&Mono, &Q)) -> Ordering {
    var_degree(b.0).cmp(&var_degree(a.0)).then_with(|| lex_cmp(b.0, a.0)).then_with(|| a.0.cmp(b.0))
}

fn rank(a: &Atom) -> u8 {
    match a {
        Atom::Num(_) => 0,
        Atom::I => 1,
        Atom::Const(_) => 2,
        Atom::Var(_) => 3,
        Atom::Func(..) => 4,
        Atom::Group(_) => 5,
        Atom::Pow(..) => 6,
    }
}

/// True if `s` can be followed by `^` or juxtaposed without parentheses.
fn is_atomic(s: &str) -> bool {
    // A plain number or a single name (not `2x`, which needs parentheses).
    if s.chars().all(|c| c.is_ascii_digit()) && !s.is_empty() {
        return true;
    }
    if s.starts_with(|c: char| c.is_alphabetic() || c == '_') && s.chars().all(|c| c.is_alphanumeric() || c == '_') {
        return true;
    }
    // A single call like `sin(x + 1)` or an already-parenthesized group.
    let Some(open) = s.find('(') else { return false };
    if !s.ends_with(')') || !s[..open].chars().all(char::is_alphanumeric) {
        return false;
    }
    let mut depth = 0;
    for (i, c) in s.char_indices().skip(open) {
        match c {
            '(' => depth += 1,
            ')' => {
                depth -= 1;
                if depth == 0 && i != s.len() - 1 {
                    return false;
                }
            }
            _ => {}
        }
    }
    true
}

fn paren(s: String) -> String {
    if is_atomic(&s) {
        s
    } else {
        format!("({s})")
    }
}

fn exponent_str(e: &Q) -> String {
    if e.is_integer() && !e.is_negative() {
        e.to_string()
    } else {
        format!("({e})")
    }
}

/// The atom's own text, without any exponent.
fn atom_str(a: &Atom) -> String {
    match a {
        Atom::Var(n) => n.clone(),
        Atom::Const(Konst::Pi) => "pi".into(),
        Atom::Const(Konst::E) => "e".into(),
        Atom::I => "i".into(),
        Atom::Num(n) => n.to_string(),
        Atom::Func(name, args) if name == "factorial" => format!("{}!", paren(args[0].to_string())),
        Atom::Func(name, args) => {
            let args: Vec<String> = args.iter().map(ToString::to_string).collect();
            format!("{name}({})", args.join(", "))
        }
        Atom::Group(b) => format!("({b})"),
        Atom::Pow(b, x) => format!("{}^{}", paren(b.to_string()), paren(x.to_string())),
    }
}

/// Inner text for `sqrt(...)` / `cbrt(...)`.
fn inner_str(a: &Atom) -> String {
    match a {
        Atom::Group(b) => b.to_string(),
        other => atom_str(other),
    }
}

fn power_str(a: &Atom, e: &Q) -> String {
    if e.is_one() {
        return atom_str(a);
    }
    if *e == Q::new(1.into(), 2.into()) {
        return format!("sqrt({})", inner_str(a));
    }
    if *e == Q::new(1.into(), 3.into()) {
        return format!("cbrt({})", inner_str(a));
    }
    let base = match a {
        Atom::Group(_) => atom_str(a),
        _ => paren(atom_str(a)),
    };
    format!("{base}^{}", exponent_str(e))
}

fn juxtapose(parts: &[String]) -> String {
    let mut out = String::new();
    for p in parts {
        let simple = p.chars().all(|c| c.is_alphanumeric() || c == '_');
        let after_name = out.ends_with(|c: char| c.is_alphabetic() || c == ')') || out.contains('^');
        // `e^x` followed by `y` would read as `e^(xy)`; `x·sin(x)` reads better than `xsin(x)`.
        let ambiguous = out.contains('^') && out.ends_with(char::is_alphabetic);
        if !out.is_empty() && (ambiguous || (!simple && after_name)) {
            out.push('·');
        }
        out.push_str(p);
    }
    out
}

/// A term with a positive coefficient `c`.
pub fn fmt_term(c: &Q, m: &Mono) -> String {
    let mut atoms: Vec<(&Atom, &Q)> = m.iter().collect();
    atoms.sort_by_key(|(a, _)| rank(a));
    let mut num = Vec::new();
    let mut den = Vec::new();
    for (a, e) in atoms {
        if e.is_positive() {
            num.push(power_str(a, e));
        } else {
            den.push(power_str(a, &-e));
        }
    }
    let mut top = Vec::new();
    if !c.numer().is_one() || num.is_empty() {
        top.push(c.numer().to_string());
    }
    top.extend(num);
    let top = juxtapose(&top);
    let mut bottom = Vec::new();
    if !c.denom().is_one() {
        bottom.push(c.denom().to_string());
    }
    bottom.extend(den);
    match bottom.len() {
        0 => top,
        1 => format!("{top}/{}", bottom[0]),
        _ => format!("{top}/({})", juxtapose(&bottom)),
    }
}

impl fmt::Display for Sym {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.is_zero() {
            return f.write_str("0");
        }
        let mut terms: Vec<(&Mono, &Q)> = self.terms.iter().collect();
        terms.sort_by(display_order);
        for (i, (m, c)) in terms.iter().enumerate() {
            let body = fmt_term(&c.abs(), m);
            match (i, c.is_negative()) {
                (0, true) => write!(f, "-{body}")?,
                (0, false) => f.write_str(&body)?,
                (_, true) => write!(f, " - {body}")?,
                (_, false) => write!(f, " + {body}")?,
            }
        }
        Ok(())
    }
}

#[allow(dead_code)]
fn _assert_zero_is_constant() {
    debug_assert!(Sym::zero().as_constant().is_some_and(|c| c.is_zero()));
}
