use num_bigint::BigInt;
use num_rational::BigRational;

use crate::builtins::{callables, is_constant, is_function, CONSTANTS, FUNCTIONS};
use crate::error::ParseError;
use crate::suggest::suggest;

#[derive(Debug, Clone, PartialEq)]
pub enum TokenKind {
    /// Numeric literal, parsed exactly (`0.1` is 1/10, not a float).
    Num(BigRational),
    Ident(String),
    Plus,
    Minus,
    Star,
    Slash,
    Caret,
    Bang,
    LParen,
    RParen,
    Comma,
    Eq,
    Semi,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Token {
    pub kind: TokenKind,
    pub start: usize,
    pub end: usize,
}

const MAX_SCI_EXPONENT: i64 = 10_000;

pub fn describe(kind: &TokenKind) -> String {
    match kind {
        TokenKind::Num(v) => format!("number {v}"),
        TokenKind::Ident(n) => format!("'{n}'"),
        TokenKind::Plus => "'+'".into(),
        TokenKind::Minus => "'-'".into(),
        TokenKind::Star => "'*'".into(),
        TokenKind::Slash => "'/'".into(),
        TokenKind::Caret => "'^'".into(),
        TokenKind::Bang => "'!'".into(),
        TokenKind::LParen => "'('".into(),
        TokenKind::RParen => "')'".into(),
        TokenKind::Comma => "','".into(),
        TokenKind::Eq => "'='".into(),
        TokenKind::Semi => "';'".into(),
    }
}

/// Split a run of letters into names: `xy` → `x`,`y`; `2pix` → `pi`,`x`;
/// `x2y` → `x2`,`y` (trailing digits are subscripts). Known names, names with
/// `_`, and likely typos of a function call (`sqr(`) are kept whole.
fn split_identifier(name: &str, paren: bool) -> Vec<(usize, String)> {
    let whole = || vec![(0, if name == "π" { "pi".to_string() } else { name.to_string() })];
    if name == "π" || is_function(name) || is_constant(name) || name.contains('_') {
        return whole();
    }
    if paren && suggest(name, &callables()).is_some() {
        return whole();
    }
    // `sinx` is almost certainly a missing-parentheses mistake; keep it whole
    // so the error can suggest `sin`.
    if FUNCTIONS.iter().any(|f| f.len() >= 3 && name.starts_with(f)) {
        return whole();
    }
    // `hello(3)`: a long unknown name directly before '(' is meant as a call,
    // not h·e·l·l·o·3 — unless it ends in a real function, like `xsin(x)`.
    if paren
        && name.chars().count() >= 4
        && name.chars().all(|c| c.is_ascii_alphabetic())
        && !(1..name.len()).any(|k| is_function(&name[k..]))
    {
        return whole();
    }
    let mut out = Vec::new();
    let mut pos = 0;
    while pos < name.len() {
        let rest = &name[pos..];
        if paren && is_function(rest) {
            out.push((pos, rest.to_string()));
            break;
        }
        if rest.starts_with('π') {
            out.push((pos, "pi".to_string()));
            pos += 'π'.len_utf8();
            continue;
        }
        if let Some(c) = CONSTANTS.iter().filter(|c| rest.starts_with(*c)).max_by_key(|c| c.len()) {
            out.push((pos, c.to_string()));
            pos += c.len();
            continue;
        }
        let first = rest.chars().next().expect("non-empty").len_utf8();
        let digits: usize = rest[first..].chars().take_while(char::is_ascii_digit).count();
        out.push((pos, rest[..first + digits].to_string()));
        pos += first + digits;
    }
    out
}

pub fn tokenize(input: &str) -> Result<Vec<Token>, ParseError> {
    let chars: Vec<(usize, char)> = input.char_indices().collect();
    let offset = |j: usize| chars.get(j).map_or(input.len(), |&(p, _)| p);
    let is_digit = |j: usize| chars.get(j).is_some_and(|&(_, c)| c.is_ascii_digit());
    let mut tokens = Vec::new();
    let mut i = 0;

    while i < chars.len() {
        let (pos, c) = chars[i];
        if c.is_whitespace() {
            i += 1;
            continue;
        }

        if c.is_ascii_digit() || (c == '.' && is_digit(i + 1)) {
            let start = i;
            let mut digits = String::new();
            while is_digit(i) {
                digits.push(chars[i].1);
                i += 1;
            }
            let mut frac_len = 0i64;
            if chars.get(i).is_some_and(|&(_, c)| c == '.') {
                i += 1;
                while is_digit(i) {
                    digits.push(chars[i].1);
                    frac_len += 1;
                    i += 1;
                }
            }
            // Scientific notation only when a digit follows the `e`, so `2e` stays 2·e.
            let mut exp = 0i64;
            if chars.get(i).is_some_and(|&(_, c)| c == 'e' || c == 'E') {
                let mut j = i + 1;
                let negative = chars.get(j).is_some_and(|&(_, c)| c == '-');
                if chars.get(j).is_some_and(|&(_, c)| c == '+' || c == '-') {
                    j += 1;
                }
                if is_digit(j) {
                    let mut ds = String::new();
                    while is_digit(j) {
                        ds.push(chars[j].1);
                        j += 1;
                    }
                    let e = ds.parse::<i64>().ok().filter(|&e| e <= MAX_SCI_EXPONENT).ok_or_else(
                        || ParseError::new("exponent too large", pos, offset(j)),
                    )?;
                    exp = if negative { -e } else { e };
                    i = j;
                }
            }
            if chars.get(i).is_some_and(|&(_, c)| c == '.') {
                return Err(ParseError::new("unexpected '.' in number", offset(i), offset(i + 1)));
            }
            let mantissa: BigInt = digits.parse().expect("non-empty digit string");
            let scale = exp - frac_len;
            let ten = BigInt::from(10);
            let value = if scale >= 0 {
                BigRational::from_integer(mantissa * num_traits::pow(ten, scale as usize))
            } else {
                BigRational::new(mantissa, num_traits::pow(ten, (-scale) as usize))
            };
            tokens.push(Token { kind: TokenKind::Num(value), start: offset(start), end: offset(i) });
            continue;
        }

        if c.is_alphabetic() || c == '_' {
            let start = offset(i);
            let mut name = String::new();
            while let Some(&(_, c)) = chars.get(i) {
                if !(c.is_alphanumeric() || c == '_') {
                    break;
                }
                name.push(c);
                i += 1;
            }
            let paren = chars[i..].iter().find(|&&(_, c)| !c.is_whitespace()).is_some_and(|&(_, c)| c == '(');
            for (at, seg) in split_identifier(&name, paren) {
                let s = start + at;
                let len = if seg == "pi" && name[at..].starts_with('π') { 'π'.len_utf8() } else { seg.len() };
                tokens.push(Token { kind: TokenKind::Ident(seg), start: s, end: s + len });
            }
            continue;
        }

        let (kind, len) = match c {
            '*' if chars.get(i + 1).is_some_and(|&(_, c)| c == '*') => (TokenKind::Caret, 2),
            '+' => (TokenKind::Plus, 1),
            '-' | '−' => (TokenKind::Minus, 1),
            '*' | '×' | '·' => (TokenKind::Star, 1),
            '/' | '÷' => (TokenKind::Slash, 1),
            '^' => (TokenKind::Caret, 1),
            '!' => (TokenKind::Bang, 1),
            '(' => (TokenKind::LParen, 1),
            ')' => (TokenKind::RParen, 1),
            ',' => (TokenKind::Comma, 1),
            '=' => (TokenKind::Eq, 1),
            ';' => (TokenKind::Semi, 1),
            _ => {
                return Err(ParseError::new(
                    format!("unexpected character '{c}'"),
                    pos,
                    offset(i + 1),
                ))
            }
        };
        tokens.push(Token { kind, start: pos, end: offset(i + len) });
        i += len;
    }
    Ok(tokens)
}
