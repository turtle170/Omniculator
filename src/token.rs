use num_bigint::BigInt;
use num_rational::BigRational;

use crate::error::ParseError;

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
            let start = i;
            let mut name = String::new();
            while let Some(&(_, c)) = chars.get(i) {
                if !(c.is_alphanumeric() || c == '_') {
                    break;
                }
                name.push(c);
                i += 1;
            }
            if name == "π" {
                name = "pi".into();
            }
            tokens.push(Token { kind: TokenKind::Ident(name), start: offset(start), end: offset(i) });
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
