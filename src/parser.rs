//! Hand-rolled Pratt parser.
//!
//! Precedence (loosest → tightest): `+ -`, `* /` and implicit multiplication,
//! unary `-`, `^` (right-assoc), postfix `!`. So `-2^2 = -4`, `2^3^2 = 512`.

use crate::ast::Expr;
use crate::builtins::{is_constant, is_function};
use crate::error::ParseError;
use crate::token::{describe, tokenize, Token, TokenKind};
use crate::value::Value;

const ADD_BP: (u8, u8) = (10, 11);
const MUL_BP: (u8, u8) = (20, 21);
const PREFIX_BP: u8 = 25;
const POW_BP: (u8, u8) = (31, 30);
const POSTFIX_BP: u8 = 40;

#[derive(Debug, Clone, PartialEq)]
pub struct Equation {
    pub lhs: Expr,
    pub rhs: Expr,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Input {
    Expr(Expr),
    /// `;`-separated equations, e.g. `x+y=60; x-y=30;`
    Equations(Vec<Equation>),
}

/// Parse a plain expression (no `=` or `;`).
pub fn parse(input: &str) -> Result<Expr, ParseError> {
    parse_tokens(&tokenize(input)?, input.len())
}

/// Parse either a plain expression or a system of equations.
pub fn parse_input(input: &str) -> Result<Input, ParseError> {
    let tokens = tokenize(input)?;
    if !tokens.iter().any(|t| matches!(t.kind, TokenKind::Eq | TokenKind::Semi)) {
        return parse_tokens(&tokens, input.len()).map(Input::Expr);
    }
    let equations = tokens
        .split(|t| t.kind == TokenKind::Semi)
        .filter(|chunk| !chunk.is_empty())
        .map(parse_equation)
        .collect::<Result<Vec<_>, _>>()?;
    if equations.is_empty() {
        return Err(ParseError::new("no equations given", 0, input.len()));
    }
    Ok(Input::Equations(equations))
}

fn parse_equation(chunk: &[Token]) -> Result<Equation, ParseError> {
    let (start, end) = (chunk[0].start, chunk[chunk.len() - 1].end);
    let eq_positions: Vec<usize> =
        chunk.iter().enumerate().filter(|(_, t)| t.kind == TokenKind::Eq).map(|(i, _)| i).collect();
    let &[at] = eq_positions.as_slice() else {
        return Err(match eq_positions.get(1) {
            Some(&i) => {
                ParseError::new("an equation can only have one '='", chunk[i].start, chunk[i].end)
            }
            None => ParseError::new("expected '=' in equation", start, end),
        });
    };
    let eq = &chunk[at];
    let (lhs, rhs) = (&chunk[..at], &chunk[at + 1..]);
    if lhs.is_empty() {
        return Err(ParseError::new("missing left-hand side before '='", eq.start, eq.end));
    }
    if rhs.is_empty() {
        return Err(ParseError::new("missing right-hand side after '='", eq.start, eq.end));
    }
    Ok(Equation { lhs: parse_tokens(lhs, eq.start)?, rhs: parse_tokens(rhs, end)? })
}

/// Parse a whole token slice as one expression; `end` is the byte offset
/// reported for "unexpected end" errors.
fn parse_tokens(tokens: &[Token], end: usize) -> Result<Expr, ParseError> {
    if tokens.is_empty() {
        return Err(ParseError::new("empty input", 0, 0));
    }
    let mut p = Parser { tokens, pos: 0, end };
    let expr = p.expr(0)?;
    if let Some(t) = p.peek() {
        return Err(p.unexpected(t));
    }
    Ok(expr)
}

struct Parser<'a> {
    tokens: &'a [Token],
    pos: usize,
    end: usize,
}

impl<'a> Parser<'a> {
    fn peek(&self) -> Option<&'a Token> {
        self.tokens.get(self.pos)
    }

    fn next(&mut self) -> Option<&'a Token> {
        let t = self.tokens.get(self.pos);
        if t.is_some() {
            self.pos += 1;
        }
        t
    }

    fn peek_is(&self, kind: &TokenKind) -> bool {
        self.peek().is_some_and(|t| t.kind == *kind)
    }

    fn unexpected(&self, t: &Token) -> ParseError {
        let message = match &t.kind {
            TokenKind::RParen => "unmatched ')'".to_string(),
            k => format!("unexpected {}", describe(k)),
        };
        ParseError::new(message, t.start, t.end)
    }

    fn expr(&mut self, min_bp: u8) -> Result<Expr, ParseError> {
        let mut lhs = self.prefix()?;
        while let Some(t) = self.peek() {
            match &t.kind {
                TokenKind::Bang => {
                    if POSTFIX_BP < min_bp {
                        break;
                    }
                    self.pos += 1;
                    lhs = Expr::Factorial(Box::new(lhs));
                }
                TokenKind::Plus | TokenKind::Minus | TokenKind::Star | TokenKind::Slash
                | TokenKind::Caret => {
                    let (l, r) = match t.kind {
                        TokenKind::Plus | TokenKind::Minus => ADD_BP,
                        TokenKind::Star | TokenKind::Slash => MUL_BP,
                        _ => POW_BP,
                    };
                    if l < min_bp {
                        break;
                    }
                    self.pos += 1;
                    let (a, b) = (Box::new(lhs), Box::new(self.expr(r)?));
                    lhs = match t.kind {
                        TokenKind::Plus => Expr::Add(a, b),
                        TokenKind::Minus => Expr::Sub(a, b),
                        TokenKind::Star => Expr::Mul(a, b),
                        TokenKind::Slash => Expr::Div(a, b),
                        _ => Expr::Pow(a, b),
                    };
                }
                // Implicit multiplication: `2x`, `3(x+1)`, `(a)(b)`.
                TokenKind::Num(_) | TokenKind::Ident(_) | TokenKind::LParen => {
                    if MUL_BP.0 < min_bp {
                        break;
                    }
                    if let TokenKind::Num(_) = t.kind {
                        return Err(ParseError::new(
                            format!("missing operator before {}", describe(&t.kind)),
                            t.start,
                            t.end,
                        ));
                    }
                    let rhs = self.expr(MUL_BP.1)?;
                    lhs = Expr::Mul(Box::new(lhs), Box::new(rhs));
                }
                _ => break,
            }
        }
        Ok(lhs)
    }

    fn prefix(&mut self) -> Result<Expr, ParseError> {
        let Some(t) = self.next() else {
            return Err(ParseError::new(
                "expected an expression, found end of input",
                self.end,
                self.end,
            ));
        };
        match &t.kind {
            TokenKind::Num(v) => Ok(Expr::Num(Value::Rational(v.clone()))),
            TokenKind::Minus => Ok(Expr::Neg(Box::new(self.expr(PREFIX_BP)?))),
            TokenKind::Plus => self.expr(PREFIX_BP),
            TokenKind::LParen => {
                let e = self.expr(0)?;
                self.close_paren(t)?;
                Ok(e)
            }
            TokenKind::Ident(name) => self.ident(name, t),
            _ => Err(self.unexpected(t)),
        }
    }

    fn close_paren(&mut self, open: &Token) -> Result<(), ParseError> {
        match self.next() {
            Some(Token { kind: TokenKind::RParen, .. }) => Ok(()),
            Some(t) => Err(self.unexpected(t)),
            None => Err(ParseError::new("unclosed '('", open.start, open.end)),
        }
    }

    fn ident(&mut self, name: &str, t: &Token) -> Result<Expr, ParseError> {
        let paren = self.peek_is(&TokenKind::LParen);
        if is_function(name) && !paren {
            return Err(ParseError::new(
                format!("function '{name}' needs parentheses, e.g. {name}(x)"),
                t.start,
                t.end,
            ));
        }
        // `name(` is a call unless `name` is a constant or single-letter
        // variable, in which case it's implicit multiplication (`e(2)`, `x(y+1)`).
        // Unknown multi-letter names become calls so eval can suggest a fix.
        let is_call = paren && (is_function(name) || !(is_constant(name) || name.chars().count() == 1));
        if !is_call {
            return Ok(Expr::Var(name.to_string()));
        }
        let open = self.next().expect("peeked '('");
        let mut args = Vec::new();
        if self.peek_is(&TokenKind::RParen) {
            self.pos += 1;
            return Ok(Expr::Call(name.to_string(), args));
        }
        loop {
            args.push(self.expr(0)?);
            match self.next() {
                Some(Token { kind: TokenKind::Comma, .. }) => {}
                Some(Token { kind: TokenKind::RParen, .. }) => break,
                Some(tk) => return Err(self.unexpected(tk)),
                None => return Err(ParseError::new("unclosed '('", open.start, open.end)),
            }
        }
        Ok(Expr::Call(name.to_string(), args))
    }
}
