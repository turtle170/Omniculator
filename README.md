# Omniculator

A command-line math tool written in Rust. Type in arithmetic, fractions,
complex numbers, algebra, equations or calculus and it computes the answer.

Every parser, evaluator and solver is hand-written and deterministic. No AI
is used at runtime. When there's no algorithm for something (general
integration, for example), Omniculator says so instead of guessing.

## Install (Windows)

```powershell
irm https://raw.githubusercontent.com/turtle170/Omniculator/main/install.ps1 | iex
```

This downloads the latest release into `%LOCALAPPDATA%\Programs\Omniculator`
and adds that folder to your PATH. To build from source instead, run
`cargo build --release`.

## What it can do

| Input | Result |
|---|---|
| `1/3 + 1/6` | `1/2 = 0.5` |
| `sqrt(8)` | `2sqrt(2) ≈ 2.82842712475` |
| `(1+2i)(3-i)` | `5 + 5i` |
| `2x + 3x` | `5x` |
| `x+y=60; x-y=30` | `x = 45`, `y = 15` (with elimination steps) |
| `x^2 - 5x + 6 = 0` | `x = 2`, `x = 3` (discriminant and quadratic formula shown) |
| `x^4 + x^2y^2 + y^4 = 21; x^2 + xy + y^2 = 7` | all four solutions, exactly |
| `simplify((x^2-1)/(x-1))` | `x + 1` |
| `expand((x+1)^3)` | `x^3 + 3x^2 + 3x + 1` |
| `factor(x^4 + x^2y^2 + y^4)` | `(x^2 - xy + y^2)(x^2 + xy + y^2)` |
| `diff(x^2 sin(x), x)` | `x^2·cos(x) + 2x·sin(x)` (rules shown) |
| `grad(x^2y + y^3)` | `∇f = (2xy, x^2 + 3y^2)` |
| `integrate(x e^x, x)` | `x·e^x - e^x + C` (by parts, verified) |
| `integrate(sin(x), x, 0, pi)` | `2` |

- **Exact arithmetic.** Fractions, radicals (`2sqrt(2)`), π and i stay exact.
  Big integers are used, so large results don't overflow. Anything that has
  to be approximated is marked with `≈`.
- **Equations.** Linear systems of any size, single polynomial equations of
  any degree, and systems of polynomial equations (solved with Gröbner
  bases). Roots are exact wherever possible and numeric (marked `≈`) only
  for polynomials of degree 3 or more that have no exact factors. Complex
  roots are labeled as complex.
- **Algebra.** `simplify`, `expand` and `factor`. Factoring is over the
  rationals and uses the rational root theorem and Kronecker's method.
- **Calculus.** `diff(f, x)` and `diff(f, x, n)`, which also work as partial
  derivatives. `grad`, `jacobian` and `hessian`. `integrate(f, x)` and
  `integrate(f, x, a, b)` use standard forms, substitution, integration by
  parts and partial fractions. Every antiderivative is checked by
  differentiating it. If no antiderivative is found, a definite integral
  falls back to a numeric value, and the output says so.
- **Functions and constants.** `sqrt cbrt root abs sin cos tan asin acos
  atan sinh cosh tanh ln log log2 exp floor ceil round re im conj arg`, and
  `pi e i tau phi`.
- **Helpful errors.** Typos get suggestions (`sqr(2)` gives "did you mean
  'sqrt'?") and syntax errors point to where the problem is.

## Usage

Pass the input as an argument for a single answer, or run `omniculator`
with no arguments to start an interactive session (REPL). Type `quit` to
exit.

```text
$ omniculator "x^2 + y^2 = 25; x + y = 7"
Steps:
  1. Move everything to one side
       x^2 + y^2 - 25 = 0
       x + y - 7 = 0
  2. Compute a Gröbner basis (lex order x > y) ...
  ...
x = 3, y = 4
x = 4, y = 3
```

## Syntax notes

- Letters next to each other multiply: `xy` is `x*y`, `2pix` is `2*pi*x`.
  Digits after a letter stay part of its name (`x1`), and so does `_`
  (`rate_1`).
- `^` binds tighter than a leading minus and groups right to left:
  `-2^2 = -4`, `2^3^2 = 512`. `**` also means `^`.
- Put `;` between equations. Each equation has exactly one `=`.
- `2e3` is 2000 (scientific notation), but `2e` is 2 × e.

## Development

```bash
cargo test
```
