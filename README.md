# Omniculator

A command-line math tool written in Rust. Type in arithmetic, fractions,
complex numbers or systems of linear equations and it computes the answer.

Every parser, evaluator and solver is hand-written and deterministic. No AI
is used at runtime.

## Features

- **Exact arithmetic.** Fractions stay exact (`0.1 + 0.2 = 3/10`), backed by
  big integers so large results don't overflow (`25!`). Results that can't be
  exact (`sqrt(2)`, `sin(1)`) are shown as decimals and marked with `≈`.
- **Complex numbers.** `sqrt(-4) = 2i` and `(1+2i)(3-i) = 5 + 5i`.
- **Scientific functions.** `sqrt cbrt root abs sin cos tan asin acos atan
  sinh cosh tanh ln log log2 exp floor ceil round re im conj arg`, and the
  constants `pi e i tau phi`.
- **Implicit multiplication.** `2x`, `3(x+1)`, `2pi`.
- **Linear systems** of any size, solved exactly with step-by-step output.
  It reports one solution, infinitely many solutions (with dependent
  variables written in terms of the free ones) or no solution.
- **Helpful errors.** Typos get suggestions (`sqr(2)` gives "did you mean
  'sqrt'?") and syntax errors point to where the problem is.

## Usage

```bash
cargo build --release
```

Pass an expression to get a single answer:

```text
$ omniculator "1/3 + 1/6"
1/2 = 0.5

$ omniculator "x+y=60; x-y=30"
Steps:
  1. Write each equation in standard form
       x + y = 60
       x - y = 30
  2. Isolate x: R2 ← R2 - R1
       x + y = 60
       -2y = -30
  3. Isolate y: R2 ← R2 ÷ (-2); R1 ← R1 - R2
       x = 45
       y = 15
x = 45
y = 15
```

Run it with no arguments to start an interactive session (REPL). Type
`quit` to exit.

## Syntax notes

- `^` binds tighter than a leading minus and groups right to left:
  `-2^2 = -4`, `2^3^2 = 512`. `**` also means `^`.
- Put `;` between equations. Each equation has exactly one `=`.
- Multi-letter names are single variables: `xy` is one variable, not `x*y`.
- `2e3` is 2000 (scientific notation), but `2e` is 2 × e.

## Development

```bash
cargo test
```

## Roadmap

1. Calculator core ✅
2. Linear system solver ✅
3. Quadratic solver
4. Simplify / expand
5. Differentiation
6. Partial derivatives and gradients
7. Integration
8. Factoring
