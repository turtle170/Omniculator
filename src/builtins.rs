//! Names of built-in functions, commands and constants. Shared by the parser
//! (to decide call vs. implicit multiplication) and by typo suggestions.

pub const FUNCTIONS: &[&str] = &[
    "sqrt", "cbrt", "root", "abs", "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh",
    "tanh", "ln", "log", "log2", "exp", "floor", "ceil", "round", "re", "im", "conj", "arg",
];

/// Symbolic commands; only valid as the whole input, e.g. `diff(x^2, x)`.
pub const COMMANDS: &[&str] = &[
    "simplify", "expand", "factor", "diff", "derivative", "integrate", "int", "grad", "gradient",
    "jacobian", "hessian",
];

pub const CONSTANTS: &[&str] = &["pi", "e", "i", "tau", "phi"];

pub fn is_function(name: &str) -> bool {
    FUNCTIONS.contains(&name) || COMMANDS.contains(&name)
}

pub fn is_command(name: &str) -> bool {
    COMMANDS.contains(&name)
}

pub fn is_constant(name: &str) -> bool {
    CONSTANTS.contains(&name)
}

/// Everything callable, for typo suggestions.
pub fn callables() -> Vec<&'static str> {
    FUNCTIONS.iter().chain(COMMANDS).copied().collect()
}
