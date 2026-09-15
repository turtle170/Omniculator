//! Names of built-in functions and constants. Shared by the parser (to decide
//! call vs. implicit multiplication) and by typo suggestions.

pub const FUNCTIONS: &[&str] = &[
    "sqrt", "cbrt", "root", "abs", "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh",
    "tanh", "ln", "log", "log2", "exp", "floor", "ceil", "round", "re", "im", "conj", "arg",
];

pub const CONSTANTS: &[&str] = &["pi", "e", "i", "tau", "phi"];

pub fn is_function(name: &str) -> bool {
    FUNCTIONS.contains(&name)
}

pub fn is_constant(name: &str) -> bool {
    CONSTANTS.contains(&name)
}
