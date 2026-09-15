use omniculator::format::format_result;
use omniculator::linear::format_solution;
use omniculator::{run, Outcome};

fn answer(input: &str) -> String {
    match run(input) {
        Ok(Outcome::Value(v)) => format_result(&v),
        Ok(Outcome::Exact(s, _)) => s.to_string(),
        Ok(Outcome::Linear(s, steps)) => format_solution(&s, &steps),
        Ok(Outcome::Text { answer, .. }) => answer,
        Err(e) => format!("error: {e}"),
    }
}

#[test]
fn implicit_products_of_letters() {
    assert_eq!(answer("2xy + xy"), "3xy");
    assert_eq!(answer("x^2y^2 - y^2x^2"), "0");
}

#[test]
fn nonlinear_system() {
    let a = answer("x^4 + x^2y^2 + y^4 = 21; x^2 + xy + y^2 = 7");
    assert_eq!(a, "x = -2, y = -1\nx = -1, y = -2\nx = 1, y = 2\nx = 2, y = 1");
    let circle = answer("x^2 + y^2 = 25; x + y = 7");
    assert_eq!(circle, "x = 3, y = 4\nx = 4, y = 3");
    assert!(answer("x^2 + y^2 = 1; x^2 + y^2 = 4").starts_with("No solution"));
}

#[test]
fn parameters_and_large_systems() {
    assert_eq!(answer("x + y = a; x - y = b"), "Solving for x, y in terms of a, b:\n  y = a/2 - b/2\n  x = a/2 + b/2");
    assert!(answer("x^2 + y^2 = a; x^3 + y^3 = b").contains("degree-6"));
    let chain = answer("x^2+y=2; y^2+z=3; z^2+w=5; w^2+x=7");
    assert_eq!(chain.lines().count(), 16, "{chain}");
    let big = answer("x^4 + y^3z + z^2w^2 = 2; y^4 + z^3w + w^2x^2 = 3; z^4 + w^3x + x^2y^2 = 5; w^4 + x^3y + y^2z^2 = 7");
    assert!(big.lines().count() >= 12, "{big}");
}

#[test]
fn non_polynomial_systems() {
    let a = answer("x^2 + sin(y) = 1; y^2 + cos(x) = 2");
    assert_eq!(a.lines().count(), 4, "{a}");
    assert!(answer("x^(3/2) + y^2 = 4; y^(3/2) + z^2 = 5; z^(3/2) + x^2 = 6").contains("x ≈ 1.83505623583"));
    let five = answer("x + 2y + 2z + 2w + 2v = 1; x^2 + 2y^2 + 2z^2 + 2w^2 + 2v^2 = x; 2xy + 2yz + 2zw + 2wv = y; 2xz + 2yw + 2zv + y^2 = z; 2xw + 2yv + 2yz = w");
    assert!(!five.contains("e-"), "{five}");
}

#[test]
fn stress_regressions() {
    assert!(answer("tan(pi/2)").starts_with("error"));
    assert_eq!(answer("e^(i*pi)"), "-1");
    assert_eq!(answer("simplify(sin(x)^2 + cos(x)^2)"), "1");
    assert_eq!(answer("integrate(1/x, x, 1, e)"), "∫ from 1 to e of 1/x dx = 1");
    assert_eq!(answer("sqrt(x) = 3"), "x = 9");
}

#[test]
fn quadratics_and_polynomials() {
    assert_eq!(answer("x^2 - 5x + 6 = 0"), "x = 2\nx = 3");
    assert!(answer("x^2 = 2").contains("x = sqrt(2)"));
    assert!(answer("x^2 + x + 1 = 0").contains("(complex)"));
    assert!(answer("x^3 = 8").starts_with("x = 2"));
    assert!(answer("x^3 = 2").starts_with("x = cbrt(2)"));
    assert_eq!(answer("x^3 - 6x^2 + 11x - 6 = 0"), "x = 1\nx = 2\nx = 3");
    assert_eq!(answer("1/x = 2"), "x = 1/2");
    // (x-1)/(x-1) cancels to 1, and x = 1 is outside the domain anyway.
    assert!(answer("x/(x-1) = 1/(x-1)").starts_with("No solution"));
}

#[test]
fn simplify_and_expand() {
    assert_eq!(answer("2x+3x"), "5x");
    assert_eq!(answer("(x+1)(x+2)"), "x^2 + 3x + 2");
    assert_eq!(answer("simplify((x^2-1)/(x-1))"), "x + 1");
    assert_eq!(answer("simplify(x^2 + 2x + 1)"), "(x + 1)^2");
    assert_eq!(answer("expand((x+1)^3)"), "x^3 + 3x^2 + 3x + 1");
    assert_eq!(answer("sqrt(8)"), "2sqrt(2)");
    assert_eq!(answer("sin(pi/6)"), "1/2 = 0.5");
}

#[test]
fn factoring() {
    assert_eq!(answer("factor(x^3 - 6x^2 + 11x - 6)"), "(x - 3)(x - 2)(x - 1)");
    assert_eq!(answer("factor(x^4 + x^2y^2 + y^4)"), "(x^2 - xy + y^2)(x^2 + xy + y^2)");
    assert_eq!(answer("factor(2x^2 - 2)"), "2(x - 1)(x + 1)");
    assert_eq!(answer("factor(x^2 + 1)"), "x^2 + 1");
    assert_eq!(answer("factor(x^2 - y^2)"), "(x - y)(x + y)");
}

#[test]
fn differentiation() {
    assert!(answer("diff(x^3, x)").ends_with("= 3x^2"));
    assert!(answer("diff(e^(2x), x)").ends_with("= 2e^(2x)"));
    assert!(answer("diff(sin(x)^2, x)").ends_with("= 2cos(x)·sin(x)"));
    assert!(answer("diff(ln(x), x)").ends_with("= 1/x"));
    assert!(answer("diff(x^4, x, 2)").ends_with("= 12x^2"));
    assert_eq!(answer("grad(x^2y + y^3)"), "∂f/∂x = 2xy\n∂f/∂y = x^2 + 3y^2\n∇f = (2xy, x^2 + 3y^2)");
    assert_eq!(answer("hessian(x^2y)"), "variables: x, y\n[2y, 2x]\n[2x, 0]");
    assert_eq!(answer("jacobian(xy, x + y)"), "variables: x, y\n[y, x]\n[1, 1]");
}

#[test]
fn integration() {
    assert_eq!(answer("integrate(x e^x, x)"), "∫ x·e^x dx = x·e^x - e^x + C");
    assert_eq!(answer("integrate(2x cos(x^2), x)"), "∫ 2x·cos(x^2) dx = sin(x^2) + C");
    assert_eq!(answer("integrate(1/(x^2+1), x)"), "∫ 1/(x^2 + 1) dx = atan(x) + C");
    assert_eq!(answer("integrate(1/x, x)"), "∫ 1/x dx = ln(abs(x)) + C");
    assert_eq!(answer("integrate(ln(x), x)"), "∫ ln(x) dx = x·ln(x) - x + C");
    assert_eq!(answer("integrate(x^2, x, 0, 3)"), "∫ from 0 to 3 of x^2 dx = 9");
    assert_eq!(answer("integrate(sin(x), x, 0, pi)"), "∫ from 0 to pi of sin(x) dx = 2");
    assert!(answer("integrate(e^(x^2), x)").starts_with("error: couldn't integrate"));
    assert!(answer("integrate(e^(x^2), x, 0, 1)").contains("≈ 1.46265174"));
}
