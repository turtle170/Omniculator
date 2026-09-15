use omniculator::linear::{format_solution, Solution};
use omniculator::{run, Error, Outcome, Value};

fn solve(input: &str) -> Solution {
    match run(input).unwrap() {
        Outcome::Linear(s, _) => s,
        o => panic!("{input}: expected a system, got {o:?}"),
    }
}

fn unique(pairs: &[(&str, Value)]) -> Solution {
    Solution::Unique(pairs.iter().map(|(n, v)| (n.to_string(), v.clone())).collect())
}

fn error(input: &str) -> String {
    match run(input) {
        Err(e) => e.to_string(),
        Ok(o) => panic!("{input}: expected an error, got {o:?}"),
    }
}

#[test]
fn two_by_two() {
    assert_eq!(solve("x+y=60;x-y=30;"), unique(&[("x", Value::int(45)), ("y", Value::int(15))]));
}

#[test]
fn single_equation() {
    assert_eq!(solve("2x+3=7"), unique(&[("x", Value::int(2))]));
    assert_eq!(solve("3(x-1) = x + 5"), unique(&[("x", Value::int(4))]));
}

#[test]
fn exact_fractions() {
    assert_eq!(
        solve("x/3 + y/4 = 1; x - y = 1/2"),
        unique(&[("x", Value::rat(27, 14)), ("y", Value::rat(10, 7))])
    );
}

#[test]
fn three_by_three() {
    assert_eq!(
        solve("x+y+z=6; 2y+5z=-4; 2x+5y-z=27"),
        unique(&[("x", Value::int(5)), ("y", Value::int(3)), ("z", Value::int(-2))])
    );
}

#[test]
fn needs_row_swap() {
    assert_eq!(solve("x - x + y = 1; x + y = 3"), unique(&[("x", Value::int(2)), ("y", Value::int(1))]));
}

#[test]
fn no_solution() {
    assert!(matches!(solve("x+y=1; x+y=2"), Solution::None { .. }));
    assert!(matches!(solve("2=3"), Solution::None { .. }));
}

#[test]
fn infinitely_many() {
    match solve("x+y=3; 2x+2y=6") {
        Solution::Infinite { dependent, free } => {
            assert_eq!(free, vec!["y".to_string()]);
            assert_eq!(dependent[0].0, "x");
            assert_eq!(dependent[0].1.to_string(), "3 - y");
        }
        s => panic!("expected infinitely many, got {s:?}"),
    }
}

#[test]
fn identity() {
    assert_eq!(solve("2=2"), Solution::Unique(vec![]));
}

#[test]
fn complex_and_float_coefficients() {
    assert_eq!(solve("i*x = 2"), unique(&[("x", Value::rat(-2, 1).mul(&solve_i()))]));
    match solve("pi*x = pi") {
        Solution::Unique(v) => assert!((v[0].1.as_real_f64().unwrap() - 1.0).abs() < 1e-12),
        s => panic!("{s:?}"),
    }
    // Constant subexpressions are evaluated exactly.
    assert_eq!(solve("sin(0)*x + x = 4"), unique(&[("x", Value::int(4))]));
}

fn solve_i() -> Value {
    omniculator::calculate("i").unwrap()
}

fn text(input: &str) -> String {
    match run(input) {
        Ok(Outcome::Text { answer, .. }) => answer,
        o => panic!("{input}: expected a text answer, got {o:?}"),
    }
}

#[test]
fn nonlinear_goes_to_polynomial_solver() {
    assert_eq!(text("x^2=4"), "x = -2\nx = 2");
    assert!(text("x*y=6").ends_with("x = 6/y"));
    assert!(text("sin(x)=0").contains("x ≈ 3.14159265359\n"));
    assert!(error("sin(x)=y").contains("isn't a polynomial equation"));
}

#[test]
fn equation_parse_errors() {
    assert!(matches!(run("x+y"), Ok(Outcome::Text { .. })));
    assert!(error("1/0 = x").contains("division by zero"));
    let _ = Error::Eval;
    assert!(error("x=1=2").contains("only have one '='"));
    assert!(error("x+1; y=2").contains("expected '='"));
    assert!(error("=5").contains("left-hand side"));
    assert!(error("x=").contains("right-hand side"));
}

#[test]
fn step_trace() {
    let Outcome::Linear(s, steps) = run("x+y=60; x-y=30").unwrap() else { panic!() };
    let text = format_solution(&s, &steps);
    assert!(text.contains("R2 ← R2 - R1"), "{text}");
    assert!(text.contains("R2 ← R2 ÷ (-2)"), "{text}");
    assert!(text.ends_with("x = 45\ny = 15"), "{text}");
}
