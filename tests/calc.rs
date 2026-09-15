use omniculator::format::format_result;
use omniculator::suggest::levenshtein;
use omniculator::{calculate, Value};

fn show(input: &str) -> String {
    match calculate(input) {
        Ok(v) => format_result(&v),
        Err(e) => format!("error: {e}"),
    }
}

fn real(input: &str) -> f64 {
    match calculate(input).unwrap() {
        Value::Real(x) => x,
        v => panic!("{input}: expected Real, got {v:?}"),
    }
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

#[test]
fn arithmetic_is_exact() {
    assert_eq!(show("2+3*4"), "14");
    assert_eq!(show("1/3+1/6"), "1/2 = 0.5");
    assert_eq!(show("0.1+0.2"), "3/10 = 0.3");
    assert_eq!(show("1/3"), "1/3 ≈ 0.333333333333");
    assert_eq!(show("2/3"), "2/3 ≈ 0.666666666667");
    assert_eq!(show("-7/2"), "-7/2 = -3.5");
    assert_eq!(show("10!/8!"), "90");
    assert_eq!(show("1e3"), "1000");
    assert_eq!(show("2.5e-1"), "1/4 = 0.25");
}

#[test]
fn precedence() {
    assert_eq!(show("-2^2"), "-4");
    assert_eq!(show("(-2)^2"), "4");
    assert_eq!(show("2^3^2"), "512");
    assert_eq!(show("2^-1"), "1/2 = 0.5");
    assert_eq!(show("5!"), "120");
    assert_eq!(show("3!^2"), "36");
    assert_eq!(show("-3!"), "-6");
    assert_eq!(show("10-4-3"), "3");
    assert_eq!(show("2**3"), "8");
}

#[test]
fn implicit_multiplication() {
    assert_eq!(show("2(3+4)"), "14");
    assert_eq!(show("(1+2)(3+4)"), "21");
    assert!(close(real("2pi"), std::f64::consts::TAU));
    assert!(close(real("2e"), 2.0 * std::f64::consts::E));
}

#[test]
fn roots_and_powers() {
    assert_eq!(show("sqrt(16)"), "4");
    assert_eq!(show("sqrt(4/9)"), "2/3 ≈ 0.666666666667");
    assert_eq!(show("8^(2/3)"), "4");
    assert_eq!(show("(-8)^(1/3)"), "-2");
    assert_eq!(show("root(81, 4)"), "3");
    assert_eq!(show("0^0"), "1");
    assert!(close(real("sqrt(2)^2"), 2.0));
}

#[test]
fn complex_numbers() {
    assert_eq!(show("sqrt(-4)"), "2i");
    assert_eq!(show("i^2"), "-1");
    assert_eq!(show("(1+2i)(3-i)"), "5 + 5i");
    assert_eq!(show("1/i"), "-i");
    assert_eq!(show("abs(3+4i)"), "5");
    assert_eq!(show("conj(1+i)"), "1 - i");
    assert_eq!(show("i/2"), "(1/2)i ≈ 0.5i");
    match calculate("sqrt(-2)").unwrap() {
        Value::Complex(c) => assert!(close(c.im, 2f64.sqrt()) && c.re.abs() < 1e-12),
        v => panic!("expected Complex, got {v:?}"),
    }
    match calculate("ln(-1)").unwrap() {
        Value::Complex(c) => assert!(close(c.im, std::f64::consts::PI)),
        v => panic!("expected Complex, got {v:?}"),
    }
}

#[test]
fn functions() {
    assert_eq!(show("log(1000)"), "3");
    assert_eq!(show("log(8, 2)"), "3");
    assert_eq!(show("log2(1/8)"), "-3");
    assert_eq!(show("sin(0)"), "0");
    assert_eq!(show("floor(7/2)"), "3");
    assert!(close(real("cos(pi)"), -1.0));
    assert!(close(real("ln(e)"), 1.0));
}

#[test]
fn errors() {
    assert_eq!(show("1/0"), "error: division by zero");
    assert!(show("sqr(4)").contains("did you mean 'sqrt'"));
    assert!(show("sqrtt(4)").contains("did you mean 'sqrt'"));
    assert!(show("sinx").contains("did you mean 'sin'"));
    assert_eq!(show("x"), "error: unknown name 'x'");
    assert!(show("2 3").contains("missing operator"));
    assert!(show("2+").contains("end of input"));
    assert!(show("(2").contains("unclosed"));
    assert!(show("2)").contains("unmatched"));
    assert!(show("sin 3").contains("parentheses"));
    assert!(show("x+1=2").contains("unexpected '='"));
    assert!(show("(-1)!").contains("non-negative integers"));
    assert!(show("sqrt(1,2)").contains("takes 1 argument"));
    assert!(show("log(0)").contains("undefined"));
}

#[test]
fn edit_distance() {
    assert_eq!(levenshtein("kitten", "sitting"), 3);
    assert_eq!(levenshtein("", "abc"), 3);
    assert_eq!(levenshtein("sqrt", "sqrt"), 0);
}
