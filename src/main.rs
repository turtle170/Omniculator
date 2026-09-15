use std::io::{self, BufRead, Write};

use omniculator::format::format_result;
use omniculator::linear::format_solution;
use omniculator::{run, Error, Outcome};

fn render(input: &str) -> Result<String, Error> {
    Ok(match run(input)? {
        Outcome::Value(v) => format_result(&v),
        Outcome::Linear(solution, steps) => format_solution(&solution, &steps),
    })
}

fn report(input: &str, e: &Error) {
    if let Error::Parse(p) = e {
        let start = p.start.min(input.len());
        let end = p.end.clamp(start, input.len());
        let col = input[..start].chars().count();
        let width = input[start..end].chars().count().max(1);
        eprintln!("  {input}");
        eprintln!("  {}{}", " ".repeat(col), "^".repeat(width));
    }
    eprintln!("error: {e}");
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();

    // Single-shot: `omniculator "1/3 + 1/6"` or `omniculator "x+y=60; x-y=30"`
    if !args.is_empty() {
        let input = args.join(" ");
        match render(&input) {
            Ok(s) => println!("{s}"),
            Err(e) => {
                report(&input, &e);
                std::process::exit(1);
            }
        }
        return;
    }

    // REPL
    println!("Omniculator. Type an expression or equations, or 'quit' to exit.");
    let stdin = io::stdin();
    loop {
        print!("> ");
        io::stdout().flush().ok();
        let mut line = String::new();
        if stdin.lock().read_line(&mut line).unwrap_or(0) == 0 {
            break;
        }
        let input = line.trim();
        if input.is_empty() {
            continue;
        }
        if matches!(input, "quit" | "exit") {
            break;
        }
        match render(input) {
            Ok(s) => println!("{s}"),
            Err(e) => report(input, &e),
        }
    }
}
