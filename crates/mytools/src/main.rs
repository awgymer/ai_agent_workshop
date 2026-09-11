//! `mytools`: argument dispatch. Each subcommand lives in `cmd/<name>.rs` and owns its
//! flags; this file only routes to it and turns its error into a message and exit code
//! (SPEC §7).

use std::process::ExitCode;

use mytools_core::error::Kind;

mod cmd;

const USAGE: &str = "usage: mytools <sort|merge|intersect|window|subtract> [options]";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let Some(command) = args.first().map(String::as_str) else {
        eprintln!("{USAGE}");
        return ExitCode::from(2);
    };
    let rest = &args[1..];
    let result = match command {
        "--version" => {
            println!("mytools {}", env!("CARGO_PKG_VERSION"));
            return ExitCode::SUCCESS;
        }
        "-h" | "--help" => {
            println!("{USAGE}");
            return ExitCode::SUCCESS;
        }
        "sort" => cmd::sort::run(rest),
        "merge" => cmd::merge::run(rest),
        "intersect" => cmd::intersect::run(rest),
        "window" => cmd::window::run(rest),
        "subtract" => cmd::subtract::run(rest),
        other => {
            eprintln!("mytools: unknown subcommand {other}");
            eprintln!("{USAGE}");
            return ExitCode::from(2);
        }
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) if e.kind() == Kind::BrokenPipe => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("mytools {command}: {}", e.message());
            if let Some(usage) = e.usage_line() {
                eprintln!("{usage}");
            }
            ExitCode::from(e.exit_code())
        }
    }
}
