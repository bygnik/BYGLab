//! Thin native CLI around `byglab-core`. Reads a JSON [`PipeCaseConfig`]
//! file, runs it, prints the resulting [`PipeCaseResult`] as JSON to
//! stdout. Useful for fast local iteration without a browser, and as a
//! skeleton for future OpenWAM reference-data diffing.

use byglab_core::{run_pipe_case, PipeCaseConfig};
use std::{env, fs, process};

fn main() {
    let mut args = env::args();
    let program = args.next().unwrap_or_else(|| "byglab-cli".into());
    let Some(path) = args.next() else {
        eprintln!("usage: {program} <case.json>");
        process::exit(2);
    };

    let text = fs::read_to_string(&path).unwrap_or_else(|e| {
        eprintln!("error: could not read {path}: {e}");
        process::exit(1);
    });

    let config: PipeCaseConfig = serde_json::from_str(&text).unwrap_or_else(|e| {
        eprintln!("error: invalid case config in {path}: {e}");
        process::exit(1);
    });

    let result = run_pipe_case(&config);
    let json = serde_json::to_string_pretty(&result).expect("PipeCaseResult is always serializable");
    println!("{json}");
}
