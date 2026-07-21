//! Thin native CLI around `byglab-core`. Reads a JSON `SimConfig` file,
//! runs the (placeholder) solver, prints the JSON `SimResult` to stdout.
//! Used for fast local iteration without a browser, and as the skeleton
//! for the future OpenWAM reference-data diffing harness.

use byglab_core::{run, SimConfig};
use std::{env, fs, process};

fn main() {
    let mut args = env::args();
    let program = args.next().unwrap_or_else(|| "byglab-cli".into());
    let Some(path) = args.next() else {
        eprintln!("usage: {program} <config.json>");
        process::exit(2);
    };

    let text = fs::read_to_string(&path).unwrap_or_else(|e| {
        eprintln!("error: could not read {path}: {e}");
        process::exit(1);
    });

    let config: SimConfig = serde_json::from_str(&text).unwrap_or_else(|e| {
        eprintln!("error: invalid config in {path}: {e}");
        process::exit(1);
    });

    let result = run(&config);
    let json = serde_json::to_string_pretty(&result).expect("SimResult is always serializable");
    println!("{json}");
}
