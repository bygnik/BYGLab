//! Thin wasm-bindgen wrapper around `byglab-core`. No solver logic lives
//! here — only (de)serialization and marshalling into JS-friendly types.

use byglab_core::{run_pipe_case, PipeCaseConfig};
use wasm_bindgen::prelude::*;

/// Call once from JS before anything else, to get Rust panics forwarded to
/// the browser console instead of an opaque "unreachable" trap.
#[wasm_bindgen(start)]
pub fn init_panic_hook() {
    console_error_panic_hook::set_once();
}

/// Runs a pipe-network simulation for the given JSON-encoded
/// [`PipeCaseConfig`] and returns the JSON-encoded [`byglab_core::PipeCaseResult`].
#[wasm_bindgen]
pub fn run(config_json: &str) -> Result<JsValue, JsValue> {
    let config: PipeCaseConfig =
        serde_json::from_str(config_json).map_err(|e| JsValue::from_str(&format!("invalid config: {e}")))?;

    let result = run_pipe_case(&config);

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsValue::from_str(&format!("could not serialize result: {e}")))
}
