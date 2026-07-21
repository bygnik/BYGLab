//! Thin wasm-bindgen wrapper around `byglab-core`. No solver logic lives
//! here — only (de)serialization and marshalling into JS-friendly types.

use byglab_core::{run as core_run, SimConfig};
use js_sys::{Float64Array, Object, Reflect};
use wasm_bindgen::prelude::*;

/// Call once from JS before anything else, to get Rust panics forwarded to
/// the browser console instead of an opaque "unreachable" trap.
#[wasm_bindgen(start)]
pub fn init_panic_hook() {
    console_error_panic_hook::set_once();
}

/// Runs the (placeholder) solver for the given JSON-encoded `SimConfig` and
/// returns `{ theta: Float64Array, value: Float64Array }`.
#[wasm_bindgen]
pub fn run(config_json: &str) -> Result<JsValue, JsValue> {
    let config: SimConfig = serde_json::from_str(config_json)
        .map_err(|e| JsValue::from_str(&format!("invalid config: {e}")))?;

    let result = core_run(&config);

    let out = Object::new();
    Reflect::set(&out, &"theta".into(), &Float64Array::from(result.theta.as_slice()))?;
    Reflect::set(&out, &"value".into(), &Float64Array::from(result.value.as_slice()))?;
    Ok(out.into())
}
