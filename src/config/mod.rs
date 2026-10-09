//! Strict, reproducible scientific configuration: the envelope that composes
//! the sections each plugin owns.

mod defaults;
mod load;
mod types;

pub(crate) use load::{load_path, resolve_path};
pub(crate) use types::Config;
