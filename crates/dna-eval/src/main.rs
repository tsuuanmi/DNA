//! `dna-eval`: the evaluation harness's command (ADR-0075). It records what
//! it was built from, which its build script stamps, and runs the harness of
//! `dna-eval-lib`. The stamps change with every commit; keeping them in this
//! thin binary relinks it without recompiling the harness.

use std::process::ExitCode;

use dna_eval_lib::Identity;

fn main() -> ExitCode {
    dna_eval_lib::main(&Identity {
        name: env!("CARGO_PKG_NAME"),
        version: env!("CARGO_PKG_VERSION"),
        revision: env!("DNA_EVAL_REVISION"),
        dirty: env!("DNA_EVAL_DIRTY"),
        read_stage: env!("DNA_EVAL_READ_STAGE"),
        placement_stage: env!("DNA_EVAL_PLACEMENT_STAGE"),
        forms_stage: env!("DNA_EVAL_FORMS_STAGE"),
    })
}
