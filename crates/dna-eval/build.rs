//! Records the source revision the harness is built from, so a run reports
//! the code that actually ran: `DNA_EVAL_REVISION` (the commit, empty without
//! git) and `DNA_EVAL_DIRTY` (`1` when tracked files differ from it, `0` when
//! they do not, empty without git). Untracked files do not count.
//!
//! It also records the read stage's code identity, `DNA_EVAL_READ_STAGE`:
//! `path=object;` for every source on the path from a trace to its read
//! evidence document, the object being its git id at `HEAD`. It is empty when one
//! of them differs from `HEAD` or without git, so the evidence cache is
//! refused rather than keyed on code that is not the commit's.
//!
//! `DNA_EVAL_PLACEMENT_STAGE` records the placement stage's identity the
//! same way: the sources from a read evidence document to its placed-read
//! document (`dna.placed_read/v1`), so the placement cache misses when the
//! placement's code changes.
//!
//! `DNA_EVAL_FORMS_STAGE` records the identity of the index of the truth's
//! written forms (`forms.tsv`) the same way, so a cached index misses when
//! the code that builds it changes.

use std::path::{Path, PathBuf};
use std::process::Command;

fn main() {
    let crate_directory = std::env::var_os("CARGO_MANIFEST_DIR")
        .map(PathBuf::from)
        .unwrap_or_default();
    let git = |arguments: &[&str]| {
        Command::new("git")
            .arg("-C")
            .arg(&crate_directory)
            .args(arguments)
            .output()
            .ok()
            .filter(|output| output.status.success())
            .map(|output| String::from_utf8_lossy(&output.stdout).trim().to_owned())
    };
    let revision = git(&["rev-parse", "HEAD"]);
    let dirty = git(&["status", "--porcelain", "--untracked-files=no"]);
    println!(
        "cargo:rustc-env=DNA_EVAL_REVISION={}",
        revision.unwrap_or_default()
    );
    println!(
        "cargo:rustc-env=DNA_EVAL_DIRTY={}",
        dirty.map_or("", |status| if status.is_empty() { "0" } else { "1" })
    );
    // The read stage: the kernel, the Sanger plugin and the trace library it
    // calls, and the facade's trace loading, configuration, evidence models,
    // and evidence writer. The read stage reads no profile or reference.
    let read_stage = [
        "Cargo.lock",
        "crates/dna-kernel",
        "crates/dna-sanger",
        "crates/dna-trace",
        "src/config",
        "src/input/mod.rs",
        "src/input/sanger",
        "src/model/basecall_result.rs",
        "src/model/evidence_result.rs",
        "src/model/read_observation.rs",
        "src/model/result.rs",
        "src/model/sample_result.rs",
        "src/model/sanger_read.rs",
        "src/pipeline/evidence.rs",
        "src/pipeline/plugins.rs",
        "src/report/evidence.rs",
        "src/report/json.rs",
    ];
    let stage_identity = |paths: &[&str]| {
        let pathspecs = paths
            .iter()
            .map(|path| format!(":(top){path}"))
            .collect::<Vec<_>>();
        let mut status = vec!["status", "--porcelain", "--untracked-files=no", "--"];
        status.extend(pathspecs.iter().map(String::as_str));
        let clean = git(&status).is_some_and(|status| status.is_empty());
        paths
            .iter()
            .map(|path| {
                git(&["rev-parse", &format!("HEAD:{path}")])
                    .filter(|object| !object.is_empty())
                    .map(|object| format!("{path}={object};"))
            })
            .collect::<Option<String>>()
            .filter(|_| clean)
            .unwrap_or_default()
    };
    let identity = stage_identity(&read_stage);
    println!("cargo:rustc-env=DNA_EVAL_READ_STAGE={identity}");
    // The placement stage: the kernel and the aligner, and the facade's
    // configuration, evidence document reading, and placed-read writer.
    let placement_stage = [
        "Cargo.lock",
        "crates/dna-kernel",
        "crates/dna-align",
        "src/config",
        "src/input/evidence.rs",
        "src/input/mod.rs",
        "src/recorded.rs",
    ];
    println!(
        "cargo:rustc-env=DNA_EVAL_PLACEMENT_STAGE={}",
        stage_identity(&placement_stage)
    );
    // The index of the truth's written forms: the truth reader, the alleles
    // and anchors it compares, and the index. The reference and the
    // profile's windows it reads are keyed by their content.
    let forms_stage = [
        "Cargo.lock",
        "crates/dna-eval/src/cause/forms.rs",
        "crates/dna-eval/src/compare.rs",
        "crates/dna-eval/src/truth.rs",
    ];
    println!(
        "cargo:rustc-env=DNA_EVAL_FORMS_STAGE={}",
        stage_identity(&forms_stage)
    );

    // Record again when the commit, the index, or a source of the workspace
    // changes.
    let mut watched = Vec::new();
    let branch = git(&["rev-parse", "--symbolic-full-name", "HEAD"]);
    for name in ["HEAD", "index", "packed-refs"]
        .into_iter()
        .map(str::to_owned)
        .chain(branch.filter(|branch| branch.starts_with("refs/")))
    {
        if let Some(path) = git(&["rev-parse", "--git-path", &name]) {
            watched.push(crate_directory.join(path));
        }
    }
    let workspace = crate_directory.join("../..");
    for name in ["src", "crates", "Cargo.toml", "Cargo.lock"] {
        watched.push(workspace.join(name));
    }
    for path in watched.iter().filter(|path| Path::exists(path)) {
        println!("cargo:rerun-if-changed={}", path.display());
    }
}
