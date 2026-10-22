//! Records the source revision the harness is built from, so a run reports
//! the code that actually ran: `DNA_EVAL_REVISION` (the commit, empty without
//! git) and `DNA_EVAL_DIRTY` (`1` when tracked files differ from it, `0` when
//! they do not, empty without git). Untracked files do not count.
//!
//! It also records the code identity of each cached stage, from the sources
//! `stages.txt` lists for it: `path=object;` for every listed path, the
//! object being its git id at `HEAD`. A stage's identity is empty when one of
//! its sources differs from `HEAD`, or without git, so its cache is refused
//! (or, for the forms index, not stored) rather than keyed on code that is
//! not the commit's: `DNA_EVAL_READ_STAGE` (from a trace to its read evidence
//! document), `DNA_EVAL_PLACEMENT_STAGE` (from a read evidence document to
//! its placed-read document, `dna.placed_read/v1`), and
//! `DNA_EVAL_FORMS_STAGE` (the index of the truth's written forms,
//! `forms.tsv`). `tools/python/scripts/validate_module_layers.py` checks the
//! lists are closed under the module graph.
//!
//! The stamps change with every commit, so they live in this thin binary:
//! the harness itself, `dna-eval-lib`, takes them at run time and is not
//! recompiled.

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
    // Each stage's sources, from `stages.txt`: a `[stage]` header, then one
    // path per line; `path: items` lists the file, and `!path` lists what
    // the stage does not depend on.
    let stages = std::fs::read_to_string(crate_directory.join("stages.txt")).unwrap_or_default();
    let mut paths: Vec<(String, Vec<String>)> = Vec::new();
    for line in stages.lines() {
        let line = line.split('#').next().unwrap_or_default().trim();
        if let Some(stage) = line
            .strip_prefix('[')
            .and_then(|line| line.strip_suffix(']'))
        {
            paths.push((stage.to_owned(), Vec::new()));
        } else if !line.is_empty() && !line.starts_with('!') {
            let path = line.split(':').next().unwrap_or_default().trim();
            if let Some((_, listed)) = paths.last_mut() {
                listed.push(path.to_owned());
            }
        }
    }
    for (stage, variable) in [
        ("read", "DNA_EVAL_READ_STAGE"),
        ("placement", "DNA_EVAL_PLACEMENT_STAGE"),
        ("forms", "DNA_EVAL_FORMS_STAGE"),
    ] {
        let listed = paths
            .iter()
            .find(|(name, _)| name == stage)
            .map(|(_, listed)| listed.iter().map(String::as_str).collect::<Vec<_>>())
            .unwrap_or_default();
        let identity = if listed.is_empty() {
            String::new()
        } else {
            stage_identity(&listed)
        };
        println!("cargo:rustc-env={variable}={identity}");
    }

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
    for name in [
        "src",
        "crates",
        "Cargo.toml",
        "Cargo.lock",
        "rust-toolchain.toml",
    ] {
        watched.push(workspace.join(name));
    }
    for path in watched.iter().filter(|path| Path::exists(path)) {
        println!("cargo:rerun-if-changed={}", path.display());
    }
}
