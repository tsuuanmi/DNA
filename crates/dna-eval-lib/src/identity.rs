//! What the `dna-eval` binary was built from, as its build script recorded
//! it. The library takes it at run time, so a new commit relinks the thin
//! binary without recompiling the harness.

/// The build of the `dna-eval` binary: its name and version, the source
/// revision, and the code identity of each cached stage. Each identity is
/// `path=object;` for every source the stage depends on (the git object id
/// at the commit), empty when one of them differed from the commit or
/// without git, so its cache is refused rather than keyed on code that is
/// not the commit's.
#[derive(Clone, Copy, Debug)]
pub struct Identity {
    /// The binary's package name, recorded in `run.json` as the harness.
    pub name: &'static str,
    /// The binary's package version.
    pub version: &'static str,
    /// The commit the binary was built from; empty without git.
    pub revision: &'static str,
    /// `1` when tracked files differed from the commit, `0` when they did
    /// not, empty without git.
    pub dirty: &'static str,
    /// The read stage's code identity: from a trace to its read evidence
    /// document.
    pub read_stage: &'static str,
    /// The placement stage's code identity: from a read evidence document to
    /// its placed-read document.
    pub placement_stage: &'static str,
    /// The code identity of the index of the truth's written forms.
    pub forms_stage: &'static str,
}
