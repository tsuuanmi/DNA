//! Validated nomenclature windows.
//!
//! Each rule carries the structure it needs, so a rule can never be paired with
//! an incompatible window structure.

/// One reference window with sequence-equivalent representation rules.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NomenclatureWindow {
    /// Window name, unique within the profile.
    pub name: String,
    /// 0-based index of the first window base.
    pub start: usize,
    /// Reference bases the window is validated against.
    pub sequence: String,
    /// Rules tried in order; the first that reproduces the haplotype wins.
    pub rules: Vec<WindowRule>,
}

impl NomenclatureWindow {
    /// 0-based exclusive end of the window.
    #[must_use]
    pub fn end(&self) -> usize {
        self.start + self.sequence.len()
    }
}

/// Two runs of one repeated base split by a single anchor base.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Anchor {
    /// Window-local index of the anchor base.
    pub index: usize,
    /// The anchor base, for example `T`.
    pub base: u8,
    /// The repeated base, for example `C`.
    pub repeat: u8,
}

/// Sequence-preserving representation rules.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WindowRule {
    /// Express anchor movement as run-length change at the run ends
    /// (`309.1C`, `315.1C`, `309DEL`).
    AnchoredRunLengths(Anchor),
    /// A duplicated anchor is the base after it plus a terminal repeat base
    /// (`311T 315.1C`).
    AnchorDuplication(Anchor),
    /// A deleted anchor becomes a repeat base with the last run base deleted
    /// (`310C 315DEL`, `16189C 16193DEL`).
    AnchorDeletion(Anchor),
    /// A one-base gain ending in the repeat base keeps the terminal insertion
    /// (EMPOP `315.1C`) and describes the rest as substitutions.
    TerminalRepeatInsertion(Anchor),
    /// A validated representation, as window-local `(index, base)` substitutions
    /// (`16183C 16184A 16189C`).
    CanonicalHaplotype(Vec<(usize, u8)>),
    /// Loss of one copy of the terminal tandem-repeat motif is that copy's
    /// deletion plus position-wise substitutions (`513A 523DEL 524DEL`).
    MotifShift {
        /// The repeated motif the window ends with.
        motif: String,
    },
}
