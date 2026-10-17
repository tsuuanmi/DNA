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
    /// The lab convention for a repeat run split by an anchor: the first anchor
    /// base splits the window; the left run's length changes at its 3' end
    /// (`309.1C`, `309DEL`), and the rest is written left-anchored against the
    /// right run, as substitutions plus repeat-base insertions or deletions at
    /// the window's 3' end (`315.1C`, `312T 315.1C 315.2C 315.3C`, `315DEL`).
    AnchoredRunLengths(Anchor),
    /// The lab convention for a lead run before an anchored repeat run, as in
    /// HVS-I's A-tract and C-tract (`16180..=16193`). While the anchor base
    /// sits before a right run of its reference length, the lead run's length
    /// changes at its 3' end (`16183DEL`), the left run is written
    /// left-anchored with repeat-base insertions or deletions at its 3' end
    /// (`16188.1C`), and the right run as substitutions. Otherwise the whole
    /// window is written left-anchored, as substitutions plus repeat-base
    /// insertions or deletions at the window's 3' end (`16182C 16183C 16189C`,
    /// `16189C 16190T 16193.1C 16193.2C`).
    AnchoredTract {
        /// The anchored repeat runs.
        anchor: Anchor,
        /// Length of the lead run that starts the window.
        lead: usize,
    },
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
    /// Loss or gain of whole copies of the terminal tandem-repeat motif is
    /// written at the window's 3' end (`523DEL 524DEL`, `524.1A 524.2C`) plus
    /// position-wise substitutions against the remaining reference
    /// (`513A 523DEL 524DEL`); an unchanged length is substitutions only.
    MotifShift {
        /// The repeated motif the window ends with.
        motif: String,
    },
}
