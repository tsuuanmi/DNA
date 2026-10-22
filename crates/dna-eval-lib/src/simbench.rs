//! `dna-eval sim-bench`: the simulated tract-decode bench (PROP-0015
//! increment 1).
//!
//! Each case renders a read through one C-tract as a mixture of molecules
//! whose run lengths follow a set distribution (the truth), with read
//! parameters drawn from `dna-trace`'s aggregate table (feature
//! `simulate`). It writes the read as ABIF, runs DNA's read evidence stage
//! on it in process, and scores the tract's decode against the truth: the
//! chosen structure, the dominant length of the marginal that carries the
//! varying run, the margin state, the ratio of the shorter neighbour to the
//! dominant, and the direction of exact ties. The case set is frozen, so the
//! numbers are the decoder's own: a baseline later decoder changes are
//! measured against. The comparison with simulated truth lives here, in the
//! harness, never in a DNA crate.
//!
//! It ports the research benchmark (`tools/python/research/trace_simulator`,
//! `bench.py` and `report.py`).

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::Path;
use std::sync::Mutex;

use dna_trace::simulate::{Rng, Table, render};
use serde_json::Value;

use crate::runner::Threads;

/// The decoder's measures on the gating set at the commit that introduced
/// the bench (PROP-0015 increment 1): the baseline later decoder changes
/// are decided against.
const BASELINE_GATE: &str = include_str!("../baselines/sim-bench.tsv");
/// The same on the smoke set.
const BASELINE_SMOKE: &str = include_str!("../baselines/sim-bench-smoke.tsv");

/// The committed baseline for a run of `reps` repetitions per cell: the
/// gating or the smoke set's; `None` for any other size.
pub(crate) const fn baseline(reps: usize) -> Option<&'static str> {
    match reps {
        GATE_REPS => Some(BASELINE_GATE),
        SMOKE_REPS => Some(BASELINE_SMOKE),
        _ => None,
    }
}

/// The frozen case set's version: raise it whenever the cases change.
pub(crate) const CASES_VERSION: u32 = 1;
/// Repetitions per cell of the gating set (3,960 reads): later decoder
/// increments decide on it, a cell's 240-300 reads resolving a rate near
/// 1 % (PROP-0015).
pub(crate) const GATE_REPS: usize = 60;
/// Repetitions per cell of the smoke set (396 reads), a quick check.
pub(crate) const SMOKE_REPS: usize = 6;
/// The pair family's neighbour ratios.
const PAIR_VALUES: [f64; 5] = [0.5, 0.8, 0.9, 1.0, 1.25];
/// Noise levels: the read's own, twice and four times it.
const NOISE: [f64; 3] = [1.0, 2.0, 4.0];
/// Pure reads (no neighbour) per tract and noise level in the smoke set;
/// half the repetitions in any set.
const PURE_REPS: usize = 3;

/// Which marginal carries the varying run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Key {
    First,
    Total,
}

impl Key {
    const fn name(self) -> &'static str {
        match self {
            Self::First => "first",
            Self::Total => "total",
        }
    }
}

/// A C-tract read by one primer: where, how its molecules are written, its
/// truth.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Tract {
    pub(crate) name: &'static str,
    primer: &'static str,
    forward: bool,
    /// The rCRS position of call 0.
    start: usize,
    /// The tract's rCRS span.
    positions: (usize, usize),
    /// Calls the read holds after the tract.
    after: usize,
    /// Whether the varying run is the read's first run (the forward
    /// reads), or the one after the mark (HVS-II reverse).
    varying_first: bool,
    /// The other run's length (0: no mark).
    fixed: usize,
    /// The scored marginal.
    pub(crate) key: Key,
    /// The structure the decode should choose.
    mark: Option<char>,
    /// The dominant length of the varying run.
    dominant: usize,
}

/// The bench's tracts: HVS-II (`C^a T C^6`, 315.1C) forward and reverse,
/// and HVS-I with T16189C (one run) forward and reverse.
pub(crate) const TRACTS: [Tract; 4] = [
    Tract {
        name: "HV2-F",
        primer: "HV2F",
        forward: true,
        start: 36,
        positions: (303, 315),
        after: 200,
        varying_first: true,
        fixed: 6,
        key: Key::Total,
        mark: Some('T'),
        dominant: 9,
    },
    Tract {
        name: "HV2-R",
        primer: "HV3R",
        forward: false,
        start: 604,
        positions: (303, 316),
        after: 200,
        varying_first: false,
        fixed: 6,
        key: Key::Total,
        mark: Some('A'),
        dominant: 9,
    },
    Tract {
        name: "HV1-F",
        primer: "HV1F",
        forward: true,
        start: 16_002,
        positions: (16_184, 16_193),
        after: 200,
        varying_first: true,
        fixed: 0,
        key: Key::First,
        mark: None,
        dominant: 11,
    },
    Tract {
        name: "HV1-R",
        primer: "HV1R",
        forward: false,
        start: 16_401,
        positions: (16_184, 16_193),
        after: 200,
        varying_first: true,
        fixed: 0,
        key: Key::First,
        mark: None,
        dominant: 11,
    },
];

impl Tract {
    /// rCRS edits (1-based position → replacing text on the reference
    /// strand) giving the varying run `length` C.
    fn edits(&self, length: usize) -> BTreeMap<usize, String> {
        let c = |count: usize| "C".repeat(count);
        if self.name.starts_with("HV2") {
            // Run 1 is 303-309 (C7), T310, run 2 is 311-315 (C5): run 1
            // varies and run 2 holds `fixed` C on both strands, so a
            // reverse read meets the fixed run first.
            BTreeMap::from([
                (309, c(length.saturating_sub(6))),
                (315, c(self.fixed.saturating_sub(4))),
            ])
        } else {
            // T16189C joins 16184-16193 into one run of ten.
            BTreeMap::from([(16_189, c(1)), (16_193, c(length.saturating_sub(9)))])
        }
    }

    /// The (first, total) lengths a molecule with varying run `length`
    /// shows.
    const fn marginals(&self, length: usize) -> (usize, usize) {
        if self.fixed == 0 {
            (length, length)
        } else if self.varying_first {
            (length, length + 1 + self.fixed)
        } else {
            (self.fixed, self.fixed + 1 + length)
        }
    }
}

/// One case of the bench.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Case {
    pub(crate) tract: usize,
    /// The shorter neighbour's ratio to the dominant; 0 for a pure read.
    pub(crate) value: f64,
    pub(crate) noise: f64,
    pub(crate) rep: usize,
}

impl Case {
    pub(crate) fn id(&self) -> String {
        format!(
            "{}_pair_{}_{}_n{}_r{}",
            TRACTS[self.tract].name, TRACTS[self.tract].dominant, self.value, self.noise, self.rep
        )
    }

    /// Molecule fractions over the varying run's lengths.
    fn fractions(&self) -> BTreeMap<usize, f64> {
        let dominant = TRACTS[self.tract].dominant;
        let total = 1.0 + self.value;
        let mut fractions = BTreeMap::from([(dominant, 1.0 / total)]);
        if self.value > 0.0 {
            fractions.insert(dominant - 1, self.value / total);
        }
        fractions
    }

    /// The case's seed: a hash of its id, so a case renders alike in every
    /// run whatever the case order.
    fn seed(&self) -> u64 {
        self.id().bytes().fold(0xCBF2_9CE4_8422_2325, |hash, byte| {
            (hash ^ u64::from(byte)).wrapping_mul(0x100_0000_01B3)
        })
    }
}

/// The case set: per tract, the pair family at each noise level and ratio,
/// `reps` each, then pure reads.
pub(crate) fn cases(reps: usize) -> Vec<Case> {
    let mut cases = Vec::new();
    for tract in 0..TRACTS.len() {
        for noise in NOISE {
            for value in PAIR_VALUES {
                for rep in 0..reps {
                    cases.push(Case {
                        tract,
                        value,
                        noise,
                        rep,
                    });
                }
            }
            for rep in 0..PURE_REPS * reps / SMOKE_REPS {
                cases.push(Case {
                    tract,
                    value: 0.0,
                    noise,
                    rep,
                });
            }
        }
    }
    cases
}

/// The read-order template from 1-based `start` on a strand, `length`
/// bases, with `edits` applied on the reference strand, and the rCRS
/// position each base comes from.
fn template(
    reference: &[u8],
    forward: bool,
    start: usize,
    length: usize,
    edits: &BTreeMap<usize, String>,
) -> (Vec<u8>, Vec<usize>) {
    let mut bases = Vec::with_capacity(length);
    let mut positions = Vec::with_capacity(length);
    let mut position = start;
    while bases.len() < length && (1..=reference.len()).contains(&position) {
        let piece = edits.get(&position).map_or_else(
            || vec![reference[position - 1]],
            |text| text.as_bytes().to_vec(),
        );
        let piece = if forward {
            piece
        } else {
            piece
                .iter()
                .rev()
                .map(|base| match base {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'T' => b'A',
                    other => *other,
                })
                .collect()
        };
        for base in piece {
            bases.push(base);
            positions.push(position);
        }
        if forward {
            position += 1;
        } else {
            position -= 1;
        }
    }
    bases.truncate(length);
    positions.truncate(length);
    (bases, positions)
}

/// One ABIF directory entry and its data.
struct Record {
    tag: [u8; 4],
    number: u32,
    kind: u16,
    size: u16,
    data: Vec<u8>,
}

/// A size or offset as an ABIF 32-bit word; a simulated read is far below
/// its range.
fn word(value: usize) -> u32 {
    u32::try_from(value).unwrap_or(u32::MAX)
}

/// An analysed-channel sample as ABIF's 16-bit signed value, rounded and
/// clipped to its range as the instrument's analysed data are.
#[allow(clippy::cast_possible_truncation)]
fn short(value: f64) -> [u8; 2] {
    (value.round_ties_even().clamp(-32_768.0, 32_767.0) as i16).to_be_bytes()
}

/// A minimal ABIF: the analysed channels as DATA.9-12 in GATC order, FWO_,
/// and the calls as PLOC.2, PBAS.2, and PCON.2.
pub(crate) fn abif(channels: &[Vec<f64>; 4], ploc: &[usize], pbas: &[u8], pcon: &[u8]) -> Vec<u8> {
    let order = *b"GATC";
    let channel = |base: u8| match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        _ => 3,
    };
    let mut records = order
        .iter()
        .zip(9..)
        .map(|(&base, number)| Record {
            tag: *b"DATA",
            number,
            kind: 4,
            size: 2,
            data: channels[channel(base)]
                .iter()
                .flat_map(|&value| short(value))
                .collect(),
        })
        .collect::<Vec<_>>();
    records.push(Record {
        tag: *b"FWO_",
        number: 1,
        kind: 2,
        size: 1,
        data: order.to_vec(),
    });
    records.push(Record {
        tag: *b"PLOC",
        number: 2,
        kind: 4,
        size: 2,
        data: ploc
            .iter()
            .flat_map(|&sample| short(sample as f64))
            .collect(),
    });
    records.push(Record {
        tag: *b"PBAS",
        number: 2,
        kind: 2,
        size: 1,
        data: pbas.to_vec(),
    });
    records.push(Record {
        tag: *b"PCON",
        number: 2,
        kind: 2,
        size: 1,
        data: pcon.to_vec(),
    });
    let directory = 128_usize;
    let mut offset = directory + 28 * records.len();
    let offsets = records
        .iter()
        .map(|record| {
            let at = offset;
            if record.data.len() > 4 {
                offset += record.data.len();
            }
            at
        })
        .collect::<Vec<_>>();
    let mut out = vec![0_u8; offset];
    out[0..4].copy_from_slice(b"ABIF");
    out[4..6].copy_from_slice(&101_u16.to_be_bytes());
    let count = word(records.len());
    let header = [
        b"tdir".to_vec(),
        1_u32.to_be_bytes().to_vec(),
        1023_u16.to_be_bytes().to_vec(),
        28_u16.to_be_bytes().to_vec(),
        count.to_be_bytes().to_vec(),
        (28 * count).to_be_bytes().to_vec(),
        word(directory).to_be_bytes().to_vec(),
    ]
    .concat();
    out[6..6 + header.len()].copy_from_slice(&header);
    for (index, (record, &at)) in records.iter().zip(&offsets).enumerate() {
        let mut entry = Vec::with_capacity(28);
        entry.extend_from_slice(&record.tag);
        entry.extend_from_slice(&record.number.to_be_bytes());
        entry.extend_from_slice(&record.kind.to_be_bytes());
        entry.extend_from_slice(&record.size.to_be_bytes());
        entry.extend_from_slice(&word(record.data.len() / usize::from(record.size)).to_be_bytes());
        entry.extend_from_slice(&word(record.data.len()).to_be_bytes());
        if record.data.len() > 4 {
            entry.extend_from_slice(&word(at).to_be_bytes());
            out[at..at + record.data.len()].copy_from_slice(&record.data);
        } else {
            let mut inline = record.data.clone();
            inline.resize(4, 0);
            entry.extend_from_slice(&inline);
        }
        entry.extend_from_slice(&[0; 4]);
        let slot = directory + index * 28;
        out[slot..slot + 28].copy_from_slice(&entry);
    }
    out
}

/// One decoded marginal as the read evidence records it.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Decoded {
    shares: BTreeMap<usize, f64>,
    dominant: Option<usize>,
    margin: String,
}

/// What DNA made of one case.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Outcome {
    pub(crate) case: Case,
    /// The truth on the scored marginal: length → share.
    truth: BTreeMap<usize, f64>,
    /// The decode at the tract, when there is one: its mark and the scored
    /// marginal.
    decode: Option<(Option<char>, Decoded)>,
    /// Why DNA produced no read evidence, if it did not.
    error: Option<String>,
}

/// The decode of the tract opening at `call` in a read evidence document.
fn decode_at(document: &Value, call: usize, key: Key) -> Option<(Option<char>, Decoded)> {
    let tracts = document.get("read")?.get("tracts")?.as_array()?;
    let tract = tracts.iter().find(|tract| {
        tract
            .pointer("/calls/start")
            .and_then(Value::as_u64)
            .and_then(|start| usize::try_from(start).ok())
            .is_some_and(|start| start.abs_diff(call) <= 3)
    })?;
    let mark = tract
        .get("mark")
        .and_then(Value::as_str)
        .and_then(|mark| mark.chars().next());
    let marginal = tract.get(key.name())?;
    let lengths = marginal.get("lengths")?.as_array()?;
    let shares = marginal.get("shares")?.as_array()?;
    Some((
        mark,
        Decoded {
            shares: lengths
                .iter()
                .zip(shares)
                .filter_map(|(length, share)| {
                    Some((usize::try_from(length.as_u64()?).ok()?, share.as_f64()?))
                })
                .collect(),
            dominant: marginal
                .get("dominant")
                .and_then(Value::as_u64)
                .and_then(|dominant| usize::try_from(dominant).ok()),
            margin: marginal
                .get("margin")
                .and_then(Value::as_str)
                .unwrap_or_default()
                .to_owned(),
        },
    ))
}

/// Renders and runs one case.
fn run_case(
    case: &Case,
    table: &Table,
    reference: &[u8],
    config: &Path,
    (scratch, keep): (&Path, bool),
) -> Outcome {
    let tract = &TRACTS[case.tract];
    let fractions = case.fractions();
    let mut truth = BTreeMap::new();
    for (&length, &share) in &fractions {
        let (first, total) = tract.marginals(length);
        *truth
            .entry(match tract.key {
                Key::First => first,
                Key::Total => total,
            })
            .or_insert(0.0) += share;
    }
    let failed = |error: String| Outcome {
        case: case.clone(),
        truth: truth.clone(),
        decode: None,
        error: Some(error),
    };
    let longest = fractions.keys().max().copied().unwrap_or(tract.dominant);
    let (_, positions) = template(
        reference,
        tract.forward,
        tract.start,
        2_000,
        &tract.edits(longest),
    );
    let Some(tract_call) = positions.iter().position(|&position| {
        if tract.forward {
            position >= tract.positions.0
        } else {
            position <= tract.positions.1
        }
    }) else {
        return failed("the template does not reach the tract".into());
    };
    let calls = tract_call + 20 + tract.after;
    let templates = fractions
        .iter()
        .map(|(&length, &share)| {
            let (bases, _) = template(
                reference,
                tract.forward,
                tract.start,
                calls + 5,
                &tract.edits(length),
            );
            (bases, share)
        })
        .collect::<Vec<_>>();
    let mut rng = Rng::new(case.seed());
    let Some(params) = table.params(tract.primer, &mut rng) else {
        return failed(format!("the table has no primer {}", tract.primer));
    };
    let params = params.scaled_noise(case.noise);
    let borrowed = templates
        .iter()
        .map(|(bases, share)| (bases.as_slice(), *share))
        .collect::<Vec<_>>();
    let rendered = render(&borrowed, &params, Some(calls), &[]);
    let path = scratch.join(format!("{}.ab1", case.id()));
    let bytes = abif(
        &rendered.channels,
        &rendered.ploc,
        &rendered.pbas,
        &rendered.pcon,
    );
    if let Err(error) = std::fs::write(&path, bytes) {
        return failed(error.to_string());
    }
    let document = dna::recorded::read_evidence_document(&path, config);
    if !keep {
        let _ = std::fs::remove_file(&path);
    }
    let document = match document
        .map_err(|error| error.to_string())
        .and_then(|bytes| {
            serde_json::from_slice::<Value>(&bytes).map_err(|error| error.to_string())
        }) {
        Ok(document) if document.get("read").is_some() => document,
        Ok(_) => return failed("the read was set aside".into()),
        Err(error) => return failed(error),
    };
    Outcome {
        case: case.clone(),
        decode: decode_at(&document, tract_call, tract.key),
        truth,
        error: None,
    }
}

/// Runs every case on `threads` workers, in case order.
pub(crate) fn run(
    cases: &[Case],
    table: &Table,
    reference: &[u8],
    config: &Path,
    scratch: (&Path, bool),
    threads: Threads,
) -> Vec<Outcome> {
    let next = Mutex::new(0_usize);
    let results = Mutex::new(vec![None; cases.len()]);
    std::thread::scope(|scope| {
        for _ in 0..threads.get().min(cases.len().max(1)) {
            scope.spawn(|| {
                loop {
                    let index = {
                        let mut next = next
                            .lock()
                            .unwrap_or_else(std::sync::PoisonError::into_inner);
                        let index = *next;
                        *next += 1;
                        index
                    };
                    let Some(case) = cases.get(index) else {
                        break;
                    };
                    let outcome = run_case(case, table, reference, config, scratch);
                    results
                        .lock()
                        .unwrap_or_else(std::sync::PoisonError::into_inner)[index] = Some(outcome);
                }
            });
        }
    });
    results
        .into_inner()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
        .into_iter()
        .flatten()
        .collect()
}

/// One decode's scores against the truth.
struct Score {
    right: bool,
    shorter: bool,
    clear: bool,
    rho_hat: Option<f64>,
}

fn score(outcome: &Outcome) -> Option<(Option<char>, Score)> {
    let (mark, decoded) = outcome.decode.as_ref()?;
    let top = outcome.truth.values().copied().fold(0.0_f64, f64::max);
    let dominants = outcome
        .truth
        .iter()
        .filter(|(_, share)| (**share - top).abs() < 1e-9)
        .map(|(length, _)| *length)
        .collect::<Vec<_>>();
    let dominant = decoded.dominant?;
    // The truth's dominant length, the shorter one on a tie.
    let big = dominants.first().copied()?;
    let rho_hat = (outcome.truth.contains_key(&(big - 1)))
        .then(|| decoded.shares.get(&big).copied())
        .flatten()
        .filter(|share| *share > 0.0)
        .map(|share| decoded.shares.get(&(big - 1)).copied().unwrap_or(0.0) / share);
    Some((
        *mark,
        Score {
            right: dominants.contains(&dominant),
            shorter: Some(&dominant) == dominants.first(),
            clear: decoded.margin == "clear",
            rho_hat,
        },
    ))
}

/// The bench's measures per tract and noise level.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Row {
    pub(crate) tract: String,
    pub(crate) marginal: String,
    pub(crate) noise: f64,
    pub(crate) values: BTreeMap<&'static str, f64>,
}

/// The measures, in column order.
pub(crate) const MEASURES: [&str; 9] = [
    "reads",
    "decoded",
    "structure_right",
    "dominant_right",
    "right_given_structure",
    "clear_wrong",
    "rho_ratio",
    "tie_shorter",
    "pure_decoded",
];

fn percent(values: &[bool]) -> f64 {
    if values.is_empty() {
        f64::NAN
    } else {
        100.0 * values.iter().filter(|value| **value).count() as f64 / values.len() as f64
    }
}

fn median(mut values: Vec<f64>) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    values.sort_by(f64::total_cmp);
    let middle = values.len() / 2;
    if values.len() % 2 == 1 {
        values[middle]
    } else {
        f64::midpoint(values[middle - 1], values[middle])
    }
}

/// The measures over the outcomes, per tract and noise level (pair family,
/// a neighbour present; pure reads only for `pure_decoded`).
pub(crate) fn measures(outcomes: &[Outcome]) -> Vec<Row> {
    let mut rows = Vec::new();
    for (index, tract) in TRACTS.iter().enumerate() {
        for noise in NOISE {
            let group = outcomes
                .iter()
                .filter(|outcome| outcome.case.tract == index && outcome.case.noise == noise)
                .collect::<Vec<_>>();
            let (pure, pair): (Vec<&Outcome>, Vec<&Outcome>) = group
                .into_iter()
                .partition(|outcome| outcome.case.value == 0.0);
            if pair.is_empty() {
                continue;
            }
            let scored = pair
                .iter()
                .filter_map(|outcome| score(outcome).map(|scored| (outcome.case.value, scored)))
                .collect::<Vec<_>>();
            let mut structure = Vec::new();
            let mut right = Vec::new();
            let mut right_given = Vec::new();
            let mut clear_wrong = Vec::new();
            let mut ratio = Vec::new();
            let mut tie = Vec::new();
            for (value, (mark, score)) in &scored {
                let ok = *mark == tract.mark;
                structure.push(ok);
                if *value == 1.0 {
                    tie.push(score.shorter);
                    continue;
                }
                right.push(score.right);
                if ok {
                    right_given.push(score.right);
                }
                clear_wrong.push(score.clear && !score.right);
                if (0.5..=0.9).contains(value)
                    && let Some(rho_hat) = score.rho_hat
                {
                    ratio.push(rho_hat / value);
                }
            }
            let values = BTreeMap::from([
                ("reads", pair.len() as f64),
                (
                    "decoded",
                    percent(
                        &pair
                            .iter()
                            .map(|outcome| outcome.decode.is_some())
                            .collect::<Vec<_>>(),
                    ),
                ),
                ("structure_right", percent(&structure)),
                ("dominant_right", percent(&right)),
                ("right_given_structure", percent(&right_given)),
                ("clear_wrong", percent(&clear_wrong)),
                ("rho_ratio", median(ratio)),
                ("tie_shorter", percent(&tie)),
                (
                    "pure_decoded",
                    percent(
                        &pure
                            .iter()
                            .map(|outcome| outcome.decode.is_some())
                            .collect::<Vec<_>>(),
                    ),
                ),
            ]);
            rows.push(Row {
                tract: tract.name.to_owned(),
                marginal: tract.key.name().to_owned(),
                noise,
                values,
            });
        }
    }
    rows
}

/// The header of `bench.tsv` and of a baseline.
pub(crate) fn header() -> String {
    format!("tract\tmarginal\tnoise\t{}\n", MEASURES.join("\t"))
}

fn cell(measure: &str, value: f64) -> String {
    if value.is_nan() {
        "-".into()
    } else if measure == "rho_ratio" {
        format!("{value:.3}")
    } else if measure == "reads" {
        format!("{value:.0}")
    } else {
        format!("{value:.1}")
    }
}

/// `bench.tsv`.
pub(crate) fn table(rows: &[Row]) -> String {
    let mut text = header();
    for row in rows {
        let _ = write!(text, "{}\t{}\t{}", row.tract, row.marginal, row.noise);
        for measure in MEASURES {
            let _ = write!(text, "\t{}", cell(measure, row.values[measure]));
        }
        text.push('\n');
    }
    text
}

/// A `bench.tsv` read back: (tract, noise) → measure → value.
pub(crate) type Baseline = BTreeMap<(String, String), BTreeMap<String, f64>>;

/// Reads a `bench.tsv` back: a cell is a finite number, or `-` where the
/// bench had nothing to measure; a cell, a row, or a repeated row of
/// another form is an error.
pub(crate) fn parse(text: &str) -> Result<Baseline, String> {
    let mut lines = text.lines();
    if lines.next().map(|line| format!("{line}\n")) != Some(header()) {
        return Err("the header is not the bench's".into());
    }
    let mut out = BTreeMap::new();
    for (number, line) in lines.enumerate() {
        let number = number + 2;
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() != 3 + MEASURES.len() {
            return Err(format!("line {number} has {} fields", fields.len()));
        }
        let values = MEASURES
            .iter()
            .zip(&fields[3..])
            .map(|(measure, value)| {
                let parsed = match *value {
                    "-" => f64::NAN,
                    value => value
                        .parse::<f64>()
                        .ok()
                        .filter(|value| value.is_finite())
                        .ok_or_else(|| {
                            format!("line {number}: {measure} {value:?} is not a number")
                        })?,
                };
                Ok(((*measure).to_owned(), parsed))
            })
            .collect::<Result<BTreeMap<_, _>, String>>()?;
        let key = (fields[0].to_owned(), fields[2].to_owned());
        if out.insert(key, values).is_some() {
            return Err(format!(
                "line {number} repeats tract {} at noise {}",
                fields[0], fields[2]
            ));
        }
    }
    Ok(out)
}

/// The measures that gate, each with whether a higher value is better.
const GATED: [(&str, bool); 5] = [
    ("decoded", true),
    ("structure_right", true),
    ("dominant_right", true),
    ("right_given_structure", true),
    ("clear_wrong", false),
];

/// The bench's verdict against `baseline`: every gated measure of a cell
/// that worsens by more than `tolerance` points, a cell the baseline lacks
/// or this run lacks, and reads without read evidence. Empty when it passes.
pub(crate) fn verdict(
    rows: &[Row],
    outcomes: &[Outcome],
    baseline: &Baseline,
    tolerance: f64,
) -> Vec<String> {
    let mut failures = Vec::new();
    let failed = outcomes
        .iter()
        .filter(|outcome| outcome.error.is_some())
        .count();
    if failed > 0 {
        failures.push(format!("{failed} reads have no read evidence"));
    }
    let mut seen = std::collections::BTreeSet::new();
    for row in rows {
        let key = (row.tract.clone(), row.noise.to_string());
        let Some(base) = baseline.get(&key) else {
            failures.push(format!(
                "{} at noise {}: not in the baseline",
                row.tract, row.noise
            ));
            continue;
        };
        seen.insert(key);
        for (measure, higher) in GATED {
            // Compared as written, as the report shows them.
            let now = cell(measure, row.values[measure]).parse::<f64>().ok();
            let before = base.get(measure).copied().filter(|value| !value.is_nan());
            let worse = match (before, now) {
                (Some(before), Some(now)) => {
                    let change = if higher { before - now } else { now - before };
                    change > tolerance + 1e-9
                }
                (Some(_), None) => true,
                (None, _) => false,
            };
            if worse {
                failures.push(format!(
                    "{} at noise {}: {measure} {} against {}",
                    row.tract,
                    row.noise,
                    now.map_or_else(|| "-".to_owned(), |now| now.to_string()),
                    before.map_or_else(|| "-".to_owned(), |before| before.to_string()),
                ));
            }
        }
    }
    for (tract, noise) in baseline.keys().filter(|key| !seen.contains(*key)) {
        failures.push(format!("{tract} at noise {noise}: not measured"));
    }
    failures
}

/// `bench.md`: the measures, and in brackets their change from the
/// baseline.
pub(crate) fn report(
    revision: &str,
    reps: usize,
    rows: &[Row],
    outcomes: &[Outcome],
    baseline: Option<&Baseline>,
) -> String {
    let failed = outcomes
        .iter()
        .filter(|outcome| outcome.error.is_some())
        .count();
    let mut text = format!(
        "# Simulated tract-decode bench\n\nRevision {revision}. Case set version {CASES_VERSION}, {} reps per cell{}: {} reads, {failed} without read evidence.\n\n",
        reps,
        match reps {
            GATE_REPS => " (the gating set)",
            SMOKE_REPS => " (the smoke set)",
            _ => " (a validation set, no committed baseline)",
        },
        outcomes.len(),
    );
    text.push_str("Pair family (the shorter neighbour at rho x the dominant), per tract and noise (x the read's own). `structure_right`: the chosen mark is the tract's. `right_given_structure`: the dominant length is right among decodes with the right structure. `clear_wrong`: a clear margin on a wrong dominant (rho 1 excluded). `rho_ratio`: median rho-hat / rho at rho 0.5-0.9. `tie_shorter`: decodes at rho 1 whose dominant is the shorter length. `pure_decoded`: pure reads (no neighbour) given a decode. Percentages.\n\n");
    text.push_str("| tract | marginal | noise |");
    for measure in MEASURES {
        let _ = write!(text, " {measure} |");
    }
    text.push_str("\n|---|---|---:|");
    text.push_str(&"---:|".repeat(MEASURES.len()));
    text.push('\n');
    for row in rows {
        let _ = write!(text, "| {} | {} | {} |", row.tract, row.marginal, row.noise);
        for measure in MEASURES {
            let value = row.values[measure];
            // Compared as written, so the baseline's rounding is no change.
            let shown = cell(measure, value);
            let delta = baseline
                .and_then(|baseline| baseline.get(&(row.tract.clone(), row.noise.to_string())))
                .and_then(|base| base.get(measure))
                .zip(shown.parse::<f64>().ok())
                .filter(|(base, _)| !base.is_nan())
                .map(|(base, now)| now - base)
                .filter(|delta| delta.abs() > 1e-9);
            let _ = write!(text, " {shown}");
            if let Some(delta) = delta {
                let delta = if measure == "rho_ratio" {
                    format!("{delta:+.3}")
                } else {
                    format!("{delta:+.1}")
                };
                let _ = write!(text, " ({delta})");
            }
            text.push_str(" |");
        }
        text.push('\n');
    }
    text
}

/// `cases.tsv`: one row per case.
pub(crate) fn cases_table(outcomes: &[Outcome]) -> String {
    let mut text =
        String::from("case\ttract\tvalue\tnoise\trep\tdecoded\tmark\tdominant\tmargin\terror\n");
    for outcome in outcomes {
        let case = &outcome.case;
        let (decoded, mark, dominant, margin) = match &outcome.decode {
            Some((mark, decoded)) => (
                "1",
                mark.map_or_else(|| "-".to_owned(), String::from),
                decoded
                    .dominant
                    .map_or_else(|| "-".to_owned(), |length| length.to_string()),
                decoded.margin.clone(),
            ),
            None => ("0", "-".into(), "-".into(), "-".into()),
        };
        let _ = writeln!(
            text,
            "{}\t{}\t{}\t{}\t{}\t{decoded}\t{mark}\t{dominant}\t{margin}\t{}",
            case.id(),
            TRACTS[case.tract].name,
            case.value,
            case.noise,
            case.rep,
            outcome
                .error
                .as_deref()
                .unwrap_or("-")
                .replace(['\t', '\n'], " "),
        );
    }
    text
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_case_sets_have_their_sizes() {
        assert_eq!(cases(SMOKE_REPS).len(), 396);
        assert_eq!(cases(GATE_REPS).len(), 3_960);
        let cases = cases(SMOKE_REPS);
        let ids = cases
            .iter()
            .map(Case::id)
            .collect::<std::collections::BTreeSet<_>>();
        assert_eq!(ids.len(), cases.len());
    }

    #[test]
    fn the_committed_baseline_covers_every_tract_and_noise_level() {
        for reps in [GATE_REPS, SMOKE_REPS] {
            let text = baseline(reps).unwrap_or_else(|| panic!("no baseline for {reps}"));
            let baseline = parse(text).unwrap_or_else(|error| panic!("{error}"));
            assert_eq!(baseline.len(), TRACTS.len() * NOISE.len());
            for tract in TRACTS {
                for noise in NOISE {
                    let row = &baseline[&(tract.name.to_owned(), noise.to_string())];
                    assert_eq!(row["reads"], (PAIR_VALUES.len() * reps) as f64);
                }
            }
        }
        assert!(baseline(7).is_none());
    }

    #[test]
    fn templates_read_the_reference_with_edits() {
        let reference = b"ACGTACGTAC";
        let (forward, positions) = template(
            reference,
            true,
            2,
            6,
            &BTreeMap::from([(4, "TT".to_owned())]),
        );
        assert_eq!(forward, b"CGTTAC");
        assert_eq!(positions, vec![2, 3, 4, 4, 5, 6]);
        let (reverse, _) = template(reference, false, 4, 4, &BTreeMap::new());
        assert_eq!(reverse, b"ACGT");
    }

    #[test]
    fn hv2_molecules_show_their_marginals() {
        let forward = &TRACTS[0];
        assert_eq!(forward.marginals(9), (9, 16));
        assert_eq!(forward.edits(9)[&309], "CCC");
        assert_eq!(forward.edits(9)[&315], "CC");
        let reverse = &TRACTS[1];
        assert_eq!(reverse.marginals(9), (6, 16));
    }

    /// The committed baseline reads strictly, and the verdict fails a gated
    /// measure that worsens, a cell gone, and nothing else.
    #[test]
    fn the_verdict_gates_a_worse_measure_against_the_baseline() {
        let text = baseline(SMOKE_REPS).unwrap_or_else(|| panic!("no smoke baseline"));
        let base = parse(text).unwrap_or_else(|error| panic!("{error}"));
        let header = header();
        let row = text.lines().nth(1).unwrap_or_default();
        for (body, message) in [
            (format!("{row}\n{row}\n"), "line 3 repeats tract"),
            (
                format!("{}\n", row.replacen("100.0", "x", 1)),
                "is not a number",
            ),
            (
                format!("{}\n", row.replacen("100.0", "NaN", 1)),
                "is not a number",
            ),
        ] {
            let error = parse(&format!("{header}{body}")).err().unwrap_or_default();
            assert!(error.contains(message), "{error}");
        }
        // Rows equal to the baseline pass; a worse dominant or a higher
        // clear_wrong fails; a better one passes.
        let rows = base
            .iter()
            .map(|((tract, noise), values)| Row {
                tract: tract.clone(),
                marginal: "total".to_owned(),
                noise: noise.parse().unwrap_or_default(),
                values: MEASURES
                    .iter()
                    .map(|measure| (*measure, values[*measure]))
                    .collect(),
            })
            .collect::<Vec<_>>();
        assert!(verdict(&rows, &[], &base, 0.0).is_empty());
        let mut worse = rows.clone();
        if let Some(value) = worse[0].values.get_mut("dominant_right") {
            *value -= 2.0;
        }
        if let Some(value) = worse[1].values.get_mut("clear_wrong") {
            *value += 2.0;
        }
        let failures = verdict(&worse, &[], &base, 0.0);
        assert_eq!(failures.len(), 2, "{failures:?}");
        assert!(verdict(&worse, &[], &base, 2.0).is_empty());
        let mut better = rows.clone();
        if let Some(value) = better[0].values.get_mut("dominant_right") {
            *value += 2.0;
        }
        assert!(verdict(&better, &[], &base, 0.0).is_empty());
        assert_eq!(verdict(&rows[1..], &[], &base, 0.0).len(), 1);
    }
}
