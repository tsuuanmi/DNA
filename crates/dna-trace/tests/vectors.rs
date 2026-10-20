//! The tract decoder against the research reference's test vectors
//! (PROP-0012 increment 2): ten HVS-II reads decoded by `ref_decoder.py`,
//! the numpy reference of the decoder spike. The downstream sequence and the
//! components must match exactly, the numbers to a QR port's tolerance.

use std::path::Path;

use dna_trace::tract::{
    LeakageParams, Marginal, ONE_HOT, TractParams, choose, decode, leakage, refit, signature,
};
use serde_json::Value;

const BASES: &str = "ACGT";
/// The calls of in-phase context the reference extracted before each tract.
const CONTEXT: usize = 200;

fn base(value: &Value) -> usize {
    BASES
        .find(value.as_str().unwrap_or_default())
        .unwrap_or_default()
}

fn index(value: &Value) -> Option<usize> {
    value.as_u64().and_then(|value| usize::try_from(value).ok())
}

fn mark(value: &Value) -> Option<usize> {
    value.as_str().and_then(|mark| BASES.find(mark))
}

fn number(value: &Value) -> f64 {
    value.as_f64().unwrap_or(f64::NAN)
}

fn close(actual: f64, expected: f64, what: &str, id: &str) {
    let tolerance = 1e-8 * expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance,
        "{id}: {what} {actual} != {expected}"
    );
}

fn check_marginal(actual: &Marginal, expected: &Value, what: &str, id: &str) {
    let keys = expected["keys"]
        .as_array()
        .map(|keys| keys.iter().filter_map(index).collect::<Vec<_>>())
        .unwrap_or_default();
    assert_eq!(actual.lengths, keys, "{id}: {what} keys");
    for (index, share) in expected["share"]
        .as_array()
        .into_iter()
        .flatten()
        .enumerate()
    {
        close(
            actual.shares[index],
            number(share),
            &format!("{what} share {index}"),
            id,
        );
    }
    for (index, error) in expected["se"].as_array().into_iter().flatten().enumerate() {
        close(
            actual.errors[index],
            number(error),
            &format!("{what} se {index}"),
            id,
        );
    }
    assert_eq!(
        actual.dominant,
        index(&expected["mode"]),
        "{id}: {what} mode"
    );
    assert_eq!(
        actual.next,
        index(&expected["second"]),
        "{id}: {what} second"
    );
    let z = number(&expected["z_raw"]);
    if z >= 1e9 {
        assert!(actual.z.is_infinite(), "{id}: {what} single length");
    } else {
        close(actual.z, z, &format!("{what} z"), id);
        close(
            actual.z_adjusted,
            number(&expected["z_adj"]),
            &format!("{what} z_adj"),
            id,
        );
    }
}

fn profiles(input: &Value) -> Vec<[f64; 4]> {
    input["Y"]
        .as_array()
        .into_iter()
        .flatten()
        .map(|row| {
            let mut profile = [0.0; 4];
            for (cell, value) in profile.iter_mut().zip(row.as_array().into_iter().flatten()) {
                *cell = number(value);
            }
            profile
        })
        .collect()
}

fn tract_params(params: &Value) -> TractParams {
    let count = |name: &str| index(&params[name]).unwrap_or_default();
    TractParams {
        span_below: count("span_below"),
        span_above: count("span_above"),
        min_length: count("min_length"),
        downstream: count("M"),
        joint_lengths: count("KMAX"),
        iterations: count("iters"),
    }
}

fn vector_files(directory: &Path) -> Result<Vec<std::path::PathBuf>, std::io::Error> {
    let mut files = std::fs::read_dir(directory)?
        .filter_map(|entry| entry.ok().map(|entry| entry.path()))
        .filter(|path| {
            path.extension()
                .is_some_and(|extension| extension == "json")
        })
        .collect::<Vec<_>>();
    files.sort();
    Ok(files)
}

// SRS-CALL-025
#[test]
fn decodes_the_reference_vectors() -> Result<(), Box<dyn std::error::Error>> {
    let directory = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/vectors");
    let files = vector_files(&directory)?;
    assert_eq!(files.len(), 10);
    for file in files {
        let vector: Value = serde_json::from_str(&std::fs::read_to_string(&file)?)?;
        let id = vector["id"].as_str().unwrap_or_default().to_owned();
        let input = &vector["input"];
        let profiles = profiles(input);
        let params = &vector["params"];
        let count = |name: &str| index(&params[name]).unwrap_or_default();
        let tract = tract_params(params);
        let decoded = decode(
            &profiles,
            index(&input["e"]).unwrap_or_default(),
            base(&input["run"]),
            Some(base(&input["mark"])),
            None,
            &ONE_HOT,
            &tract,
        )
        .ok_or("no decode")?;
        let derived = &vector["derived"];
        assert_eq!(
            decoded.extent as u64,
            derived["X"].as_u64().unwrap_or_default(),
            "{id}: X"
        );
        assert_eq!(
            decoded.window as u64,
            derived["W"].as_u64().unwrap_or_default(),
            "{id}: W"
        );
        let expected = &vector["decode"];
        let downstream = expected["D"]
            .as_array()
            .into_iter()
            .flatten()
            .filter_map(index)
            .collect::<Vec<_>>();
        assert_eq!(decoded.downstream, downstream, "{id}: D");
        close(decoded.residual, number(&expected["res"]), "res", &id);
        close(decoded.sigma, number(&expected["sig_fit"]), "sig_fit", &id);
        close(
            decoded.autocorrelation,
            number(&expected["rho1"]),
            "rho",
            &id,
        );
        check_marginal(&decoded.first, &expected["first"], "first", &id);
        check_marginal(&decoded.total, &expected["total"], "total", &id);
        // SL: the read's own signature and leakage null, a signature decode,
        // and the refit with leakage over its downstream sequence.
        let sl = &vector["sl"];
        let start = index(&input["e"]).unwrap_or_default();
        let gap = count("signature_gap");
        let null = LeakageParams {
            gap,
            min_run: count("leak_minrun"),
            pad: count("leak_pad"),
            min_rows: count("leak_min_rows"),
            context: CONTEXT,
        };
        let rows = signature(&profiles, start, &null);
        for (row, wanted) in rows
            .iter()
            .zip(sl["signature"].as_array().into_iter().flatten())
        {
            for (value, cell) in row.iter().zip(wanted.as_array().into_iter().flatten()) {
                close(*value, number(cell), "signature", &id);
            }
        }
        let shares = leakage(&profiles, start, &rows, &null);
        match (&shares, sl["leak"].as_object()) {
            (Some(shares), Some(wanted)) => {
                for (offset, share) in shares {
                    close(*share, number(&wanted[&offset.to_string()]), "leak", &id);
                }
            }
            (None, None) => {}
            _ => panic!("{id}: leak presence differs"),
        }
        let run = base(&input["run"]);
        let mark = Some(base(&input["mark"]));
        let signed = decode(&profiles, start, run, mark, None, &rows, &tract)
            .ok_or("no signature decode")?;
        let wanted_downstream = sl["signature_decode"]["D"]
            .as_array()
            .into_iter()
            .flatten()
            .filter_map(index)
            .collect::<Vec<_>>();
        assert_eq!(signed.downstream, wanted_downstream, "{id}: signature D");
        close(
            signed.residual,
            number(&sl["signature_decode"]["res"]),
            "signature res",
            &id,
        );
        let refitted = refit(
            &profiles,
            start,
            run,
            mark,
            &rows,
            shares.as_deref().unwrap_or(&[]),
            &signed.downstream,
            &tract,
        )
        .ok_or("no refit")?;
        close(
            refitted.residual,
            number(&sl["refit"]["res"]),
            "refit res",
            &id,
        );
        close(
            refitted.autocorrelation,
            number(&sl["refit"]["rho1"]),
            "refit rho",
            &id,
        );
        check_marginal(&refitted.first, &sl["refit"]["first"], "refit first", &id);
    }
    Ok(())
}

/// The structure choice against the research reference's three vectors:
/// every candidate's residual, decoded sequence, and marginal modes over the
/// shared window, and the choice.
// SRS-CALL-025
#[test]
fn chooses_the_reference_structures() -> Result<(), Box<dyn std::error::Error>> {
    let directory = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/vectors/structure");
    let files = vector_files(&directory)?;
    assert_eq!(files.len(), 3);
    for file in files {
        let vector: Value = serde_json::from_str(&std::fs::read_to_string(&file)?)?;
        let id = vector["id"].as_str().unwrap_or_default().to_owned();
        let input = &vector["input"];
        let profiles = profiles(input);
        let tract = tract_params(&vector["params"]);
        let start = index(&input["e"]).unwrap_or_default();
        let run = base(&input["run"]);
        let window = index(&vector["shared_W"]);
        let chosen = choose(&profiles, start, run, &ONE_HOT, &tract).ok_or("no choice")?;
        assert_eq!(chosen.mark, mark(&vector["chosen"]), "{id}: chosen");
        let candidates = vector["candidates"].as_array().ok_or("no candidates")?;
        assert_eq!(chosen.residuals.len(), candidates.len(), "{id}: candidates");
        for ((candidate_mark, residual), expected) in chosen.residuals.iter().zip(candidates) {
            let wanted = mark(&expected["mark"]);
            assert_eq!(*candidate_mark, wanted, "{id}: candidate order");
            close(*residual, number(&expected["res"]), "candidate res", &id);
            let decoded = decode(&profiles, start, run, wanted, window, &ONE_HOT, &tract)
                .ok_or("no candidate decode")?;
            assert_eq!(
                decoded.extent as u64,
                expected["X"].as_u64().unwrap_or_default(),
                "{id}: X"
            );
            assert_eq!(decoded.window, chosen.decode.window, "{id}: shared W");
            let downstream = expected["D"]
                .as_array()
                .into_iter()
                .flatten()
                .filter_map(index)
                .collect::<Vec<_>>();
            assert_eq!(decoded.downstream, downstream, "{id}: D");
            for (marginal, name) in [(&decoded.first, "first"), (&decoded.total, "total")] {
                assert_eq!(
                    marginal.dominant,
                    index(&expected[name]["mode"]),
                    "{id}: {name} mode"
                );
                for (share, wanted) in marginal
                    .shares
                    .iter()
                    .zip(expected[name]["share"].as_array().into_iter().flatten())
                {
                    close(*share, number(wanted), name, &id);
                }
            }
            if wanted == chosen.mark {
                assert_eq!(decoded, chosen.decode, "{id}: chosen decode");
            }
        }
    }
    Ok(())
}

/// The locator, from the call after each reference tract's run, finds the
/// tract's first call and run base. The decode vectors' tracts open where
/// the reference put them; a structure vector's hand-placed run may extend
/// back over a single-call mark (`AA C AAAAAA` opens two calls earlier).
// SRS-CALL-025
#[test]
fn locates_the_reference_tracts() -> Result<(), Box<dyn std::error::Error>> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/vectors");
    let params = dna_trace::locate::LocateParams {
        min_run: 4,
        min_before: 2,
    };
    for (directory, exact) in [(root.clone(), true), (root.join("structure"), false)] {
        for file in vector_files(&directory)? {
            let vector: Value = serde_json::from_str(&std::fs::read_to_string(&file)?)?;
            let id = vector["id"].as_str().unwrap_or_default().to_owned();
            let input = &vector["input"];
            let primaries = profiles(input)
                .iter()
                .map(|profile| {
                    (1..4).fold(0, |best, channel| {
                        if profile[channel] > profile[best] {
                            channel
                        } else {
                            best
                        }
                    })
                })
                .collect::<Vec<_>>();
            let start = index(&input["e"]).unwrap_or_default();
            let run = base(&input["run"]);
            let extent = primaries[start..]
                .iter()
                .take_while(|&&call| call == run)
                .count();
            let (located, base) =
                dna_trace::locate::locate(&primaries, start + extent, &params).ok_or("none")?;
            assert_eq!(base, run, "{id}: run");
            if exact {
                assert_eq!(located, start, "{id}: start");
            } else {
                assert!(located <= start, "{id}: start");
            }
        }
    }
    Ok(())
}

/// A read that ends within the candidate lengths: the decode fits only the
/// lengths its window holds, never reading past the read's end, and finds
/// nothing when the window holds no candidate.
// SRS-CALL-025
#[test]
fn fits_only_the_lengths_a_short_window_holds() -> Result<(), Box<dyn std::error::Error>> {
    let file = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/vectors/01_easy_F_clean.json");
    let vector: Value = serde_json::from_str(&std::fs::read_to_string(file)?)?;
    let input = &vector["input"];
    let tract = tract_params(&vector["params"]);
    let start = index(&input["e"]).unwrap_or_default();
    let run = base(&input["run"]);
    let full = profiles(input);
    for tail in [9, 5] {
        let short = &full[..start + tail];
        let decoded = decode(
            short,
            start,
            run,
            Some(base(&input["mark"])),
            None,
            &ONE_HOT,
            &tract,
        )
        .ok_or("no decode")?;
        assert_eq!(decoded.window, tail);
        assert!(decoded.total.lengths.iter().all(|&length| length <= tail));
        let chosen = choose(short, start, run, &ONE_HOT, &tract).ok_or("no choice")?;
        assert!(
            chosen
                .decode
                .total
                .lengths
                .iter()
                .all(|&length| length <= tail)
        );
    }
    let shortest = &full[..start + 2];
    assert!(decode(shortest, start, run, None, None, &ONE_HOT, &tract).is_none());
    Ok(())
}
