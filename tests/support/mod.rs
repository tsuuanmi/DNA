//! Synthetic ABIF, FASTA, and configuration fixtures shared by integration tests.

//!
//! Every test binary declares `pub mod support` and these helpers are `pub`, so a
//! binary that uses only some of them reports no dead code.

use std::fs;
use std::path::{Path, PathBuf};

/// Writes a clean synthetic ABIF trace whose vendor calls equal `sequence`.
///
/// # Errors
///
/// Fails when the fixture parameters are inconsistent or the file cannot be written.
pub fn write_abif(path: &Path, sequence: &str) -> Result<(), Box<dyn std::error::Error>> {
    write_abif_fixture(
        path,
        sequence,
        sequence.as_bytes(),
        1,
        *b"ACGT",
        None,
        None,
        None,
    )
}

/// Writes a synthetic ABIF trace with one primary peak height per call.
///
/// # Errors
///
/// Fails when the fixture parameters are inconsistent or the file cannot be written.
pub fn write_abif_with_peak_heights(
    path: &Path,
    sequence: &str,
    peak_heights: Vec<i16>,
) -> Result<(), Box<dyn std::error::Error>> {
    write_abif_fixture(
        path,
        sequence,
        sequence.as_bytes(),
        1,
        *b"ACGT",
        None,
        None,
        Some(peak_heights),
    )
}

/// Writes a synthetic ABIF trace with separate vendor base calls and PCON element type.
///
/// # Errors
///
/// Fails when the fixture parameters are inconsistent or the file cannot be written.
pub fn write_abif_with_vendor(
    path: &Path,
    sequence: &str,
    vendor_primary: &str,
    pcon_element_type: u16,
) -> Result<(), Box<dyn std::error::Error>> {
    if vendor_primary.len() != sequence.len() {
        return Err("synthetic vendor sequence length must equal signal sequence length".into());
    }
    write_abif_fixture(
        path,
        sequence,
        vendor_primary.as_bytes(),
        pcon_element_type,
        *b"ACGT",
        None,
        None,
        None,
    )
}

/// Writes a synthetic ABIF trace with a custom `FWO_` channel order.
///
/// # Errors
///
/// Fails when the fixture parameters are inconsistent or the file cannot be written.
pub fn write_abif_with_channel_order(
    path: &Path,
    sequence: &str,
    channel_order: [u8; 4],
) -> Result<(), Box<dyn std::error::Error>> {
    write_abif_fixture(
        path,
        sequence,
        sequence.as_bytes(),
        1,
        channel_order,
        None,
        None,
        None,
    )
}

/// Writes a synthetic ABIF trace with explicit `PLOC.2` peak locations.
///
/// # Errors
///
/// Fails when the fixture parameters are inconsistent or the file cannot be written.
pub fn write_abif_with_ploc(
    path: &Path,
    sequence: &str,
    ploc: Vec<usize>,
) -> Result<(), Box<dyn std::error::Error>> {
    write_abif_fixture(
        path,
        sequence,
        sequence.as_bytes(),
        1,
        *b"ACGT",
        Some(ploc),
        None,
        None,
    )
}

/// Writes a synthetic ABIF trace carrying an unused `P2BA` record.
///
/// # Errors
///
/// Fails when the fixture parameters are inconsistent or the file cannot be written.
pub fn write_abif_with_unused_p2ba(
    path: &Path,
    sequence: &str,
    p2ba: Vec<u8>,
) -> Result<(), Box<dyn std::error::Error>> {
    write_abif_fixture(
        path,
        sequence,
        sequence.as_bytes(),
        1,
        *b"ACGT",
        None,
        Some(p2ba),
        None,
    )
}

/// Writes a synthetic ABIF trace whose vendor `PBAS` is one base short.
///
/// # Errors
///
/// Fails when the fixture parameters are inconsistent or the file cannot be written.
pub fn write_abif_with_short_pbas(
    path: &Path,
    sequence: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let short = &sequence.as_bytes()[..sequence.len() - 1];
    write_abif_fixture(path, sequence, short, 1, *b"ACGT", None, None, None)
}

#[expect(
    clippy::too_many_arguments,
    reason = "synthetic ABIF fixtures expose every record field the tests vary"
)]
fn write_abif_fixture(
    path: &Path,
    sequence: &str,
    vendor_primary: &[u8],
    pcon_element_type: u16,
    channel_order: [u8; 4],
    ploc_override: Option<Vec<usize>>,
    p2ba: Option<Vec<u8>>,
    peak_heights: Option<Vec<i16>>,
) -> Result<(), Box<dyn std::error::Error>> {
    write_abif_fixture_options(
        path,
        sequence,
        vendor_primary,
        pcon_element_type,
        channel_order,
        ploc_override,
        p2ba,
        peak_heights,
        None,
        &[],
    )
}

/// Adds one co-located secondary peak at `call_index`.
///
/// # Errors
///
/// Fails when the fixture parameters are inconsistent or the file cannot be written.
pub fn write_abif_with_secondary_signal(
    path: &Path,
    sequence: &str,
    call_index: usize,
    secondary_base: u8,
    height: i16,
) -> Result<(), Box<dyn std::error::Error>> {
    write_abif_with_secondary_signals(path, sequence, &[(call_index, secondary_base, height)])
}

/// Adds co-located secondary peaks as `(call index, base, height)` triples.
///
/// # Errors
///
/// Fails when the fixture parameters are inconsistent or the file cannot be written.
pub fn write_abif_with_secondary_signals(
    path: &Path,
    sequence: &str,
    secondary_signals: &[(usize, u8, i16)],
) -> Result<(), Box<dyn std::error::Error>> {
    write_abif_fixture_options(
        path,
        sequence,
        sequence.as_bytes(),
        1,
        *b"ACGT",
        None,
        None,
        None,
        None,
        secondary_signals,
    )
}

/// Adds alternating background signal on every channel across `noisy_calls`.
///
/// # Errors
///
/// Fails when the fixture parameters are inconsistent or the file cannot be written.
pub fn write_abif_with_background_noise(
    path: &Path,
    sequence: &str,
    noisy_calls: std::ops::Range<usize>,
    amplitude: i16,
) -> Result<(), Box<dyn std::error::Error>> {
    write_abif_fixture_options(
        path,
        sequence,
        sequence.as_bytes(),
        1,
        *b"ACGT",
        None,
        None,
        None,
        Some((noisy_calls, amplitude)),
        &[],
    )
}

/// Adds a slippage "shadow" ladder: from `onset_call` on, every call also
/// carries, for each `(offset, fraction)`, that fraction of the primary height
/// on the channel of the call `offset` positions away; shadows on the same
/// channel add up, and a shadow on the call's own channel is skipped.
///
/// # Errors
///
/// Fails when the fixture parameters are inconsistent or the file cannot be written.
pub fn write_abif_with_shadow_ladder(
    path: &Path,
    sequence: &str,
    onset_call: usize,
    shadows: &[(isize, f64)],
) -> Result<(), Box<dyn std::error::Error>> {
    let bases = sequence.as_bytes();
    let mut heights = std::collections::BTreeMap::<(usize, u8), f64>::new();
    for call_index in onset_call..bases.len() {
        for &(offset, fraction) in shadows {
            let Some(&shadow) = call_index
                .checked_add_signed(offset)
                .and_then(|source| bases.get(source))
            else {
                continue;
            };
            if shadow != bases[call_index] {
                *heights.entry((call_index, shadow)).or_default() += fraction;
            }
        }
    }
    let mut secondary_signals = Vec::new();
    for ((call_index, shadow), fraction) in heights {
        let height = scaled_height(fraction)?;
        if height >= 1 {
            secondary_signals.push((call_index, shadow, height));
        }
    }
    write_abif_with_secondary_signals(path, sequence, &secondary_signals)
}

/// Adds double peaks across `calls` whose secondary channel copies a
/// neighbour's primary at an offset that changes from call to call, so no
/// single slippage offset explains them.
///
/// # Errors
///
/// Fails when the fixture parameters are inconsistent or the file cannot be written.
pub fn write_abif_with_incoherent_doubles(
    path: &Path,
    sequence: &str,
    calls: std::ops::Range<usize>,
    height: i16,
) -> Result<(), Box<dyn std::error::Error>> {
    let bases = sequence.as_bytes();
    let offsets: [isize; 6] = [1, -1, 2, -2, 3, -3];
    let mut secondary_signals = Vec::new();
    for call_index in calls {
        let own = *bases
            .get(call_index)
            .ok_or("synthetic incoherent-double call index is out of range")?;
        let secondary = offsets
            .iter()
            .cycle()
            .skip(call_index % offsets.len())
            .take(offsets.len())
            .filter_map(|&offset| call_index.checked_add_signed(offset))
            .filter_map(|source| bases.get(source).copied())
            .find(|&base| base != own)
            .ok_or("no neighbouring primary differs from the call's own channel")?;
        secondary_signals.push((call_index, secondary, height));
    }
    write_abif_with_secondary_signals(path, sequence, &secondary_signals)
}

/// Scales the primary peak heights from `from_call` on by `factor`.
///
/// # Errors
///
/// Fails when the fixture parameters are inconsistent or the file cannot be written.
pub fn write_abif_with_amplitude_decay(
    path: &Path,
    sequence: &str,
    from_call: usize,
    factor: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    let decayed = scaled_height(factor)?;
    let heights = (0..sequence.len())
        .map(|index| if index >= from_call { decayed } else { 1000 })
        .collect();
    write_abif_with_peak_heights(path, sequence, heights)
}

/// The synthetic full peak height of 1000 scaled by `fraction`, rounded to the
/// nearest integer; fractions outside `0..=1` are rejected.
fn scaled_height(fraction: f64) -> Result<i16, Box<dyn std::error::Error>> {
    if !(0.0..=1.0).contains(&fraction) {
        return Err("synthetic height fraction must be within 0..=1".into());
    }
    // Thousandths are exact in i32 once the fraction is in range.
    let thousandths = (fraction * 1000.0).round().to_string().parse::<i32>()?;
    Ok(i16::try_from(thousandths)?)
}

#[expect(
    clippy::too_many_arguments,
    reason = "synthetic ABIF fixtures expose every record field the tests vary"
)]
fn write_abif_fixture_options(
    path: &Path,
    sequence: &str,
    vendor_primary: &[u8],
    pcon_element_type: u16,
    channel_order: [u8; 4],
    ploc_override: Option<Vec<usize>>,
    p2ba: Option<Vec<u8>>,
    peak_heights: Option<Vec<i16>>,
    background_noise: Option<(std::ops::Range<usize>, i16)>,
    secondary_signals: &[(usize, u8, i16)],
) -> Result<(), Box<dyn std::error::Error>> {
    let spacing = 4_usize;
    let signal_locations: Vec<usize> = (0..sequence.len())
        .map(|index| 2 + index * spacing)
        .collect();
    let sample_count = signal_locations.last().copied().unwrap_or(0) + 3;
    let peak_heights = peak_heights.unwrap_or_else(|| vec![1000; sequence.len()]);
    if peak_heights.len() != sequence.len() {
        return Err("synthetic peak-height count must equal signal sequence length".into());
    }
    let mut channels: [Vec<i16>; 4] = std::array::from_fn(|_| vec![0; sample_count]);
    if let Some((noisy_calls, amplitude)) = background_noise {
        if noisy_calls.start >= noisy_calls.end || noisy_calls.end > signal_locations.len() {
            return Err("synthetic noisy-call range is invalid".into());
        }
        if amplitude <= 0 {
            return Err("synthetic noise amplitude must be positive".into());
        }
        let sample_start = signal_locations[noisy_calls.start].saturating_sub(1);
        let sample_end = signal_locations[noisy_calls.end - 1]
            .saturating_add(2)
            .min(sample_count);
        for channel in &mut channels {
            for (offset, value) in channel[sample_start..sample_end].iter_mut().enumerate() {
                if offset % 2 == 0 {
                    *value = amplitude;
                }
            }
        }
    }
    for (index, base) in sequence.bytes().enumerate() {
        let channel = channel_index(base)?;
        channels[channel][signal_locations[index]] = peak_heights[index];
    }
    for &(call_index, secondary_base, height) in secondary_signals {
        let position = *signal_locations
            .get(call_index)
            .ok_or("synthetic secondary-signal call index is out of range")?;
        if height <= 0 {
            return Err("synthetic secondary-signal height must be positive".into());
        }
        let secondary_channel = channel_index(secondary_base)?;
        let primary_channel = channel_index(sequence.as_bytes()[call_index])?;
        if secondary_channel == primary_channel {
            return Err("synthetic secondary signal must use a different channel".into());
        }
        channels[secondary_channel][position] = height;
    }
    let mut records = Vec::new();
    for (index, base) in channel_order.iter().enumerate() {
        let channel = &channels[channel_index(*base)?];
        let mut payload = Vec::with_capacity(channel.len() * 2);
        for value in channel {
            payload.extend_from_slice(&value.to_be_bytes());
        }
        records.push(Record::new(
            *b"DATA",
            9 + u32::try_from(index)?,
            4,
            2,
            payload,
        ));
    }
    records.push(Record::new(*b"FWO_", 1, 2, 1, channel_order.to_vec()));
    let ploc_locations = ploc_override.unwrap_or(signal_locations);
    let mut ploc = Vec::with_capacity(ploc_locations.len() * 2);
    for location in &ploc_locations {
        ploc.extend_from_slice(&u16::try_from(*location)?.to_be_bytes());
    }
    records.push(Record::new(*b"PLOC", 2, 4, 2, ploc));
    records.push(Record::new(*b"PBAS", 2, 2, 1, vendor_primary.to_vec()));
    if let Some(p2ba) = p2ba {
        records.push(Record::new(*b"P2BA", 1, 2, 1, p2ba));
    }
    records.push(Record::new(
        *b"PCON",
        2,
        pcon_element_type,
        1,
        vec![40; sequence.len()],
    ));

    let directory_offset = 128_usize;
    let directory_size = records.len() * 28;
    let mut payload_offset = directory_offset + directory_size;
    for record in &mut records {
        if record.payload.len() > 4 {
            record.offset = payload_offset;
            payload_offset += record.payload.len();
        }
    }
    let mut bytes = vec![0_u8; payload_offset];
    bytes[0..4].copy_from_slice(b"ABIF");
    bytes[4..6].copy_from_slice(&101_u16.to_be_bytes());
    write_entry(
        &mut bytes,
        6,
        *b"tdir",
        1,
        1023,
        28,
        records.len(),
        directory_size,
        directory_offset,
        &[],
    )?;
    for (index, record) in records.iter().enumerate() {
        let offset = directory_offset + index * 28;
        write_entry(
            &mut bytes,
            offset,
            record.tag,
            record.number,
            record.element_type,
            record.element_size,
            record.payload.len() / record.element_size,
            record.payload.len(),
            record.offset,
            &record.payload,
        )?;
        if record.payload.len() > 4 {
            bytes[record.offset..record.offset + record.payload.len()]
                .copy_from_slice(&record.payload);
        }
    }
    fs::write(path, bytes)?;
    Ok(())
}

fn channel_index(base: u8) -> Result<usize, Box<dyn std::error::Error>> {
    match base {
        b'A' => Ok(0),
        b'C' => Ok(1),
        b'G' => Ok(2),
        b'T' => Ok(3),
        _ => Err(format!("unsupported synthetic base {}", char::from(base)).into()),
    }
}

/// Writes a single-record FASTA named `synthetic`.
///
/// # Errors
///
/// Fails when the fixture parameters are inconsistent or the file cannot be written.
pub fn write_reference(path: &Path, sequence: &str) -> Result<(), Box<dyn std::error::Error>> {
    fs::write(path, format!(">synthetic\n{sequence}\n"))?;
    Ok(())
}

/// Writes a valid strict configuration plus a synthetic target profile of the
/// given reference topology, reporting every position, next to it.
///
/// # Errors
///
/// Returns an error when either file cannot be written.
pub fn write_config(path: &Path, topology: &str) -> Result<(), Box<dyn std::error::Error>> {
    let profile = path.with_extension("profile.toml");
    fs::write(
        &profile,
        format!(
            "schema_version=1\nid='synthetic-{topology}'\n[reference]\ntopology='{topology}'\n[variant_calling]\nregions=[[1, 50000]]\n"
        ),
    )?;
    // Relative to the configuration directory, so the configuration bytes and
    // their recorded SHA-256 do not depend on the temporary directory.
    let relative = profile
        .file_name()
        .ok_or("configuration path has no file name")?;
    write_config_with_profile(path, Path::new(relative))
}

/// Writes a valid strict configuration that references `profile`, resolved
/// against the configuration directory when relative.
///
/// # Errors
///
/// Returns an error when the file cannot be written or `profile` is not UTF-8.
pub fn write_config_with_profile(
    path: &Path,
    profile: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let profile = profile.to_str().ok_or("profile path must be UTF-8")?;
    fs::write(
        path,
        format!(
            "schema_version=7\nprofile='{profile}'\n[basecalling]\nsecondary_peak_ratio=0.33\n[signal_processing]\nwindow_size_bases=10\nminimum_primary_snr=3.0\nminimum_noisy_windows=2\n[quality_control]\npenalty_window_size=10\nmax_relative_quality_score=60\n[callability]\nwindow_calls=16\nonset_defect_fraction=0.375\nexit_defect_fraction=0.125\nminimum_main_share=0.35\nmaximum_far_share=0.12\nminimum_shadow_share=0.1\nweak_amplitude_fraction=0.1\nrepeat_min_length=8\nminimum_callable_calls=20\n[alignment]\nmatch_score=3\nmismatch_score=-5\nambiguous_score=0\ngap_open_score=-10\ngap_extension_score=-4\nminimum_callable_bases=20\nminimum_identity=0.80\n[sample_reconciliation]\nminimum_comparable_bases=25\nminimum_overlap_agreement=0.50\n[variant_calling]\nmax_indel_length=50\nminimum_peak_height=150\nrelative_quality_threshold=30\nread_end_margin=0\n"
        ),
    )?;
    Ok(())
}

/// Path of the shipped human-mtDNA target profile.
#[must_use]
pub fn human_mtdna_profile() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("config/profiles/human-mtdna-rcrs.toml")
}

/// Path of the analysis result the CLI publishes for `trace` under `workdir`.
#[must_use]
pub fn analysis_output_path(workdir: &Path, trace: &Path) -> PathBuf {
    workdir
        .join("results")
        .join(format!("{}.json", trace_stem(trace)))
}

/// Path of the basecall result the CLI publishes for `trace` under `workdir`.
#[must_use]
pub fn basecall_output_path(workdir: &Path, trace: &Path) -> PathBuf {
    workdir
        .join("results")
        .join(format!("{}.basecalls.json", trace_stem(trace)))
}

fn trace_stem(trace: &Path) -> &str {
    trace
        .file_stem()
        .and_then(|value| value.to_str())
        .unwrap_or("invalid")
}

struct Record {
    tag: [u8; 4],
    number: u32,
    element_type: u16,
    element_size: usize,
    payload: Vec<u8>,
    offset: usize,
}

impl Record {
    fn new(
        tag: [u8; 4],
        number: u32,
        element_type: u16,
        element_size: usize,
        payload: Vec<u8>,
    ) -> Self {
        Self {
            tag,
            number,
            element_type,
            element_size,
            payload,
            offset: 0,
        }
    }
}

#[expect(
    clippy::too_many_arguments,
    reason = "synthetic ABIF fixtures expose every record field the tests vary"
)]
fn write_entry(
    bytes: &mut [u8],
    offset: usize,
    tag: [u8; 4],
    number: u32,
    element_type: u16,
    element_size: usize,
    element_count: usize,
    data_size: usize,
    data_offset: usize,
    inline: &[u8],
) -> Result<(), Box<dyn std::error::Error>> {
    bytes[offset..offset + 4].copy_from_slice(&tag);
    bytes[offset + 4..offset + 8].copy_from_slice(&number.to_be_bytes());
    bytes[offset + 8..offset + 10].copy_from_slice(&element_type.to_be_bytes());
    bytes[offset + 10..offset + 12].copy_from_slice(&u16::try_from(element_size)?.to_be_bytes());
    bytes[offset + 12..offset + 16].copy_from_slice(&u32::try_from(element_count)?.to_be_bytes());
    bytes[offset + 16..offset + 20].copy_from_slice(&u32::try_from(data_size)?.to_be_bytes());
    if data_size <= 4 {
        bytes[offset + 20..offset + 24].fill(0);
        bytes[offset + 20..offset + 20 + inline.len()].copy_from_slice(inline);
    } else {
        bytes[offset + 20..offset + 24].copy_from_slice(&u32::try_from(data_offset)?.to_be_bytes());
    }
    Ok(())
}
