#![allow(dead_code)]

use std::fs;
use std::path::{Path, PathBuf};

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

pub fn write_abif_with_short_pbas(
    path: &Path,
    sequence: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let short = &sequence.as_bytes()[..sequence.len() - 1];
    write_abif_fixture(path, sequence, short, 1, *b"ACGT", None, None, None)
}

#[allow(clippy::too_many_arguments)]
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
        None,
    )
}

pub fn write_abif_with_secondary_signal(
    path: &Path,
    sequence: &str,
    call_index: usize,
    secondary_base: u8,
    height: i16,
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
        Some((call_index, secondary_base, height)),
    )
}

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
        None,
    )
}

#[allow(clippy::too_many_arguments)]
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
    secondary_signal: Option<(usize, u8, i16)>,
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
    if let Some((call_index, secondary_base, height)) = secondary_signal {
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
        records.push(Record::new(*b"DATA", 9 + index as u32, 4, 2, payload));
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

pub fn write_reference(path: &Path, sequence: &str) -> Result<(), Box<dyn std::error::Error>> {
    fs::write(path, format!(">synthetic\n{sequence}\n"))?;
    Ok(())
}

pub fn write_config(path: &Path, topology: &str) -> Result<(), Box<dyn std::error::Error>> {
    fs::write(
        path,
        format!(
            "schema_version=5\n[reference]\ntopology='{topology}'\n[basecalling]\nsecondary_peak_ratio=0.33\n[signal_processing]\nwindow_size_bases=10\nminimum_primary_snr=3.0\nminimum_noisy_windows=2\n[quality_control]\ntrim_window_size=10\nbest_section_fraction=0.10\nmax_relative_quality_score=60\ntrim_stringency=7.0\nminimum_retained_bases=20\n[alignment]\nmatch_score=3\nmismatch_score=-5\nambiguous_score=0\ngap_open_score=-10\ngap_extension_score=-4\nminimum_callable_bases=20\nminimum_identity=0.80\n[sample_reconciliation]\nminimum_comparable_bases=25\nminimum_overlap_agreement=0.50\n[variant_calling]\nmax_indel_length=50\nminimum_peak_height=150\nrelative_quality_threshold=30\nregions=[[1, 50000]]\n"
        ),
    )?;
    Ok(())
}

pub fn analysis_output_path(workdir: &Path, trace: &Path) -> PathBuf {
    workdir
        .join("results")
        .join(format!("{}.json", trace_stem(trace)))
}

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

#[allow(clippy::too_many_arguments)]
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
