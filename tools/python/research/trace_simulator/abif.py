"""Minimal ABIF reader and writer for the trace simulator (research only).

The reader returns the tags the forward-model fit needs from a real trace;
the writer emits the smallest container `dna` accepts: DATA.9-12, FWO_.1,
PLOC.2, PBAS.2 and PCON.2. Neither is a DNA contract: production decoding is
`dna_sanger::abif`.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

BASES = "ACGT"

Entry = tuple[int, int, int, int, bytes]


@dataclass
class Trace:
    """The analysed channels in A/C/G/T order, the vendor calls and metadata."""

    channels: np.ndarray  # shape (4, samples), float64, A/C/G/T
    ploc: np.ndarray  # vendor call positions (samples)
    pbas: str  # vendor calls
    pcon: np.ndarray  # vendor qualities
    meta: dict[str, str] = field(default_factory=dict)
    raw: np.ndarray | None = None  # DATA.1-4 in A/C/G/T order, when read


def _entries(data: bytes) -> dict[tuple[str, int], Entry]:
    if data[:4] != b"ABIF":
        raise ValueError("not an ABIF container")
    count, offset = (
        struct.unpack(">I", data[18:22])[0],
        struct.unpack(">I", data[26:30])[0],
    )
    entries: dict[tuple[str, int], Entry] = {}
    for index in range(count):
        base = offset + index * 28
        tag, number, etype, esize, n, dsize, doff = struct.unpack(
            ">4sIHHIII", data[base : base + 24]
        )
        payload = (
            data[base + 20 : base + 20 + dsize]
            if dsize <= 4
            else data[doff : doff + dsize]
        )
        entries[(tag.decode("latin-1"), number)] = (etype, esize, n, dsize, payload)
    return entries


def _text(entry: Entry) -> str:
    etype, _, _, _, payload = entry
    if etype == 18 and payload:  # pString: length byte first
        return payload[1 : 1 + payload[0]].decode("latin-1", "replace")
    return payload.decode("latin-1", "replace").rstrip("\x00")


def _ordered(
    entries: dict[tuple[str, int], Entry], first: int, order: str
) -> np.ndarray:
    rows = [
        np.frombuffer(entries[("DATA", first + index)][4], dtype=">i2").astype(
            np.float64
        )
        for index in range(4)
    ]
    channels = np.zeros((4, len(rows[0])))
    for index, base in enumerate(order):
        channels[BASES.index(base)] = rows[index]
    return channels


META_TAGS = {
    "instrument": ("HCFG", 4),
    "machine": ("MCHN", 1),
    "run": ("RunN", 1),
    "mix_lot": ("SMLt", 1),
    "dye_set": ("DySN", 1),
    "mobility": ("PDMF", 1),
    "basecaller": ("SPAC", 2),
}


def read(path: Path | str, raw: bool = False) -> Trace:
    """Reads the analysed channels, vendor calls and run metadata of `path`."""
    data = Path(path).read_bytes()
    entries = _entries(data)
    order = entries[("FWO_", 1)][4][:4].decode()
    channels = _ordered(entries, 9, order)
    raw_channels = _ordered(entries, 1, order) if raw else None
    ploc = np.frombuffer(entries[("PLOC", 2)][4], dtype=">i2").astype(np.int64)
    pbas = entries[("PBAS", 2)][4].decode("latin-1")
    pcon = np.frombuffer(entries[("PCON", 2)][4], dtype=np.uint8).astype(np.int64)
    meta: dict[str, str] = {}
    for key, tag in META_TAGS.items():
        if tag in entries:
            meta[key] = _text(entries[tag])
    if ("SPAC", 1) in entries:
        meta["spacing"] = f"{struct.unpack('>f', entries[('SPAC', 1)][4][:4])[0]:.3f}"
    if ("RUND", 1) in entries:
        year, month, day = struct.unpack(">HBB", entries[("RUND", 1)][4][:4])
        meta["run_date"] = f"{year:04d}-{month:02d}-{day:02d}"
    return Trace(channels, ploc, pbas, pcon, meta, raw_channels)


def write(
    path: Path | str,
    channels: np.ndarray,
    ploc: np.ndarray,
    pbas: str,
    pcon: np.ndarray,
) -> None:
    """Writes a minimal ABIF: `channels` (A/C/G/T rows) as DATA.9-12 in GATC order."""
    order = "GATC"
    clipped = np.clip(np.rint(channels), -32768, 32767).astype(">i2")
    records: list[tuple[bytes, int, int, int, bytes]] = []
    for index, base in enumerate(order):
        records.append((b"DATA", 9 + index, 4, 2, clipped[BASES.index(base)].tobytes()))
    records.append((b"FWO_", 1, 2, 1, order.encode()))
    records.append((b"PLOC", 2, 4, 2, np.asarray(ploc, dtype=">i2").tobytes()))
    records.append((b"PBAS", 2, 2, 1, pbas.encode()))
    records.append((b"PCON", 2, 2, 1, np.asarray(pcon, dtype=np.uint8).tobytes()))
    directory_offset = 128
    payload_offset = directory_offset + 28 * len(records)
    offsets = []
    for record in records:
        if len(record[4]) > 4:
            offsets.append(payload_offset)
            payload_offset += len(record[4])
        else:
            offsets.append(0)
    out = bytearray(payload_offset)
    out[0:4] = b"ABIF"
    out[4:6] = struct.pack(">H", 101)
    out[6:34] = struct.pack(
        ">4sIHHIII",
        b"tdir",
        1,
        1023,
        28,
        len(records),
        28 * len(records),
        directory_offset,
    ) + bytes(4)
    for index, (tag, number, etype, esize, payload) in enumerate(records):
        base = directory_offset + index * 28
        head = struct.pack(
            ">4sIHHII", tag, number, etype, esize, len(payload) // esize, len(payload)
        )
        if len(payload) > 4:
            tail = struct.pack(">I", offsets[index])
            out[offsets[index] : offsets[index] + len(payload)] = payload
        else:
            tail = payload.ljust(4, b"\x00")
        out[base : base + 28] = head + tail + bytes(4)
    Path(path).write_bytes(bytes(out))
