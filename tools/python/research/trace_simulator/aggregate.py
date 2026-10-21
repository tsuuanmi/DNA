# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy"]
# ///
"""Reduces `fit.py` output to parameter distributions and pooled tables (research only).

Usage:
    uv run --script aggregate.py FIT.jsonl OUT_DIR

Writes `parameters.md` (q10/median/q90 per group: all, primer, instrument, period)
and `pooled.json` (context tables, per-trace parameter sets for the simulator).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reference  # noqa: E402

BASES = "ACGT"


def q(values: list[float]) -> str:
    v = np.array([x for x in values if x is not None and math.isfinite(x)], dtype=float)
    if len(v) == 0:
        return "-"
    lo, mid, hi = np.percentile(v, [10, 50, 90])
    return f"{lo:.3g} / {mid:.3g} / {hi:.3g}"


def ratio(record: dict, key: str, at: str) -> float | None:
    width = record.get(key, {}).get(at)
    spacing = record.get("spacing_at", {}).get(at)
    return width / spacing if width and spacing else None


def white_rel(record: dict) -> float | None:
    white = record.get("floor_white_sd")
    height = record.get("height_median")
    if not white or not height:
        return None
    finite = [w for w in white if w is not None and math.isfinite(w)]
    return float(np.median(finite)) / height if finite else None


def floor_rel(record: dict) -> float | None:
    sd = record.get("floor_sd")
    height = record.get("height_median")
    if not sd or not height:
        return None
    return float(np.median(sd)) / height


def env_slope(record: dict) -> float | None:
    """Log-height change per 100 calls from call 100 to 250 (signal decay)."""
    env = record.get("envelope_poly")
    span = record.get("span")
    if not env or not span or span[1] < 250 or span[0] > 100:
        return None
    a = np.polyval(env, 0.100)
    b = np.polyval(env, 0.250)
    return float((b - a) / 1.5)


ROWS: list[tuple[str, callable]] = [
    ("spacing at call 100 (samples)", lambda r: r.get("spacing_at", {}).get("100")),
    ("spacing at call 300 (samples)", lambda r: r.get("spacing_at", {}).get("300")),
    (
        "spacing drift, 100->300 (samples per 100 calls)",
        lambda r: (r["spacing_at"]["300"] - r["spacing_at"]["100"]) / 2
        if "300" in r.get("spacing_at", {}) and "100" in r.get("spacing_at", {})
        else None,
    ),
    (
        "Gaussian sigma at call 100 (samples)",
        lambda r: r.get("width_at", {}).get("100"),
    ),
    ("sigma / spacing at call 100", lambda r: ratio(r, "width_at", "100")),
    ("sigma / spacing at call 300", lambda r: ratio(r, "width_at", "300")),
    ("peak asymmetry (R-L)/(R+L) at half max", lambda r: r.get("asymmetry")),
    ("EMG tau (samples)", lambda r: r.get("emg_tau")),
    ("EMG residual gain over Gaussian", lambda r: r.get("emg_rss_gain")),
    (
        "dye offset spread, max-min (samples)",
        lambda r: max(r["dye_offset"]) - min(r["dye_offset"])
        if r.get("dye_offset")
        else None,
    ),
    ("dye offset A (samples)", lambda r: r.get("dye_offset", [None] * 4)[0]),
    ("dye offset C (samples)", lambda r: r.get("dye_offset", [None] * 4)[1]),
    ("dye offset G (samples)", lambda r: r.get("dye_offset", [None] * 4)[2]),
    ("dye offset T (samples)", lambda r: r.get("dye_offset", [None] * 4)[3]),
    ("peak position residual SD (samples)", lambda r: r.get("position_jitter_sd")),
    ("position residual lag-1 autocorrelation", lambda r: r.get("position_resid_lag1")),
    ("vendor PLOC - centre SD (samples)", lambda r: r.get("vendor_jitter_sd")),
    ("envelope peak call", lambda r: r.get("envelope_peak_call")),
    ("envelope peak height (RFU)", lambda r: r.get("envelope_peak_height")),
    ("log-height slope, calls 100-250 (per 100 calls)", env_slope),
    ("dye gain A (log)", lambda r: r.get("dye_gain", [None] * 4)[0]),
    ("dye gain C (log)", lambda r: r.get("dye_gain", [None] * 4)[1]),
    ("dye gain G (log)", lambda r: r.get("dye_gain", [None] * 4)[2]),
    ("dye gain T (log)", lambda r: r.get("dye_gain", [None] * 4)[3]),
    ("log-height residual SD", lambda r: r.get("height_resid_sd")),
    ("log-height residual lag-1", lambda r: r.get("height_resid_lag1")),
    ("off-peak floor SD / median height", floor_rel),
    ("white noise SD / median height", white_rel),
    ("off-peak floor lag-1 (samples)", lambda r: r.get("floor_lag1")),
    ("samples at or below zero (fraction)", lambda r: r.get("floor_zero_fraction")),
    (
        "clean-call second/primary, median",
        lambda r: r.get("second_ratio_q", [None] * 3)[1],
    ),
    (
        "clean-call second/primary, q90",
        lambda r: r.get("second_ratio_q", [None] * 3)[2],
    ),
    ("clean calls", lambda r: r.get("clean_calls")),
]


def table(groups: dict[str, list[dict]], title: str) -> list[str]:
    names = list(groups)
    lines = [
        f"### {title}",
        "",
        "| parameter | "
        + " | ".join(f"{n} (n={len(groups[n])})" for n in names)
        + " |",
    ]
    lines.append("|---|" + "---|" * len(names))
    for label, fn in ROWS:
        cells = []
        for name in names:
            values = []
            for record in groups[name]:
                try:
                    values.append(fn(record))
                except (KeyError, TypeError, IndexError):
                    values.append(None)
            cells.append(q(values))
        lines.append(f"| {label} | " + " | ".join(cells) + " |")
    lines.append("")
    return lines


def pooled_context(records: list[dict], key: str) -> dict[str, list[float]]:
    sums: dict[str, list[float]] = defaultdict(lambda: [0.0, 0.0, 0])
    for record in records:
        for context, (s, ss, n) in record.get(key, {}).items():
            cell = sums[context]
            cell[0] += s
            cell[1] += ss
            cell[2] += n
    out = {}
    for context, (s, ss, n) in sums.items():
        if n >= 30:
            mean = s / n
            out[context] = [mean, math.sqrt(max(ss / n - mean * mean, 0)), n]
    return out


def explained(records: list[dict], key: str, table_: dict[str, list[float]]) -> float:
    """Share of the residual variance the pooled context means explain."""
    total = between = count = 0.0
    for record in records:
        for context, (s, ss, n) in record.get(key, {}).items():
            total += ss
            count += n
            if context in table_:
                between += n * table_[context][0] ** 2
    return between / total if total else float("nan")


def reproducibility(
    records: list[dict], key: str, table_: dict[str, list[float]]
) -> dict:
    """Intraclass correlation of the context-free residual at one reference position.

    Calls 20-180 of HV2F (forward) and HV1R (reverse) reads, mapped to rCRS by the
    trace's placement (no indel assumed so early in these reads). The share of the
    residual left after the trinucleotide table that recurs at the same position in
    other samples' reads is sequence-determined; the rest is per-read noise.
    """
    text = reference.rcrs()
    groups: dict[tuple[str, int], list[float]] = defaultdict(list)
    for record in records:
        placement = record.get("placement")
        if not placement or record.get("primer") not in ("HV2F", "HV1R"):
            continue
        strand, start = placement["strand"], placement["start"]
        for i, value in record.get(key, []):
            if not 20 <= i <= 180:
                continue
            pos = start + i if strand == "F" else start - i
            if not 2 <= pos <= len(text) - 1:
                continue
            tri = text[pos - 2 : pos + 1]
            if strand == "R":
                tri = reference.revcomp(tri)
            mean = table_.get(tri, [0.0])[0]
            groups[(strand, pos)].append(value - mean)
    data = [np.array(v) for v in groups.values() if len(v) >= 5]
    if len(data) < 2:
        return {}
    sizes = [len(v) for v in data]
    grand = np.concatenate(data).mean()
    k = len(data)
    n_total = sum(sizes)
    msb = sum(len(v) * (v.mean() - grand) ** 2 for v in data) / (k - 1)
    msw = sum(((v - v.mean()) ** 2).sum() for v in data) / (n_total - k)
    n0 = (n_total - sum(s * s for s in sizes) / n_total) / (k - 1)
    icc = (msb - msw) / (msb + (n0 - 1) * msw)
    total_sd = float(np.concatenate(data).std())
    return {
        "positions": k,
        "values": n_total,
        "icc": float(icc),
        "residual_sd": total_sd,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fit", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    records = [json.loads(line) for line in args.fit.open()]
    ok = [r for r in records if r.get("status") == "ok"]
    args.out.mkdir(parents=True, exist_ok=True)
    lines = [
        f"# Forward-model parameters ({len(ok)} of {len(records)} traces fitted)",
        "",
        "Cells are q10 / median / q90 over traces.",
        "",
    ]
    lines += table({"all": ok}, "All traces")
    by_primer: dict[str, list[dict]] = defaultdict(list)
    by_inst: dict[str, list[dict]] = defaultdict(list)
    by_period: dict[str, list[dict]] = defaultdict(list)
    for r in ok:
        by_primer[r.get("primer") or "?"].append(r)
        by_inst[r.get("instrument") or "?"].append(r)
        by_period[
            "2026" if r.get("batch_date", "") >= "2026-01" else "before 2026"
        ].append(r)
    lines += table(dict(sorted(by_primer.items())), "By primer")
    lines += table(dict(sorted(by_period.items())), "By period (batch date)")
    big = {
        k: v
        for k, v in sorted(by_inst.items(), key=lambda kv: -len(kv[1]))
        if len(v) >= 20
    }
    lines += table(big, "By instrument (serial, at least 20 traces)")
    height_ctx = pooled_context(ok, "context")
    position_ctx = pooled_context(ok, "position_context")
    lines.append("### Context")
    lines.append("")
    lines.append(
        f"- Height: {len(height_ctx)} trinucleotide cells (n >= 30) explain "
        f"{explained(ok, 'context', height_ctx):.1%} of the log-height residual variance."
    )
    lines.append(
        f"- Position: {len(position_ctx)} cells explain "
        f"{explained(ok, 'position_context', position_ctx):.1%} of the position residual variance."
    )
    pooled_rep: dict[str, dict] = {}
    for label, key, table_ in (
        ("Height", "height_resid_calls", height_ctx),
        ("Position", "position_resid_calls", position_ctx),
    ):
        rep = reproducibility(ok, key, table_)
        if rep:
            lines.append(
                f"- {label} residual after the context table, same rCRS position across "
                f"samples (HV2F, HV1R calls 20-180; {rep['positions']} positions, "
                f"{rep['values']} values): SD {rep['residual_sd']:.3f}, intraclass "
                f"correlation {rep['icc']:.2f} (the sequence-determined share)."
            )
        pooled_rep[label.lower()] = rep
    ranked = sorted(height_ctx.items(), key=lambda kv: kv[1][0])
    lines.append(
        "- Weakest heights (prev-own-next, mean log residual): "
        + ", ".join(f"{k} {v[0]:+.2f}" for k, v in ranked[:6])
    )
    lines.append(
        "- Strongest heights: "
        + ", ".join(f"{k} {v[0]:+.2f}" for k, v in ranked[-6:][::-1])
    )
    di: dict[str, list[float]] = defaultdict(lambda: [0.0, 0])
    for k, (mean, _sd, n) in height_ctx.items():
        di[k[:2]][0] += mean * n
        di[k[:2]][1] += n
    lines.append(
        "- By previous and own base (mean log residual): "
        + ", ".join(f"{k} {s / n:+.2f}" for k, (s, n) in sorted(di.items()))
    )
    lines.append("")
    (args.out / "parameters.md").write_text("\n".join(lines))
    keep = [
        "path",
        "primer",
        "instrument",
        "batch",
        "batch_date",
        "placement",
        "calls",
        "span",
        "spacing_at",
        "width_at",
        "width_line",
        "dye_offset",
        "position_jitter_sd",
        "vendor_jitter_sd",
        "envelope_poly",
        "envelope_at",
        "dye_gain",
        "height_resid_sd",
        "height_resid_lag1",
        "floor",
        "floor_sd",
        "floor_white_sd",
        "floor_lag1",
        "crosstalk",
        "height_median",
        "emg_tau",
        "qv20_first",
        "qv20_last",
    ]
    pooled = {
        "height_context": height_ctx,
        "position_context": position_ctx,
        "reproducibility": pooled_rep,
        "traces": [{k: r.get(k) for k in keep} for r in ok],
    }
    (args.out / "pooled.json").write_text(json.dumps(pooled))
    print(f"wrote {args.out}/parameters.md and pooled.json")


if __name__ == "__main__":
    main()
