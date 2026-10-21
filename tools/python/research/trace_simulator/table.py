# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy"]
# ///
"""Reduces the fitted reads to the simulator's aggregate parameter table (research only).

Each fitted read's parameters are those `model.from_fit` gives the research
simulator. The table keeps per primer only quantiles over reads (levels 0.025 to
0.975), never a read's own values, plus the pooled context tables and the
calibration constants. `dna-trace`'s renderer (feature `simulate`) draws a read's
parameters from it, so the committed data hold aggregates only (PROP-0015
increment 1).

Usage:
    uv run --script table.py POOLED.json OUT.json [--min-span 150]
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import model  # noqa: E402

LEVELS = [round(0.025 + 0.05 * i, 3) for i in range(20)]
ENVELOPE_CALLS = list(range(0, 801, 25))
SCALARS = (
    "spacing0",
    "spacing_slope",
    "sigma0",
    "sigma_slope",
    "height_det_sd",
    "height_iid_sd",
    "position_det_sd",
    "position_iid_sd",
    "floor_sd",
    "white_sd",
)
VECTORS = ("dye_offset", "dye_gain", "floor")
CROSSTALK = [(i, j) for i in range(4) for j in range(4) if i != j]


def norm_ppf(p: np.ndarray) -> np.ndarray:
    """The standard normal quantile (Acklam's rational approximation)."""
    from statistics import NormalDist

    return np.array([NormalDist().inv_cdf(float(x)) for x in p])


def copula_columns(reads: list[model.Params]) -> dict[str, list[float]]:
    """Each read's parameters the copula ties together, by name."""
    columns: dict[str, list[float]] = {}
    for name in ("spacing0", "spacing_slope", "sigma0", "sigma_slope",
                 "height_det_sd", "position_det_sd", "floor_sd", "white_sd"):
        columns[name] = [getattr(p, name) for p in reads]
    for name in VECTORS:
        for c in range(4):
            columns[f"{name}.{c}"] = [getattr(p, name)[c] for p in reads]
    calls = np.asarray(ENVELOPE_CALLS, dtype=float)
    curves = np.array([p.log_envelope(calls) for p in reads])
    for k, call in enumerate(ENVELOPE_CALLS):
        columns[f"envelope.{call}"] = list(curves[:, k])
    for i, j in CROSSTALK:
        columns[f"crosstalk.{i}{j}"] = [p.crosstalk[i][j] for p in reads]
    return columns


def quantiles(values: list[float]) -> list[float]:
    return [round(float(v), 6) for v in np.quantile(np.asarray(values, float), LEVELS)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pooled", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--min-span", type=int, default=150)
    args = parser.parse_args()
    pooled = model.load_pooled(args.pooled)
    by_primer: dict[str, list[model.Params]] = {}
    for trace in pooled["traces"]:
        if trace["span"][1] < args.min_span:
            continue
        by_primer.setdefault(trace["primer"], []).append(model.from_fit(trace, pooled))
    primers = {}
    for primer, reads in sorted(by_primer.items()):
        entry: dict = {"reads": len(reads)}
        for name in SCALARS:
            entry[name] = quantiles([getattr(p, name) for p in reads])
        for name in VECTORS:
            entry[name] = [
                quantiles([getattr(p, name)[c] for p in reads]) for c in range(4)
            ]
        calls = np.asarray(ENVELOPE_CALLS, dtype=float)
        curves = np.array([p.log_envelope(calls) for p in reads])
        entry["envelope"] = [quantiles(list(curves[:, k])) for k in range(len(calls))]
        # The dependence between a read's parameters, as a Gaussian copula:
        # the correlation of the parameters' normal scores (ranks) over reads.
        columns = copula_columns(reads)
        scores = []
        for values in columns.values():
            ranks = np.argsort(np.argsort(values, kind="stable"), kind="stable")
            scores.append(norm_ppf((ranks + 0.5) / len(values)))
        correlation = np.corrcoef(np.array(scores))
        entry["copula"] = {
            "parameters": list(columns),
            "correlation": [[round(float(v), 4) for v in row] for row in correlation],
        }
        entry["crosstalk"] = {
            f"{i}{j}": quantiles([p.crosstalk[i][j] for p in reads]) for i, j in CROSSTALK
        }
        primers[primer] = entry
    defaults = model.Params()
    sample = next(iter(by_primer.values()))[0]
    table = {
        "version": 1,
        "provenance": {
            "fit": "tools/python/research/trace_simulator fit.py, dev v5, 300 samples, seed 20261021",
            "reads": sum(len(v) for v in by_primer.values()),
            "min_span": args.min_span,
            "levels": "quantiles over reads; no read's own values",
        },
        "levels": LEVELS,
        "envelope_calls": ENVELOPE_CALLS,
        "first_sample": sample.first_sample,
        "minor_mu": sample.minor_mu,
        "minor_sd": sample.minor_sd,
        "floor_len": sample.floor_len,
        "position_smooth": sample.position_smooth,
        "start_calls": defaults.start_calls,
        "start_noise": defaults.start_noise,
        "qv_share": list(model.QV_SHARE),
        "qv_value": list(model.QV_VALUE),
        "height_context": {k: round(v, 6) for k, v in sorted(sample.height_context.items())},
        "position_context": {
            k: round(v, 6) for k, v in sorted(sample.position_context.items())
        },
        "primers": primers,
    }
    assert all(
        f.name in {*SCALARS, *VECTORS} or f.name in table or f.name in {
            "envelope", "crosstalk", "height_context", "position_context", "seed",
            "pattern_seed",
        }
        for f in dataclasses.fields(model.Params)
    )
    # One line per list of numbers, so the committed table diffs by parameter.
    text = json.dumps(table, indent=1)
    text = re.sub(
        r"\[\s*(-?[0-9.e+-]+(?:,\s*-?[0-9.e+-]+)*)\s*\]",
        lambda m: "[" + ", ".join(v.strip() for v in m.group(1).split(",")) + "]",
        text,
    )
    args.out.write_text(text + "\n")
    print(f"wrote {args.out}: " + ", ".join(f"{k} {v['reads']}" for k, v in primers.items()))


if __name__ == "__main__":
    main()
