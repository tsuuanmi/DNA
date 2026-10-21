# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy", "scipy"]
# ///
"""Realism checks: simulated vs real reads on the same summary statistics (research only).

For a sample of fitted real traces, renders the same read (rCRS at the trace's
placement, its own fitted parameters), then measures both with the same code:

- the forward-model fit itself (`fit.fit_trace`), a round trip that also checks
  the fitter recovers the parameters it was given;
- what `dna evidence` sees over each primer's pre-tract window: the second-peak
  share at the call (q50/q90/q99), the primary's profile share, the fraction of
  calls kept in phase, and, for simulated reads, the calls that differ from the
  template.

Usage:
    uv run --script realism.py POOLED.json OUT_DIR [--per-primer N] [--calibrate]
"""

from __future__ import annotations

import argparse
import json
import math
import random
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.stats import ks_2samp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import abif  # noqa: E402
import dnarun  # noqa: E402
import fit  # noqa: E402
import model  # noqa: E402
import reference  # noqa: E402

# Calls before each primer's first C-tract (the window both kinds of read are compared on).
WINDOW = {"HV1F": (20, 170), "HV1R": (20, 195), "HV2F": (20, 250), "HV3R": (50, 275)}


def read_stats(ev: dict, primer: str, template: str | None = None) -> dict:
    lo, hi = WINDOW[primer]
    hi = min(hi, ev["call_count"])
    second = np.array(ev["second"][lo:hi])
    profiles = np.array(ev["profiles"][lo:hi], dtype=float)
    primary = profiles.max(axis=1) if len(profiles) else np.array([np.nan])
    in_phase = np.zeros(ev["call_count"], dtype=bool)
    for segment in ev["segments"]:
        if segment["state"] == "in_phase":
            in_phase[segment["calls"]["start"] : segment["calls"]["end"]] = True
    out = {
        "second_q50": float(np.percentile(second, 50)),
        "second_q90": float(np.percentile(second, 90)),
        "second_q99": float(np.percentile(second, 99)),
        "primary_share_median": float(np.median(primary)),
        "in_phase_fraction": float(in_phase[lo:hi].mean()),
        "withheld_fraction": float(np.mean([bool(w) for w in ev["withheld"][lo:hi]])),
        "stands_out_fraction": float(np.mean(ev["stands_out"][lo:hi])),
    }
    if template is not None:
        out["mismatches"] = sum(
            1 for i in range(lo, hi) if ev["bases"][i] != template[i]
        )
    return out


FIT_KEYS = {
    "spacing@100": lambda r: r.get("spacing_at", {}).get("100"),
    "sigma@100": lambda r: r.get("width_at", {}).get("100"),
    "dye spread": lambda r: max(r["dye_offset"]) - min(r["dye_offset"]),
    "position residual SD": lambda r: r.get("position_jitter_sd"),
    "position residual lag-1": lambda r: r.get("position_resid_lag1"),
    "log-height residual SD": lambda r: r.get("height_resid_sd"),
    "log-height residual lag-1": lambda r: r.get("height_resid_lag1"),
    "floor SD / height": lambda r: float(np.median(r["floor_sd"])) / r["height_median"],
    "floor lag-1": lambda r: r.get("floor_lag1"),
    "zero fraction": lambda r: r.get("floor_zero_fraction"),
    "clean second/primary q50": lambda r: r["second_ratio_q"][1],
    "clean second/primary q90": lambda r: r["second_ratio_q"][2],
    "envelope peak height": lambda r: r.get("envelope_peak_height"),
}


def simulate_one(job: tuple[dict, dict, dict, str, int]) -> dict:
    trace, pooled, calibration, out_dir, seed = job
    placement = trace["placement"]
    params = model.from_fit(trace, pooled, seed=seed, calibration=calibration)
    calls = trace["calls"]
    template = reference.template(placement["strand"], placement["start"], calls)
    rendered = model.render([(template, 1.0)], params, calls=len(template))
    path = Path(out_dir) / f"sim_{seed:05d}_{trace['primer']}.ab1"
    abif.write(path, rendered.channels, rendered.ploc, rendered.pbas, rendered.pcon)
    refit = fit.fit_trace(str(path))
    refit.pop("context", None)
    refit.pop("position_context", None)
    refit.pop("height_resid_calls", None)
    refit.pop("position_resid_calls", None)
    return {
        "path": str(path),
        "template": template,
        "refit": refit,
        "primer": trace["primer"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pooled", type=Path)
    parser.add_argument("fit", type=Path, help="fit.py JSONL of the same traces")
    parser.add_argument("out", type=Path)
    parser.add_argument("--per-primer", type=int, default=50)
    parser.add_argument("--calibration", type=str, default="{}")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--workers", type=int, default=32)
    args = parser.parse_args()
    pooled = model.load_pooled(args.pooled)
    calibration = json.loads(args.calibration)
    fits = {}
    for line in args.fit.open():
        record = json.loads(line)
        fits[record["path"]] = record
    rng = random.Random(args.seed)
    chosen = []
    for primer in WINDOW:
        pool = [
            t
            for t in pooled["traces"]
            if t["primer"] == primer
            and t.get("placement")
            and t["calls"] >= WINDOW[primer][1] + 50
        ]
        chosen += rng.sample(pool, min(args.per_primer, len(pool)))
    sim_dir = args.out / "sim"
    sim_dir.mkdir(parents=True, exist_ok=True)
    jobs = [(t, pooled, calibration, str(sim_dir), k) for k, t in enumerate(chosen)]
    with ProcessPoolExecutor(args.workers) as pool:
        sims = list(pool.map(simulate_one, jobs))
    run_root = args.out / "run"
    real_ev = dnarun.evidence_many(
        [Path(t["path"]) for t in chosen], run_root, args.workers
    )
    sim_ev = dnarun.evidence_many(
        [Path(s["path"]) for s in sims], run_root, args.workers
    )
    rows: dict[str, dict[str, list[float]]] = {}

    def put(primer: str, key: str, kind: str, value: float | None) -> None:
        if value is None or not math.isfinite(value):
            return
        for group in ("all", primer):
            rows.setdefault(key, {}).setdefault(f"{group}:{kind}", []).append(value)

    for trace, sim, rev, sev in zip(chosen, sims, real_ev, sim_ev, strict=True):
        primer = trace["primer"]
        real_fit = fits.get(trace["path"], {})
        for key, fn in FIT_KEYS.items():
            for kind, record in (("real", real_fit), ("sim", sim["refit"])):
                try:
                    put(primer, key, kind, fn(record))
                except (KeyError, TypeError, ValueError, ZeroDivisionError):
                    pass
        if "error" not in rev:
            for key, value in read_stats(rev, primer).items():
                put(primer, "dna " + key, "real", value)
        if "error" not in sev:
            for key, value in read_stats(sev, primer, sim["template"]).items():
                put(primer, "dna " + key, "sim", value)
    groups = ["all", *WINDOW]
    lines = [
        f"# Realism: simulated vs real ({len(chosen)} reads, calibration {calibration})",
        "",
        "Cells: real median / simulated median (KS distance). Window: each primer's calls "
        "before its first C-tract.",
        "",
        "| statistic | " + " | ".join(groups) + " |",
        "|---|" + "---|" * len(groups),
    ]
    for key in rows:
        cells = []
        for group in groups:
            real = rows[key].get(f"{group}:real", [])
            sim = rows[key].get(f"{group}:sim", [])
            if not sim:
                cells.append("-")
                continue
            real_m = f"{np.median(real):.3g}" if real else "-"
            ks = (
                f" ({ks_2samp(real, sim).statistic:.2f})"
                if real and len(sim) > 3
                else ""
            )
            cells.append(f"{real_m} / {np.median(sim):.3g}{ks}")
        lines.append(f"| {key} | " + " | ".join(cells) + " |")
    mismatches = rows.get("dna mismatches", {}).get("all:sim", [])
    lines += [
        "",
        f"Simulated reads whose `dna evidence` calls all equal the template over the window: "
        f"{sum(1 for m in mismatches if m == 0)} of {len(mismatches)}; total differing calls "
        f"{int(sum(mismatches))}.",
        f"`dna evidence` failures: real {sum('error' in e for e in real_ev)}, simulated "
        f"{sum('error' in e for e in sim_ev)}.",
    ]
    (args.out / "realism.md").write_text("\n".join(lines) + "\n")
    (args.out / "realism.json").write_text(json.dumps(rows))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
