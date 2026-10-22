# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy", "scipy"]
# ///
"""Checks the Rust port of the simulator against this research simulator (research only).

PROP-0015 increment 1 accepts the Rust renderer (`dna-trace`, feature `simulate`)
and its bench (`dna-eval sim-bench`) when they reproduce this simulator:

- realism: on pure reads of each bench tract (noise x1), the fitted sigma and
  spacing at call 100, the clean-call second/primary quantiles, and `dna`'s
  second-share quantiles over the pre-tract window, within 5 % of this
  simulator's;
- decoder: the bench measures per tract and noise within 3 points.

Usage:
    uv run --script port_check.py cases OUT.jsonl --reps N
        (the sim-bench case grid, for `bench.py --cases`)
    uv run --script port_check.py realism POOLED.json RUST_TRACES OUT_DIR [--reads N]
        (RUST_TRACES: `dna-eval sim-bench --keep-traces` output's `traces`)
    uv run --script port_check.py compare PY_BENCH_DIR RUST_BENCH.tsv
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
import abif  # noqa: E402
import bench  # noqa: E402
import dnarun  # noqa: E402
import fit  # noqa: E402
import model  # noqa: E402
import realism  # noqa: E402
import reference  # noqa: E402

# The sim-bench grid (crates/dna-eval-lib/src/simbench.rs, case set version 1).
TRACTS = {"HV2-F": 9, "HV2-R": 9, "HV1-F": 11, "HV1-R": 11}
PAIR_VALUES = (0.5, 0.8, 0.9, 1.0, 1.25)
NOISE = (1.0, 2.0, 4.0)
FROZEN_REPS, PURE_REPS = 6, 3
EXPECTED_MARK = {"HV2-F": "T", "HV2-R": "A", "HV1-F": None, "HV1-R": None}
KEY = {"HV2-F": "total", "HV2-R": "total", "HV1-F": "first", "HV1-R": "first"}


def grid(reps: int) -> list[dict]:
    cases = []
    for tract, dominant in TRACTS.items():
        for noise in NOISE:
            for value in PAIR_VALUES:
                for rep in range(reps):
                    cases.append(dict(tract=tract, family="pair", dominant=dominant,
                                      value=value, noise=noise, rep=rep))
            for rep in range(PURE_REPS * reps // FROZEN_REPS):
                cases.append(dict(tract=tract, family="pair", dominant=dominant,
                                  value=0.0, noise=noise, rep=rep))
    return cases


def measures(rows: list[dict]) -> dict[tuple[str, float], dict[str, float]]:
    """The sim-bench measures (simbench.rs `measures`) over `bench.py` rows."""
    out = {}
    for tract in TRACTS:
        for noise in NOISE:
            group = [r for r in rows if r["tract"] == tract and r["noise"] == noise
                     and "error" not in r]
            pure = [r for r in group if r["value"] == 0]
            pair = [r for r in group if r["value"] > 0]
            if not pair:
                continue
            structure, right, right_s, clear_wrong, ratio, tie = [], [], [], [], [], []
            for r in pair:
                if not r.get("decoded"):
                    continue
                s = score(r, KEY[tract])
                if s is None:
                    continue
                ok = r.get("mark") == EXPECTED_MARK[tract]
                structure.append(ok)
                if r["value"] == 1.0:
                    tie.append(s["shorter"])
                    continue
                right.append(s["right"])
                if ok:
                    right_s.append(s["right"])
                clear_wrong.append(s["clear"] and not s["right"])
                if 0.5 <= r["value"] <= 0.9 and s.get("rho_hat") is not None:
                    ratio.append(s["rho_hat"] / r["value"])
            pct = lambda v: 100 * float(np.mean(v)) if v else math.nan  # noqa: E731
            out[(tract, noise)] = {
                "reads": len(pair),
                "decoded": pct([bool(r.get("decoded")) for r in pair]),
                "structure_right": pct(structure),
                "dominant_right": pct(right),
                "right_given_structure": pct(right_s),
                "clear_wrong": pct(clear_wrong),
                "rho_ratio": float(np.median(ratio)) if ratio else math.nan,
                "tie_shorter": pct(tie),
                "pure_decoded": pct([bool(r.get("decoded")) for r in pure]),
            }
    return out


def score(row: dict, key: str) -> dict | None:
    truth = {int(k): float(v) for k, v in row[f"truth_{key}"].items()}
    decoded = row.get(key) or {}
    if decoded.get("dominant") is None:
        return None
    top = max(truth.values())
    dominants = sorted(k for k, v in truth.items() if abs(v - top) < 1e-9)
    shares = {int(k): float(v) for k, v in decoded["shares"].items()}
    big = dominants[0]
    rho_hat = None
    if (big - 1) in truth and shares.get(big, 0) > 0:
        rho_hat = shares.get(big - 1, 0.0) / shares[big]
    return {
        "right": decoded["dominant"] in dominants,
        "shorter": decoded["dominant"] == dominants[0],
        "clear": decoded.get("margin") == "clear",
        "rho_hat": rho_hat,
    }


def compare(py_dir: Path, rust_tsv: Path) -> None:
    rows = [json.loads(line) for line in (py_dir / "cases.jsonl").open()]
    py = measures(rows)
    lines = rust_tsv.read_text().splitlines()
    header = lines[0].split("\t")
    rust = {}
    for line in lines[1:]:
        cells = dict(zip(header, line.split("\t"), strict=True))
        rust[(cells["tract"], float(cells["noise"]))] = {
            k: (float(v) if v != "-" else math.nan) for k, v in cells.items()
            if k not in ("tract", "marginal", "noise")
        }
    names = ["decoded", "structure_right", "dominant_right", "right_given_structure",
             "clear_wrong", "tie_shorter", "pure_decoded", "rho_ratio"]
    print("| tract | noise | " + " | ".join(names) + " |")
    print("|---|---:|" + "---|" * len(names))
    worst: dict[str, float] = {}
    for key in sorted(py):
        cells = []
        for name in names:
            a, b = py[key][name], rust.get(key, {}).get(name, math.nan)
            if name == "rho_ratio":
                cells.append(f"{a:.3f} / {b:.3f}")
                d = abs(a - b) * 100
            else:
                cells.append(f"{a:.0f} / {b:.0f}")
                d = abs(a - b)
            if math.isfinite(d):
                worst[name] = max(worst.get(name, 0.0), d)
        print(f"| {key[0]} | {key[1]:g} | " + " | ".join(cells) + " |")
    print("\nResearch / Rust. Largest absolute difference (points; rho_ratio x100):",
          ", ".join(f"{k} {v:.1f}" for k, v in worst.items()))
    # Pooled over noise per tract, the measures the acceptance reads.
    print("\nPooled per tract (reads-weighted means of the cells):")
    for tract in TRACTS:
        parts = []
        for name in names:
            vals = [(py[(tract, n)][name], rust[(tract, n)][name], py[(tract, n)]["reads"])
                    for n in NOISE if (tract, n) in py and (tract, n) in rust]
            vals = [(a, b, w) for a, b, w in vals if math.isfinite(a) and math.isfinite(b)]
            if not vals:
                continue
            w = sum(v[2] for v in vals)
            a = sum(v[0] * v[2] for v in vals) / w
            b = sum(v[1] * v[2] for v in vals) / w
            scale = 100 if name == "rho_ratio" else 1
            parts.append(f"{name} {a * scale:.1f}/{b * scale:.1f}")
        print(f"- {tract}: " + "; ".join(parts))


def render_python(job: tuple[dict, dict, str]) -> str:
    """One pure research-simulator read of a bench case, written as ABIF."""
    case, pooled, out_dir = job
    tract = bench.TRACTS[case["tract"]]
    rng = random.Random(bench.case_id(case))
    pool = [t for t in pooled["traces"] if t["primer"] == tract.primer and t["span"][1] >= 150]
    trace = rng.choice(pool)
    params = model.from_fit(trace, pooled, seed=rng.randrange(2**31)).scaled_noise(1.0)
    length = case["dominant"]
    probe = reference.template_map(tract.strand, tract.start, 2000, tract.edits(length))
    if tract.strand == "F":
        call = next(i for i, p in enumerate(probe) if p >= tract.positions[0])
    else:
        call = next(i for i, p in enumerate(probe) if p <= tract.positions[1])
    calls = call + 20 + tract.after
    template = reference.template(tract.strand, tract.start, calls + 5, tract.edits(length))
    rendered = model.render([(template, 1.0)], params, calls=calls)
    path = Path(out_dir) / f"py_{bench.case_id(case)}.ab1"
    abif.write(path, rendered.channels, rendered.ploc, rendered.pbas, rendered.pcon)
    return str(path)


def measure(path: str) -> dict:
    record = fit.fit_trace(path)
    return {
        "sigma@100": record.get("width_at", {}).get("100"),
        "spacing@100": record.get("spacing_at", {}).get("100"),
        "clean second q50": (record.get("second_ratio_q") or [None] * 3)[1],
        "clean second q90": (record.get("second_ratio_q") or [None] * 3)[2],
    }


def realism_check(pooled_path: Path, rust_dir: Path, out: Path, reads: int) -> None:
    pooled = model.load_pooled(pooled_path)
    out.mkdir(parents=True, exist_ok=True)
    primers = {"HV2-F": "HV2F", "HV2-R": "HV3R", "HV1-F": "HV1F", "HV1-R": "HV1R"}
    cases = [dict(tract=t, family="pair", dominant=d, value=0.0, noise=1.0, rep=r)
             for t, d in TRACTS.items() for r in range(reads)]
    with ProcessPoolExecutor(48) as pool:
        py_paths = list(pool.map(render_python, [(c, pooled, str(out)) for c in cases]))
    rust_paths = []
    for c in cases:
        name = f"{c['tract']}_pair_{c['dominant']}_0_n1_r{c['rep']}.ab1"
        path = rust_dir / name
        rust_paths.append(str(path) if path.exists() else None)
    kinds = {"research": py_paths, "rust": rust_paths}
    stats: dict[str, dict[str, dict[str, list[float]]]] = {}
    for kind, paths in kinds.items():
        valid = [p for p in paths if p]
        with ProcessPoolExecutor(48) as pool:
            fitted = list(pool.map(measure, valid))
        evidence = dnarun.evidence_many([Path(p) for p in valid], out / "run", 48)
        for path, f, ev in zip(valid, fitted, evidence, strict=True):
            tract = Path(path).name.removeprefix("py_").split("_")[0]
            primer = primers[tract]
            values = dict(f)
            if "error" not in ev:
                s = realism.read_stats(ev, primer)
                values.update({f"dna second {q}": s[f"second_{q}"] for q in ("q50", "q90", "q99")})
            for key, value in values.items():
                if value is None or not math.isfinite(value):
                    continue
                for group in ("all", primer):
                    stats.setdefault(group, {}).setdefault(key, {}).setdefault(kind, []).append(value)
    print("| group | statistic | research median | Rust median | difference |")
    print("|---|---|---:|---:|---:|")
    for group in sorted(stats):
        for key, kinds_ in stats[group].items():
            a = float(np.median(kinds_.get("research", [math.nan])))
            b = float(np.median(kinds_.get("rust", [math.nan])))
            print(f"| {group} | {key} | {a:.4g} | {b:.4g} | {100 * (b - a) / a:+.1f} % |")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    c = sub.add_parser("cases")
    c.add_argument("out", type=Path)
    c.add_argument("--reps", type=int, default=FROZEN_REPS)
    r = sub.add_parser("realism")
    r.add_argument("pooled", type=Path)
    r.add_argument("rust", type=Path)
    r.add_argument("out", type=Path)
    r.add_argument("--reads", type=int, default=30)
    m = sub.add_parser("compare")
    m.add_argument("python", type=Path)
    m.add_argument("rust", type=Path)
    args = parser.parse_args()
    if args.command == "cases":
        with args.out.open("w") as handle:
            for case in grid(args.reps):
                handle.write(json.dumps(case) + "\n")
    elif args.command == "realism":
        realism_check(args.pooled, args.rust, args.out, args.reads)
    else:
        compare(args.python, args.rust)


if __name__ == "__main__":
    main()
