# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy", "scipy"]
# ///
"""Replays real tract decodes through the simulator (research only).

For each real read whose C-tract `dna evidence` decodes, the decoded length
distribution becomes a simulated read's truth, rendered with that read's own
fitted parameters. Decoding the replay shows (a) whether the simulator's dephased
signal fits the decoder as well as the real one does (residual, autocorrelation,
margins), and (b) how well the decoder recovers a distribution shaped like real
ones.

Usage:
    uv run --script replay.py cases POOLED.json OUT_DIR [--per-primer N]
    uv run --script bench.py POOLED.json OUT_DIR --cases OUT_DIR/replay-cases.jsonl
    uv run --script replay.py compare OUT_DIR
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bench  # noqa: E402
import dnarun  # noqa: E402
import model  # noqa: E402
import reference  # noqa: E402

TRACT_OF = {"HV2F": "HV2-F", "HV3R": "HV2-R", "HV1F": "HV1-F", "HV1R": "HV1-R"}
MIN_SHARE = 0.03


def expected_call(tract: bench.Tract, start: int) -> int:
    probe = reference.template_map(
        tract.strand, start, 2000, tract.edits(9 if tract.fixed else 11)
    )
    if tract.strand == "F":
        return next(i for i, p in enumerate(probe) if p >= tract.positions[0])
    return next(i for i, p in enumerate(probe) if p <= tract.positions[1])


def shares(marginal: dict) -> dict[int, float]:
    kept = {
        int(k): v
        for k, v in zip(marginal["lengths"], marginal["shares"], strict=True)
        if v >= MIN_SHARE
    }
    total = sum(kept.values())
    return {k: v / total for k, v in kept.items()} if total else {}


def build_cases(pooled_path: Path, out: Path, per_primer: int, workers: int) -> None:
    pooled = model.load_pooled(pooled_path)
    rng = random.Random(11)
    chosen = []
    for primer in TRACT_OF:
        pool = [
            t for t in pooled["traces"] if t["primer"] == primer and t.get("placement")
        ]
        chosen += rng.sample(pool, min(per_primer, len(pool)))
    evs = dnarun.evidence_many(
        [Path(t["path"]) for t in chosen], out / "run-real", workers
    )
    cases, real = [], []
    for trace, ev in zip(chosen, evs, strict=True):
        if "error" in ev:
            continue
        name = TRACT_OF[trace["primer"]]
        tract = bench.TRACTS[name]
        start = trace["placement"]["start"]
        call = expected_call(tract, start)
        matched = [t for t in ev["tracts"] if abs(t["calls"]["start"] - call) <= 8]
        if not matched:
            continue
        decode = matched[0]
        first, total = decode["first"], decode["total"]
        if tract.fixed:
            if decode.get("mark") != ("T" if tract.strand == "F" else "A"):
                continue
            if tract.varying == "first":
                fixed = total["dominant"] - first["dominant"] - 1
                fractions = shares(first)
            else:
                fixed = first["dominant"]
                fractions = {k - fixed - 1: v for k, v in shares(total).items()}
            if not 4 <= fixed <= 8:
                continue
        else:
            if decode.get("mark") is not None:
                continue  # not a T16189C read: no single run
            fixed = 0
            fractions = shares(first)
        lo = 6 if name.startswith("HV2") else 9
        if not fractions or min(fractions) < lo:
            continue
        case = {
            "tract": name,
            "family": "replay",
            "dominant": max(fractions, key=fractions.get),
            "value": 0.0,
            "noise": 1.0,
            "rep": len(cases),
            "fixed": fixed,
            "fractions": fractions,
            "source": trace["path"],
        }
        cases.append(case)
        real.append(
            {
                "case": bench.case_id(case),
                "residual": decode["residual"],
                "autocorrelation": decode["autocorrelation"],
                "first_z": first.get("z_adjusted"),
                "total_z": total.get("z_adjusted"),
                "first_margin": first.get("margin"),
                "total_margin": total.get("margin"),
                "start": start,
                "call": call,
                "decode_start": decode["calls"]["start"],
            }
        )
    with (out / "replay-cases.jsonl").open("w") as handle:
        for case in cases:
            handle.write(json.dumps(case) + "\n")
    (out / "replay-real.json").write_text(json.dumps(real))
    print(f"{len(cases)} replay cases from {len(chosen)} real reads")


def tv(a: dict, b: dict) -> float:
    keys = set(a) | set(b)
    return 0.5 * sum(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in keys)


def compare(out: Path) -> None:
    real = {r["case"]: r for r in json.loads((out / "replay-real.json").read_text())}
    rows = [json.loads(line) for line in (out / "cases.jsonl").open()]
    lines = [
        "| tract | reads | decoded | residual real / sim | autocorr real / sim | "
        "varying-marginal z real / sim | sim dominant right | sim TV to truth (median) |",
        "|---|---:|---:|---|---|---|---:|---:|",
    ]
    for name in ("HV2-F", "HV2-R", "HV1-F", "HV1-R"):
        group = [r for r in rows if r["tract"] == name]
        dec = [r for r in group if r.get("decoded")]
        if not group:
            continue
        key = "total" if bench.TRACTS[name].varying == "total" else "first"
        rr = [real[bench.case_id(r)] for r in dec]
        res = (
            np.median([x["residual"] for x in rr]),
            np.median([r["residual"] for r in dec]),
        )
        ac = (
            np.median([x["autocorrelation"] for x in rr]),
            np.median([r["autocorrelation"] for r in dec]),
        )
        zr = [x[f"{key}_z"] for x in rr if x[f"{key}_z"] is not None]
        zs = [r[key]["z"] for r in dec if r[key]["z"] is not None]
        right = 0
        tvs = []
        for r in dec:
            truth = {int(k): v for k, v in r[f"truth_{key}"].items()}
            dom = max(truth, key=truth.get)
            right += r[key]["dominant"] == dom
            tvs.append(tv({int(k): v for k, v in r[key]["shares"].items()}, truth))
        lines.append(
            f"| {name} | {len(group)} | {len(dec)} | {res[0]:.4f} / {res[1]:.4f} | "
            f"{ac[0]:.2f} / {ac[1]:.2f} | {np.median(zr):.1f} / {np.median(zs):.1f} | "
            f"{right}/{len(dec)} | {np.median(tvs):.3f} |"
        )
    text = "\n".join(lines) + "\n"
    (out / "replay.md").write_text(text)
    print(text)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("cases", "compare"))
    parser.add_argument("args", nargs="+")
    parser.add_argument("--per-primer", type=int, default=150)
    parser.add_argument("--workers", type=int, default=48)
    args = parser.parse_args()
    if args.command == "cases":
        pooled, out = Path(args.args[0]), Path(args.args[1])
        out.mkdir(parents=True, exist_ok=True)
        build_cases(pooled, out.resolve(), args.per_primer, args.workers)
    else:
        compare(Path(args.args[0]))


if __name__ == "__main__":
    main()
