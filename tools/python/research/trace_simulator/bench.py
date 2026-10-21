# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy", "scipy"]
# ///
"""Benchmarks DNA's tract decode on simulated C-tract length mixtures (research only).

Each case renders a read through one C-tract as a mixture of molecules whose run
lengths follow a set distribution (the truth), with a real read's fitted
parameters, writes it as ABIF, runs `dna evidence`, and scores the decode of the
tract against the truth: its first-run and total marginals' dominant lengths, the
share vector, the ratio of the shorter neighbour to the dominant, and the margin
state. The comparison lives here, in research, never in DNA (ADR-0065, ADR-0075).

Usage:
    uv run --script bench.py POOLED.json OUT_DIR [--grid main|ablation|smoke] [--reps N]
    uv run --script bench.py POOLED.json OUT_DIR --cases CASES.jsonl
    uv run --script report.py OUT_DIR   (the tables)
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import abif  # noqa: E402
import dnarun  # noqa: E402
import model  # noqa: E402
import reference  # noqa: E402


@dataclass(frozen=True)
class Tract:
    """A C-tract read by one primer: where, how its molecules are written, its truth."""

    name: str
    primer: str
    strand: str
    start: int  # rCRS position of call 0
    positions: tuple[int, int]  # the tract's rCRS span
    after: int  # calls the read holds after the tract
    varying: str  # "first" or "total": which marginal carries the varying run
    fixed: int  # the other run's length (0: no mark)

    def edits(self, length: int) -> dict[int, str]:
        """rCRS edits giving the varying run `length` C (on the reference strand)."""
        if self.name.startswith("HV2"):
            # Run 1 is 303-309 (C7), T310, run 2 is 311-315 (C5); `fixed` sets run 2.
            return {309: "C" * (length - 6), 315: "C" * (self.fixed - 4)}
        if self.name.startswith("HV1"):
            # T16189C joins 16184-16193 into one run of ten.
            return {16189: "C", 16193: "C" * (length - 9)}
        return {573: "C" * (length - 5)}  # HV3: 568-573 (C6)

    def marginals(self, length: int) -> tuple[int, int]:
        """The (first, total) lengths a molecule with varying run `length` shows."""
        if not self.fixed:
            return length, length
        if self.varying == "first":
            return length, length + 1 + self.fixed
        return self.fixed, self.fixed + 1 + length


TRACTS = {
    "HV2-F": Tract("HV2-F", "HV2F", "F", 36, (303, 315), 200, "first", 6),
    "HV2-R": Tract("HV2-R", "HV3R", "R", 604, (303, 316), 200, "total", 6),
    "HV1-F": Tract("HV1-F", "HV1F", "F", 16002, (16184, 16193), 200, "first", 0),
    "HV1-R": Tract("HV1-R", "HV1R", "R", 16401, (16184, 16193), 200, "first", 0),
    "HV3-F": Tract("HV3-F", "HV2F", "F", 36, (568, 573), 120, "first", 0),
    "HV3-R": Tract("HV3-R", "HV3R", "R", 604, (568, 573), 200, "first", 0),
}


def mixture(family: str, dominant: int, value: float) -> dict[int, float]:
    """Molecule fractions over the varying run's lengths."""
    if family == "pair":  # shorter neighbour at `value` x the dominant
        f = {dominant: 1.0, dominant - 1: value}
    elif family == "spread":  # geometric on both sides, ratio `value` per step
        f = {dominant + d: value ** abs(d) for d in range(-2, 3)}
    elif family == "stutter":  # n-1 at `value`, n-2 at value^2/2, n+1 at 0.1
        f = {
            dominant: 1.0,
            dominant - 1: value,
            dominant - 2: value * value / 2,
            dominant + 1: 0.1,
        }
    else:
        raise ValueError(family)
    total = sum(f.values())
    return {k: v / total for k, v in f.items() if v > 0}


def case_id(case: dict) -> str:
    ablate = "-".join(case.get("ablate", []))
    return "{tract}_{family}_{dominant}_{value:g}_n{noise:g}_r{rep}".format(**case) + (
        f"_x{ablate}" if ablate else ""
    )


def run_case(job: tuple[dict, dict, str, str]) -> dict:
    case, pooled, out_dir, run_root = job
    tract = TRACTS[case["tract"]]
    if "fixed" in case:  # the other run's length, when the case sets it
        tract = replace(tract, fixed=case["fixed"])
    rng = random.Random(case_id(case))
    if case.get("source"):  # replay: the fitted parameters of one named real read
        trace = next(t for t in pooled["traces"] if t["path"] == case["source"])
    else:
        pool = [
            t
            for t in pooled["traces"]
            if t["primer"] == tract.primer and t["span"][1] >= 150
        ]
        trace = rng.choice(pool)
    params = model.from_fit(trace, pooled, seed=rng.randrange(2**31))
    params = params.scaled_noise(case["noise"]).ablated(case.get("ablate", []))
    fractions = case.get("fractions") or mixture(
        case["family"], case["dominant"], case["value"]
    )
    fractions = {int(k): v for k, v in fractions.items()}
    reference_length = max(fractions)
    probe_map = reference.template_map(
        tract.strand, tract.start, 2000, tract.edits(reference_length)
    )
    if tract.strand == "F":
        tract_call = next(i for i, p in enumerate(probe_map) if p >= tract.positions[0])
    else:
        tract_call = next(i for i, p in enumerate(probe_map) if p <= tract.positions[1])
    calls = tract_call + 20 + tract.after
    templates = [
        (
            reference.template(
                tract.strand, tract.start, calls + 5, tract.edits(length)
            ),
            share,
        )
        for length, share in fractions.items()
    ]
    rendered = model.render(templates, params, calls=calls)
    path = Path(out_dir) / f"{case_id(case)}.ab1"
    abif.write(path, rendered.channels, rendered.ploc, rendered.pbas, rendered.pcon)
    ev = dnarun.evidence(path, Path(run_root))
    path.unlink()
    row = {
        **case,
        "fractions": fractions,
        "tract_call": tract_call,
        "source": trace["path"],
    }
    if "error" in ev:
        row["error"] = ev["error"]
        return row
    first_truth: dict[int, float] = defaultdict(float)
    total_truth: dict[int, float] = defaultdict(float)
    for length, share in fractions.items():
        a, t = tract.marginals(length)
        first_truth[a] += share
        total_truth[t] += share
    row["truth_first"] = dict(first_truth)
    row["truth_total"] = dict(total_truth)
    # The run's length in the read's own primaries (what an in-phase read reports).
    bases = ev["bases"]
    run_base = bases[tract_call] if tract_call < len(bases) else "N"
    k = tract_call
    while k < len(bases) and bases[k] == run_base:
        k += 1
    row["primary_run"] = k - tract_call
    in_phase = all(
        s["state"] == "in_phase"
        for s in ev["segments"]
        if s["calls"]["start"] < tract_call + 40 and s["calls"]["end"] > tract_call - 5
    )
    row["in_phase"] = in_phase
    matched = [t for t in ev["tracts"] if abs(t["calls"]["start"] - tract_call) <= 3]
    row["decoded"] = bool(matched)
    if matched:
        t = matched[0]
        row["decode_start"] = t["calls"]["start"]
        row["mark"] = t.get("mark")
        row["residual"] = t.get("residual")
        row["autocorrelation"] = t.get("autocorrelation")
        for key in ("first", "total"):
            m = t[key]
            row[key] = {
                "shares": dict(zip(m["lengths"], m["shares"], strict=True)),
                "dominant": m.get("dominant"),
                "next": m.get("next"),
                "z": m.get("z_adjusted"),
                "margin": m.get("margin"),
            }
        cleaned = t.get("cleaned") or {}
        if cleaned.get("first"):
            m = cleaned["first"]
            row["cleaned_first"] = {
                "shares": dict(zip(m["lengths"], m["shares"], strict=True)),
                "dominant": m.get("dominant"),
                "margin": m.get("margin"),
            }
    row["second_q50"] = float(np.median(ev["second"][20 : max(tract_call - 10, 30)]))
    return row


def grid(name: str, reps: int) -> list[dict]:
    cases = []
    if name == "ablation":
        # Which model part biases the decoded shares near a tie: the same cases with
        # noise, then height variation, then position variation switched off.
        for ablate in (
            [],
            ["noise"],
            ["noise", "heights"],
            ["noise", "heights", "position"],
        ):
            for tract in ("HV2-F", "HV2-R", "HV1-F", "HV1-R"):
                for value in (0.5, 0.8, 0.9, 1.0, 1.1, 1.25):
                    for rep in range(reps):
                        cases.append(
                            {
                                "tract": tract,
                                "family": "pair",
                                "dominant": 9 if "HV2" in tract else 11,
                                "value": value,
                                "noise": 1.0,
                                "rep": rep,
                                "ablate": ablate,
                            }
                        )
        return cases
    if name == "smoke":
        for tract in ("HV2-F", "HV2-R", "HV1-F", "HV1-R"):
            for value in (0.5, 1.0):
                cases.append(
                    {
                        "tract": tract,
                        "family": "pair",
                        "dominant": 9 if "HV2" in tract else 11,
                        "value": value,
                        "noise": 1.0,
                        "rep": 0,
                    }
                )
        return cases
    pair_values = (0.0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.25, 1.5, 2.0)
    for tract in ("HV2-F", "HV2-R", "HV1-F", "HV1-R"):
        dominants = (8, 9) if "HV2" in tract else (11, 12)
        for dominant in dominants:
            for noise in (1.0, 2.0, 4.0):
                for value in pair_values:
                    for rep in range(reps):
                        cases.append(
                            {
                                "tract": tract,
                                "family": "pair",
                                "dominant": dominant,
                                "value": value,
                                "noise": noise,
                                "rep": rep,
                            }
                        )
                for value in (0.1, 0.25, 0.4, 0.6):
                    for rep in range(reps):
                        cases.append(
                            {
                                "tract": tract,
                                "family": "spread",
                                "dominant": dominant,
                                "value": value,
                                "noise": noise,
                                "rep": rep,
                            }
                        )
                for value in (0.3, 0.6, 0.9):
                    for rep in range(reps):
                        cases.append(
                            {
                                "tract": tract,
                                "family": "stutter",
                                "dominant": dominant,
                                "value": value,
                                "noise": noise,
                                "rep": rep,
                            }
                        )
    for tract in ("HV3-F", "HV3-R"):
        for value in pair_values:
            for rep in range(reps):
                cases.append(
                    {
                        "tract": tract,
                        "family": "pair",
                        "dominant": 7,
                        "value": value,
                        "noise": 1.0,
                        "rep": rep,
                    }
                )
    return cases


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pooled", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--grid", default="main")
    parser.add_argument(
        "--cases", type=Path, help="JSONL of explicit cases (with fractions)"
    )
    parser.add_argument("--reps", type=int, default=8)
    parser.add_argument("--workers", type=int, default=48)
    args = parser.parse_args()
    pooled = model.load_pooled(args.pooled)
    args.out.mkdir(parents=True, exist_ok=True)
    sim_dir = args.out / "traces"
    sim_dir.mkdir(exist_ok=True)
    run_root = dnarun.workdir((args.out / "run").resolve())
    if args.cases:
        cases = [json.loads(line) for line in args.cases.open()]
    else:
        cases = grid(args.grid, args.reps)
    jobs = [(case, pooled, str(sim_dir), str(run_root)) for case in cases]
    done = 0
    with (
        ProcessPoolExecutor(args.workers) as pool,
        (args.out / "cases.jsonl").open("w") as handle,
    ):
        for row in pool.map(run_case, jobs, chunksize=2):
            handle.write(json.dumps(row) + "\n")
            done += 1
            if done % 200 == 0:
                print(f"{done}/{len(jobs)}", flush=True)
    print(f"{len(jobs)} cases -> {args.out / 'cases.jsonl'}")


if __name__ == "__main__":
    main()
