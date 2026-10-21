# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy", "scipy"]
# ///
"""What `dna evidence` makes of each simulated artefact (research only).

Renders clean HV2F reads, each with one artefact at a known call (spike, dye blob,
noise burst, a point heteroplasmy as a second molecule), runs `dna evidence`, and
counts, within three calls of the artefact: calls that differ from the dominant
template, calls withheld and why, the second-peak causes, and second peaks that
stand out. A control read without the artefact gives the baseline.

Usage:
    uv run --script artefacts.py POOLED.json OUT_DIR [--reads N]
"""

from __future__ import annotations

import argparse
import random
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import abif  # noqa: E402
import dnarun  # noqa: E402
import model  # noqa: E402
import reference  # noqa: E402

START, CALLS, AT = 36, 300, 120
# rCRS position read at call AT by an HV2F read starting at 36, for the point heteroplasmy.
SNP_POSITION = START + AT
KINDS = {
    "control": [],
    "spike x3": [model.Artefact("spike", call=AT, height=3.0)],
    "spike x6": [model.Artefact("spike", call=AT, height=6.0)],
    "spike x12": [model.Artefact("spike", call=AT, height=12.0)],
    "dye blob G x0.6": [
        model.Artefact("blob", call=AT, channel=2, height=0.6, width=4.0)
    ],
    "dye blob G x0.25": [
        model.Artefact("blob", call=AT, channel=2, height=0.25, width=4.0)
    ],
    "noise burst x0.3": [model.Artefact("noise", call=AT - 5, height=0.3, width=10.0)],
    "point heteroplasmy 30 %": [],
    "point heteroplasmy 15 %": [],
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pooled", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--reads", type=int, default=30)
    parser.add_argument("--workers", type=int, default=32)
    args = parser.parse_args()
    pooled = model.load_pooled(args.pooled)
    args.out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(5)
    traces = rng.sample(
        [t for t in pooled["traces"] if t["primer"] == "HV2F"], args.reads
    )
    clean_read = reference.template("F", START, CALLS)
    other = next(b for b in "ACGT" if b not in clean_read[AT - 1 : AT + 2])
    jobs = []
    for kind, artefacts in KINDS.items():
        for k, trace in enumerate(traces):
            params = model.from_fit(trace, pooled, seed=k)
            clean = reference.template("F", START, CALLS)
            templates = [(clean, 1.0)]
            if kind.startswith("point"):
                share = 0.30 if "30" in kind else 0.15
                variant = reference.template("F", START, CALLS, {SNP_POSITION: other})
                templates = [(clean, 1 - share), (variant, share)]
            rendered = model.render(templates, params, calls=CALLS, artefacts=artefacts)
            path = args.out / f"{kind.replace(' ', '_').replace('%', 'pc')}_{k}.ab1"
            abif.write(
                path, rendered.channels, rendered.ploc, rendered.pbas, rendered.pcon
            )
            jobs.append((kind, path, clean))
    evs = dnarun.evidence_many(
        [p for _, p, _ in jobs], (args.out / "run").resolve(), args.workers
    )
    lines = [
        f"| artefact (at call {AT}) | reads | reads with a miscall | withheld calls (reasons) | "
        "second-peak causes | second peaks standing out |",
        "|---|---:|---:|---|---|---:|",
    ]
    for kind in KINDS:
        rows = [
            (ev, clean)
            for (k, _, clean), ev in zip(jobs, evs, strict=True)
            if k == kind
        ]
        miscall, reasons, causes, stand = 0, Counter(), Counter(), 0
        for ev, clean in rows:
            if "error" in ev:
                continue
            window = range(AT - 3, AT + 4)
            miscall += any(ev["bases"][i] != clean[i] for i in window)
            for i in window:
                reasons.update(ev["withheld"][i])
                if ev["cause"][i]:
                    causes[ev["cause"][i]] += 1
                stand += ev["stands_out"][i]
        fmt = lambda c: ", ".join(f"{k} {v}" for k, v in c.most_common(4)) or "-"  # noqa: E731
        lines.append(
            f"| {kind} | {len(rows)} | {miscall} | {fmt(reasons)} | {fmt(causes)} | {stand} |"
        )
    text = "\n".join(lines) + "\n"
    (args.out / "artefacts.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
