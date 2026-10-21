# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy", "scipy"]
# ///
"""Renders a synthetic Sanger trace (ABIF) with known truth from a JSON spec (research only).

Spec (all keys but `read` optional):

    {
      "read": {"strand": "F", "start": 36, "calls": 600},     # rCRS read, or
      "template": "ACGT...",                                   # an explicit sequence
      "molecules": [{"edits": {"309": "CC"}, "fraction": 0.6},
                    {"edits": {"309": "CCC"}, "fraction": 0.4}],
      "params": {"primer": "HV2F", "seed": 1, "noise": 1.0},  # a fitted read of that primer,
                                                               # or {"trace": "<path>"} for one
      "artefacts": [{"kind": "spike", "call": 120, "height": 3},
                    {"kind": "blob", "call": 70, "channel": "G", "height": 1.5, "width": 4},
                    {"kind": "noise", "call": 300, "width": 20, "height": 0.2}]
    }

`edits` map a 1-based rCRS position to the reference-strand text replacing its base
(see `reference.template`); a point heteroplasmy is a molecule with a base edit and
its fraction. Writes OUT.ab1 (DATA.9-12, FWO_, PLOC, PBAS, PCON) and OUT.truth.json
(the molecules' read-order sequences and fractions, the peak grid, the parameters).

Usage:
    uv run --script simulate.py POOLED.json SPEC.json OUT_STEM
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import abif  # noqa: E402
import model  # noqa: E402
import reference  # noqa: E402


def build(pooled: dict, spec: dict) -> tuple[model.Rendered, dict]:
    settings = spec.get("params", {})
    seed = int(settings.get("seed", 0))
    rng = random.Random(seed)
    if "trace" in settings:
        trace = next(t for t in pooled["traces"] if t["path"] == settings["trace"])
    else:
        primer = settings.get("primer", "HV2F")
        trace = rng.choice([t for t in pooled["traces"] if t["primer"] == primer])
    params = model.from_fit(trace, pooled, seed=seed).scaled_noise(
        float(settings.get("noise", 1.0))
    )
    molecules = spec.get("molecules") or [{"edits": {}, "fraction": 1.0}]
    if "template" in spec:
        templates = [(spec["template"], 1.0)]
        calls = len(spec["template"])
    else:
        read = spec["read"]
        calls = int(read["calls"])
        templates = [
            (
                reference.template(
                    read["strand"],
                    int(read["start"]),
                    calls,
                    {int(k): v for k, v in m.get("edits", {}).items()},
                ),
                float(m["fraction"]),
            )
            for m in molecules
        ]
    artefacts = [
        model.Artefact(
            kind=a["kind"],
            call=float(a.get("call", 0)),
            channel="ACGT".index(a.get("channel", "A")),
            height=float(a.get("height", 1.0)),
            width=float(a.get("width", 3.0)),
        )
        for a in spec.get("artefacts", [])
    ]
    rendered = model.render(templates, params, calls=calls, artefacts=artefacts)
    truth = {
        "spec": spec,
        "source_trace": trace["path"],
        "molecules": [{"sequence": t, "fraction": f} for t, f in templates],
        "grid": [round(float(g), 3) for g in rendered.truth_grid],
        "params": {
            k: v
            for k, v in dataclasses.asdict(params).items()
            if k not in ("height_context", "position_context")
        },
    }
    return rendered, truth


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pooled", type=Path)
    parser.add_argument("spec", type=Path)
    parser.add_argument(
        "out", type=Path, help="output stem (writes .ab1 and .truth.json)"
    )
    args = parser.parse_args()
    pooled = model.load_pooled(args.pooled)
    rendered, truth = build(pooled, json.loads(args.spec.read_text()))
    abif.write(
        args.out.with_suffix(".ab1"),
        rendered.channels,
        rendered.ploc,
        rendered.pbas,
        rendered.pcon,
    )
    args.out.with_suffix(".truth.json").write_text(json.dumps(truth))
    print(f"wrote {args.out.with_suffix('.ab1')} ({len(rendered.ploc)} calls)")


if __name__ == "__main__":
    main()
