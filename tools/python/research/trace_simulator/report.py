# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy"]
# ///
"""Summarises `bench.py` cases into the decoder benchmark tables (research only).

Usage:
    uv run --script report.py OUT_DIR   (reads OUT_DIR/cases.jsonl, writes OUT_DIR/report.md)
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bench  # noqa: E402


def as_int(d: dict) -> dict[int, float]:
    return {int(k): float(v) for k, v in d.items()}


def tv(a: dict[int, float], b: dict[int, float]) -> float:
    return 0.5 * sum(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in set(a) | set(b))


def varying(row: dict) -> str:
    return bench.TRACTS[row["tract"]].varying


def score(row: dict, key: str) -> dict:
    """One decode's scores on marginal `key` against the truth."""
    truth = as_int(row[f"truth_{key}"])
    top = max(truth.values())
    dominants = sorted(k for k, v in truth.items() if abs(v - top) < 1e-9)
    out = {"decoded": row.get("decoded", False), "tie": len(dominants) > 1}
    if not out["decoded"]:
        # An undecoded read keeps its primaries: the in-phase run length.
        out["primary_right"] = row.get("primary_run") in (
            dominants if key == "first" or not bench.TRACTS[row["tract"]].fixed else []
        )
        return out
    decoded = row[key]
    shares = as_int(decoded["shares"])
    out["dominant"] = decoded["dominant"]
    out["right"] = decoded["dominant"] in dominants
    out["shorter"] = decoded["dominant"] == dominants[0]
    out["margin"] = decoded["margin"]
    out["tv"] = tv(shares, truth)
    big = max(truth, key=lambda k: (truth[k], -k))  # the longer one on a tie
    if (big - 1) in truth and shares.get(big, 0) > 0:
        out["rho_hat"] = shares.get(big - 1, 0.0) / shares[big]
    return out


def pct(values: list[bool]) -> str:
    return f"{100 * np.mean(values):.0f}" if values else "-"


def pair_table(rows: list[dict], tract: str, key: str) -> list[str]:
    lines = [
        f"#### {tract}, {key} marginal",
        "",
        "| noise | rho | n | decoded % | dominant right % | clear % | near-equal % | "
        "clear and wrong % | rho-hat median (q10-q90) | TV median |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---|---:|",
    ]
    by = defaultdict(list)
    for r in rows:
        if r["tract"] == tract and r["family"] == "pair":
            by[(r["noise"], r["value"])].append(score(r, key))
    for (noise, rho), scores in sorted(by.items()):
        dec = [s for s in scores if s["decoded"]]
        right = [s["right"] for s in dec]
        if rho == 1.0:
            right_txt = f"shorter {pct([s['shorter'] for s in dec])}"
        else:
            right_txt = pct(right)
        clear = [s["margin"] == "clear" for s in dec]
        near = [s["margin"] == "near_equal" for s in dec]
        clear_wrong = [s["margin"] == "clear" and not s["right"] for s in dec]
        if rho == 1.0:
            clear_wrong = [False for _ in dec]
        rh = [s["rho_hat"] for s in dec if "rho_hat" in s]
        rh_txt = (
            f"{np.median(rh):.2f} ({np.percentile(rh, 10):.2f}-{np.percentile(rh, 90):.2f})"
            if rh and rho > 0
            else "-"
        )
        tvs = [s["tv"] for s in dec]
        lines.append(
            f"| {noise:g} | {rho:g} | {len(scores)} | {pct([s['decoded'] for s in scores])} | "
            f"{right_txt} | {pct(clear)} | {pct(near)} | {pct(clear_wrong)} | {rh_txt} | "
            f"{np.median(tvs) if tvs else float('nan'):.3f} |"
        )
    lines.append("")
    return lines


def family_table(rows: list[dict], family: str) -> list[str]:
    lines = [
        f"#### {family}",
        "",
        "| tract | marginal | noise | value | n | decoded % | dominant right % | clear % | "
        "clear and wrong % | TV median |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    by = defaultdict(list)
    for r in rows:
        if r["family"] == family:
            by[(r["tract"], r["noise"], r["value"])].append(r)
    for (tract, noise, value), group in sorted(by.items()):
        keys = ["first", "total"] if tract == "HV2-F" else [varying(group[0])]
        for key in keys:
            scores = [score(r, key) for r in group]
            dec = [s for s in scores if s["decoded"]]
            lines.append(
                f"| {tract} | {key} | {noise:g} | {value:g} | {len(scores)} | "
                f"{pct([s['decoded'] for s in scores])} | {pct([s['right'] for s in dec])} | "
                f"{pct([s['margin'] == 'clear' for s in dec])} | "
                f"{pct([s['margin'] == 'clear' and not s['right'] for s in dec])} | "
                f"{np.median([s['tv'] for s in dec]) if dec else float('nan'):.3f} |"
            )
    lines.append("")
    return lines


def resolution(rows: list[dict]) -> list[str]:
    """The smallest rho step from a tie the decoder resolves, per tract and noise."""
    lines = [
        "#### Resolution of the shortest-of-equally-dominant regime",
        "",
        "For each tract and noise level: the decoded ratio rho-hat = share(L-1)/share(L) "
        "against the true rho (pair family, both dominants pooled); its bias and spread "
        "near the tie; and the largest rho < 1 at which the dominant is still right in at "
        "least 95 % of decodes (`rho_95`), i.e. how close to equal the decoder separates.",
        "",
        "| tract | marginal | noise | rho-hat bias at rho 0.5-0.9 (median) | rho-hat SD at rho 0.8-1.25 "
        "(median) | rho_95 | right % at rho 0.9 | near-equal % at rho 0.9 |",
        "|---|---|---|---:|---:|---:|---:|---:|",
    ]
    for tract in ("HV2-F", "HV2-R", "HV1-F", "HV1-R", "HV3-F", "HV3-R"):
        keys = (
            ["first", "total"]
            if tract == "HV2-F"
            else ["total" if tract == "HV2-R" else "first"]
        )
        for key in keys:
            for noise in (1.0, 2.0, 4.0):
                by = defaultdict(list)
                for r in rows:
                    if (
                        r["tract"] == tract
                        and r["family"] == "pair"
                        and r["noise"] == noise
                    ):
                        by[r["value"]].append(score(r, key))
                if not by:
                    continue
                bias, sds = [], []
                rho95 = None
                for rho in sorted(by):
                    dec = [s for s in by[rho] if s["decoded"]]
                    rh = [s["rho_hat"] for s in dec if "rho_hat" in s]
                    if 0.5 <= rho <= 0.9 and rh:
                        bias.append(np.median(rh) - rho)
                    if 0.8 <= rho <= 1.25 and len(rh) > 2:
                        sds.append(np.std(rh) / rho)
                    if (
                        0 < rho < 1
                        and dec
                        and np.mean([s["right"] for s in dec]) >= 0.95
                    ):
                        rho95 = rho
                d9 = [s for s in by.get(0.9, []) if s["decoded"]]
                lines.append(
                    f"| {tract} | {key} | {noise:g} | {np.median(bias) if bias else float('nan'):+.3f} | "
                    f"{np.median(sds) if sds else float('nan'):.3f} | "
                    f"{rho95 if rho95 is not None else '-'} | {pct([s['right'] for s in d9])} | "
                    f"{pct([s['margin'] == 'near_equal' for s in d9])} |"
                )
    lines.append("")
    return lines


EXPECTED_MARK = {
    "HV2-F": "T",
    "HV2-R": "A",
    "HV1-F": None,
    "HV1-R": None,
    "HV3-F": "A",  # C6 A C: the decoder may take 574A as a mark with a run of one
    "HV3-R": "T",
}


def headline(rows: list[dict]) -> list[str]:
    """One row per tract and noise level over the pair family (rho != 0)."""
    lines = [
        "#### Headline (pair family, rho in 0.2-2, both dominants)",
        "",
        "`structure right`: the chosen mark is the tract's (T for HV2-F, A for HV2-R, none "
        "for the single runs). `right | structure`: the dominant length is right among "
        "decodes whose structure is right. `clear, wrong`: a clear margin on a wrong "
        "dominant (rho 1 excluded). `rho-hat/rho`: median at rho 0.5-0.9. `tie`: decodes "
        "at rho 1 whose dominant is the shorter length. `undecoded, primaries right`: "
        "reads with no decode whose own run length is the dominant (single runs and HV2-F).",
        "",
        "| tract | marginal | noise | decoded % | structure right % | dominant right % | "
        "right given structure % | clear, wrong % | rho-hat/rho | tie -> shorter % | "
        "undecoded, primaries right |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for tract in ("HV2-F", "HV2-R", "HV1-F", "HV1-R", "HV3-F", "HV3-R"):
        key = "total" if tract in ("HV2-F", "HV2-R") else "first"
        for noise in (1.0, 2.0, 4.0):
            group = [
                r
                for r in rows
                if r["tract"] == tract
                and r["family"] == "pair"
                and r["noise"] == noise
                and r["value"] > 0
            ]
            if not group:
                continue
            dec = [r for r in group if r.get("decoded")]
            structure = [r.get("mark") == EXPECTED_MARK[tract] for r in dec]
            right, right_s, clear_wrong, ratio, tie = [], [], [], [], []
            for r, ok in zip(dec, structure, strict=True):
                s = score(r, key)
                if r["value"] == 1.0:
                    tie.append(s["shorter"])
                    continue
                right.append(s["right"])
                if ok:
                    right_s.append(s["right"])
                clear_wrong.append(s["margin"] == "clear" and not s["right"])
                if 0.5 <= r["value"] <= 0.9 and "rho_hat" in s:
                    ratio.append(s["rho_hat"] / r["value"])
            undecoded = [r for r in group if not r.get("decoded")]
            primaries = "-"
            if undecoded and tract != "HV2-R":
                hits = 0
                for r in undecoded:
                    truth = as_int(r["truth_first"])
                    hits += r.get("primary_run") == max(truth, key=truth.get)
                primaries = f"{hits}/{len(undecoded)}"
            lines.append(
                f"| {tract} | {key} | {noise:g} | {pct([bool(r.get('decoded')) for r in group])} | "
                f"{pct(structure)} | {pct(right)} | {pct(right_s)} | {pct(clear_wrong)} | "
                f"{np.median(ratio) if ratio else float('nan'):.2f} | {pct(tie)} | {primaries} |"
            )
    lines.append("")
    return lines


def ablation_table(rows: list[dict]) -> list[str]:
    lines = [
        "#### Ablation: what biases the shares near a tie",
        "",
        "Pair family, noise x1, dominant 9 (HV2) or 11 (HV1). rho-hat/rho is the median "
        "decoded ratio over the true one at rho 0.5-0.9; `tie -> shorter` is the share of "
        "decodes at rho 1 whose dominant is the shorter length; `wrong mark` the share of "
        "decodes whose chosen structure is not the tract's.",
        "",
        "| tract | marginal | switched off | rho-hat / rho | tie -> shorter % | right % (rho != 1) "
        "| wrong mark % |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    by = defaultdict(list)
    for r in rows:
        by[(r["tract"], "+".join(r.get("ablate", [])) or "nothing")].append(r)
    for (tract, ablate), group in sorted(by.items()):
        key = "total" if tract in ("HV2-R", "HV2-F") else "first"
        dec = [r for r in group if r.get("decoded")]
        ratios, ties, right = [], [], []
        for r in dec:
            s = score(r, key)
            if 0.5 <= r["value"] <= 0.9 and "rho_hat" in s:
                ratios.append(s["rho_hat"] / r["value"])
            if r["value"] == 1.0:
                ties.append(s["shorter"])
            else:
                right.append(s["right"])
        wrong_mark = [r.get("mark") != EXPECTED_MARK[tract] for r in dec]
        lines.append(
            f"| {tract} | {key} | {ablate} | {np.median(ratios):.3f} | {pct(ties)} | "
            f"{pct(right)} | {pct(wrong_mark)} |"
        )
    lines.append("")
    return lines


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    rows = [json.loads(line) for line in (args.out / "cases.jsonl").open()]
    errors = [r for r in rows if "error" in r]
    rows = [r for r in rows if "error" not in r]
    noise_levels = defaultdict(list)
    for r in rows:
        noise_levels[r["noise"]].append(r["second_q50"])
    lines = [
        f"# Tract decode benchmark ({len(rows)} cases, {len(errors)} failed)",
        "",
        "Noise x1 is the fitted read's own noise; the pre-tract median second-peak share "
        "per level: "
        + ", ".join(
            f"x{k:g}: {np.median(v):.3f}" for k, v in sorted(noise_levels.items())
        ),
        "",
    ]
    if any("ablate" in r for r in rows):
        lines += ablation_table(rows)
        (args.out / "report.md").write_text("\n".join(lines) + "\n")
        print("\n".join(lines))
        return
    lines += headline(rows)
    lines += resolution(rows)
    for tract in ("HV2-F", "HV2-R", "HV1-F", "HV1-R", "HV3-F", "HV3-R"):
        keys = (
            ["first", "total"]
            if tract == "HV2-F"
            else ["total" if tract == "HV2-R" else "first"]
        )
        for key in keys:
            lines += pair_table(rows, tract, key)
    lines += family_table(rows, "spread")
    lines += family_table(rows, "stutter")
    # HV2-R's first run (fixed C6 of 315.1C) must read 6 everywhere.
    fixed = [r for r in rows if r["tract"] == "HV2-R" and r.get("decoded")]
    if fixed:
        lines.append(
            f"HV2-R first-run marginal (fixed at {bench.TRACTS['HV2-R'].fixed}): dominant right in "
            f"{pct([r['first']['dominant'] == bench.TRACTS['HV2-R'].fixed for r in fixed])} % of "
            f"{len(fixed)} decodes."
        )
    marks = defaultdict(lambda: defaultdict(int))
    for r in rows:
        if r.get("decoded"):
            marks[r["tract"]][str(r.get("mark"))] += 1
    lines.append("")
    lines.append(
        "Structure chosen (mark): "
        + "; ".join(
            f"{t}: " + ", ".join(f"{m} {n}" for m, n in sorted(c.items()))
            for t, c in sorted(marks.items())
        )
    )
    (args.out / "report.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
