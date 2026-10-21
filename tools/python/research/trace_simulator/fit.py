# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy", "scipy"]
# ///
"""Fits the electropherogram forward model's parameters on real traces (research only).

For each trace, over its clean calls (vendor QV >= 30, one channel holding at least
`CLEAN_SHARE` of the four at the call position), it measures:

- peak shape: a Gaussian per isolated peak (no neighbour of the same base) by a
  weighted log-parabola, its width vs call index, its asymmetry, and on a subset an
  exponentially modified Gaussian (EMG) tail;
- spacing: a cubic of peak centre vs call index, its derivative the local spacing;
- dye mobility: each dye's median centre residual from that cubic;
- placement: the vendor position minus the fitted centre, per dye;
- amplitude: the log-height envelope vs call index, each dye's gain, the context
  (previous, own and next base) effect pooled later, the residual's spread and
  lag-1 autocorrelation;
- noise: each channel's floor and robust spread away from its own peaks, and the
  floor's sample autocorrelation;
- cross-talk: a dye's value at another dye's isolated apex.

Usage:
    uv run --script fit.py MANIFEST.tsv OUT.jsonl [--samples N] [--seed S] [--workers W]

Writes one JSON object per trace; `aggregate.py` reduces them to distributions.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import exponnorm

sys.path.insert(0, str(Path(__file__).resolve().parent))
import abif  # noqa: E402
import reference  # noqa: E402

BASES = "ACGT"
CLEAN_QV = 30
CLEAN_SHARE = 0.75
EMG_PEAKS = 40
PRIMERS = ("HV1F", "HV1R", "HV2F", "HV3R")


def log_parabola(
    y: np.ndarray, apex: int, half: int
) -> tuple[float, float, float] | None:
    """Centre, sigma, height of the Gaussian through the samples above 25 % of the apex."""
    top = y[apex]
    if top <= 0:
        return None
    lo, hi = apex, apex
    while lo > apex - half and lo > 0 and y[lo - 1] > 0.25 * top and y[lo - 1] <= y[lo]:
        lo -= 1
    while (
        hi < apex + half
        and hi + 1 < len(y)
        and y[hi + 1] > 0.25 * top
        and y[hi + 1] <= y[hi]
    ):
        hi += 1
    if hi - lo < 2:
        return None
    t = np.arange(lo, hi + 1, dtype=float) - apex
    v = y[lo : hi + 1]
    w = v  # Caruana weighting: weights y^2 on log y
    design = np.stack([np.ones_like(t), t, t * t], axis=1) * w[:, None]
    coef, *_ = np.linalg.lstsq(design, np.log(v) * w, rcond=None)
    a, b, c = coef
    if c >= 0:
        return None
    sigma = math.sqrt(-1.0 / (2.0 * c))
    centre = -b / (2.0 * c)
    if abs(centre) > 2.0:
        return None
    height = math.exp(a - b * b / (4.0 * c))
    return apex + centre, sigma, height


def half_widths(y: np.ndarray, apex: int, half: int) -> tuple[float, float] | None:
    """Left and right half widths at half maximum, linearly interpolated."""
    top = y[apex] / 2.0
    left = right = None
    for j in range(apex, max(apex - half, 0), -1):
        if y[j - 1] < top <= y[j]:
            left = apex - (j - (y[j] - top) / (y[j] - y[j - 1]))
            break
    for j in range(apex, min(apex + half, len(y) - 1)):
        if y[j + 1] < top <= y[j]:
            right = (j + (y[j] - top) / (y[j] - y[j + 1])) - apex
            break
    if left is None or right is None:
        return None
    return left, right


def emg(
    t: np.ndarray, area: float, centre: float, sigma: float, tau: float
) -> np.ndarray:
    return area * exponnorm.pdf(t, tau / sigma, loc=centre, scale=sigma)


def gauss(t: np.ndarray, height: float, centre: float, sigma: float) -> np.ndarray:
    return height * np.exp(-0.5 * ((t - centre) / sigma) ** 2)


def robust_poly(x: np.ndarray, y: np.ndarray, degree: int) -> np.ndarray:
    keep = np.ones(len(x), dtype=bool)
    coef = np.polyfit(x, y, degree)
    for _ in range(3):
        coef = np.polyfit(x[keep], y[keep], degree)
        resid = y - np.polyval(coef, x)
        scale = 1.4826 * np.median(np.abs(resid[keep] - np.median(resid[keep]))) + 1e-9
        keep = np.abs(resid) < 4 * scale
    return coef


def mad_sd(values: np.ndarray) -> float:
    return float(1.4826 * np.median(np.abs(values - np.median(values))))


def lag1(values: np.ndarray) -> float:
    if len(values) < 10:
        return float("nan")
    v = values - values.mean()
    den = float(np.dot(v, v))
    return float(np.dot(v[1:], v[:-1]) / den) if den > 0 else float("nan")


def fit_trace(path: str) -> dict:
    trace = abif.read(path)
    y = trace.channels
    calls = trace.pbas
    n = min(len(calls), len(trace.ploc), len(trace.pcon))
    ploc = trace.ploc[:n]
    pcon = trace.pcon[:n]
    bases = np.array([BASES.index(c) if c in BASES else -1 for c in calls[:n]])
    at = y[:, np.clip(ploc, 0, y.shape[1] - 1)].T  # (n, 4)
    total = at.sum(axis=1)
    share = np.where(
        total > 0, at[np.arange(n), np.clip(bases, 0, 3)] / np.maximum(total, 1e-9), 0
    )
    clean = (pcon >= CLEAN_QV) & (bases >= 0) & (share >= CLEAN_SHARE)
    clean[:10] = False
    out: dict = {"path": path, "calls": int(n), "meta": trace.meta}
    serial = re.search(r"SerialNumber=([^;]+)", trace.meta.get("instrument", ""))
    out["instrument"] = serial.group(1) if serial else ""
    primer = next((p for p in PRIMERS if p in Path(path).name), "")
    out["primer"] = primer
    placement = reference.place(calls[:n])
    if placement:
        out["placement"] = {
            "strand": placement.strand,
            "start": placement.start,
            "votes": placement.votes,
        }
    good = np.flatnonzero(pcon >= 20)
    if len(good) < 50 or clean.sum() < 80:
        out["status"] = "too_few_clean_calls"
        return out
    out["qv20_first"], out["qv20_last"] = int(good[0]), int(good[-1])
    out["clean_calls"] = int(clean.sum())
    spacing0 = float(np.median(np.diff(ploc)))
    half = max(3, int(round(0.8 * spacing0)))

    # Isolated peaks: clean, neighbours clean and of other bases.
    peaks = []
    for i in range(2, n - 2):
        b = bases[i]
        if not (clean[i] and clean[i - 1] and clean[i + 1]):
            continue
        if bases[i - 1] == b or bases[i + 1] == b:
            continue
        channel = y[b]
        p = int(ploc[i])
        lo, hi = max(p - 2, 1), min(p + 3, len(channel) - 1)
        apex = lo + int(np.argmax(channel[lo:hi]))
        fitted = log_parabola(channel, apex, half)
        if fitted is None:
            continue
        hw = half_widths(channel, apex, half)
        peaks.append((i, b, apex, *fitted, *(hw if hw else (np.nan, np.nan))))
    if len(peaks) < 60:
        out["status"] = "too_few_isolated_peaks"
        return out
    arr = np.array(peaks, dtype=float)
    idx, base, apexes, centre, sigma, height, left, right = arr.T
    base = base.astype(int)
    x = idx / 1000.0
    pos_coef = robust_poly(x, centre, 3)
    # Position residual from a local line (+-20 calls): sequence-dependent
    # mobility shows as peak-to-peak spacing irregularity, not as drift.
    resid = np.empty(len(idx))
    for k in range(len(idx)):
        near = np.abs(idx - idx[k]) <= 20
        line = np.polyfit(idx[near], centre[near], 1)
        resid[k] = centre[k] - np.polyval(line, idx[k])
    dye = np.array(
        [np.median(resid[base == b]) if np.any(base == b) else 0.0 for b in range(4)]
    )
    dye_centered = dye - dye.mean()
    jitter = resid - dye[base]
    deriv = np.polyder(pos_coef)
    lo_i, hi_i = int(idx.min()), int(idx.max())
    out["status"] = "ok"
    out["isolated_peaks"] = len(peaks)
    out["span"] = [lo_i, hi_i]
    out["position_poly"] = pos_coef.tolist()  # centre (samples) vs call/1000
    out["spacing_at"] = {
        str(k): float(np.polyval(deriv, k / 1000.0) / 1000.0)
        for k in (50, 100, 200, 300, 400, 500, 600)
        if lo_i <= k <= hi_i
    }
    out["dye_offset"] = dye_centered.tolist()  # samples
    out["position_jitter_sd"] = mad_sd(jitter)
    out["position_resid_lag1"] = lag1(jitter)
    position_context: dict[str, list[float]] = {}
    for k in range(len(idx)):
        i = int(idx[k])
        if bases[i - 1] >= 0 and bases[i + 1] >= 0:
            cell = position_context.setdefault(
                calls[i - 1] + calls[i] + calls[i + 1], [0.0, 0.0, 0]
            )
            cell[0] += float(jitter[k])
            cell[1] += float(jitter[k]) ** 2
            cell[2] += 1
    out["position_context"] = position_context
    out["position_resid_calls"] = [
        [int(i), round(float(v), 3)] for i, v in zip(idx, jitter, strict=True)
    ]
    vendor = ploc[idx.astype(int)] - centre
    out["vendor_offset"] = [
        float(np.median(vendor[base == b])) if np.any(base == b) else 0.0
        for b in range(4)
    ]
    out["vendor_jitter_sd"] = mad_sd(vendor - np.array(out["vendor_offset"])[base])
    # Width vs index, in samples.
    width_coef = robust_poly(x, sigma, 1)
    out["width_line"] = width_coef.tolist()  # sigma (samples) vs call/1000
    out["width_at"] = {
        k: float(np.polyval(width_coef, int(k) / 1000.0)) for k in out["spacing_at"]
    }
    out["width_dye"] = [
        float(np.median(sigma[base == b] - np.polyval(width_coef, x[base == b])))
        if np.any(base == b)
        else 0.0
        for b in range(4)
    ]
    ok = np.isfinite(left) & np.isfinite(right)
    out["asymmetry"] = float(np.median((right[ok] - left[ok]) / (right[ok] + left[ok])))
    # EMG on a spread subset.
    order = np.argsort(height)[::-1]
    chosen = sorted(order[: max(EMG_PEAKS * 3, 1)], key=lambda k: idx[k])[::3][
        :EMG_PEAKS
    ]
    taus, gains = [], []
    for k in chosen:
        b, apex = base[k], int(apexes[k])
        span = int(round(1.5 * spacing0))
        t = np.arange(
            max(apex - span, 0), min(apex + span + 1, y.shape[1]), dtype=float
        )
        v = y[b, t.astype(int)]
        try:
            g, _ = curve_fit(gauss, t, v, p0=[v.max(), apex, sigma[k]], maxfev=2000)
            gr = float(np.sum((v - gauss(t, *g)) ** 2))
            area0 = v.max() * sigma[k] * math.sqrt(2 * math.pi)
            e, _ = curve_fit(
                emg,
                t,
                v,
                p0=[area0, apex - 0.5, sigma[k] * 0.9, 0.5],
                bounds=([0, apex - 5, 0.3, 0.02], [np.inf, apex + 5, 20, 20]),
                maxfev=4000,
            )
            er = float(np.sum((v - emg(t, *e)) ** 2))
        except (RuntimeError, ValueError):
            continue
        taus.append((e[2], e[3]))
        gains.append(1.0 - er / gr if gr > 0 else 0.0)
    if taus:
        t_arr = np.array(taus)
        out["emg_sigma"] = float(np.median(t_arr[:, 0]))
        out["emg_tau"] = float(np.median(t_arr[:, 1]))
        out["emg_rss_gain"] = float(np.median(gains))
    # Heights: every clean call's own channel apex near its position.
    ci = np.flatnonzero(clean)
    h = np.array(
        [y[bases[i], max(int(ploc[i]) - 2, 0) : int(ploc[i]) + 3].max() for i in ci],
        dtype=float,
    )
    keep = h > 0
    ci, h = ci[keep], h[keep]
    logh = np.log(h)
    env_coef = robust_poly(ci / 1000.0, logh, 4)
    r = logh - np.polyval(env_coef, ci / 1000.0)
    gain = np.array(
        [
            np.median(r[bases[ci] == b]) if np.any(bases[ci] == b) else 0
            for b in range(4)
        ]
    )
    gain -= gain.mean()
    r2 = r - gain[bases[ci]]
    out["envelope_poly"] = env_coef.tolist()  # log height vs call/1000
    grid = np.arange(ci.min(), ci.max() + 1)
    env = np.polyval(env_coef, grid / 1000.0)
    out["envelope_peak_call"] = int(grid[np.argmax(env)])
    out["envelope_peak_height"] = float(math.exp(env.max()))
    out["height_median"] = float(np.median(h))

    def env_at(k: int) -> float | None:
        return (
            float(math.exp(np.polyval(env_coef, k / 1000.0)))
            if ci.min() <= k <= ci.max()
            else None
        )

    out["envelope_at"] = {str(k): env_at(k) for k in (50, 100, 200, 300, 400, 500, 600)}
    out["dye_gain"] = gain.tolist()  # log
    # Context sums keyed by previous, own, next base, for pooling.
    context: dict[str, list[float]] = {}
    for i, value in zip(ci, r2, strict=True):
        if 0 < i < n - 1 and bases[i - 1] >= 0 and bases[i + 1] >= 0:
            key = calls[i - 1] + calls[i] + calls[i + 1]
            cell = context.setdefault(key, [0.0, 0.0, 0])
            cell[0] += float(value)
            cell[1] += float(value) ** 2
            cell[2] += 1
    out["context"] = context
    out["height_resid_sd"] = mad_sd(r2)
    out["height_resid_calls"] = [
        [int(i), round(float(v), 4)] for i, v in zip(ci, r2, strict=True)
    ]
    consecutive = np.diff(ci) == 1
    pairs = np.flatnonzero(consecutive)
    if len(pairs) > 20:
        a_, b_ = r2[pairs], r2[pairs + 1]
        out["height_resid_lag1"] = float(np.corrcoef(a_, b_)[0, 1])
    # Noise floor per channel away from its own peaks, inside the clean span.
    s_lo, s_hi = int(ploc[ci.min()]), int(ploc[ci.max()])
    floors, sds, lags, white = [], [], [], []
    for c in range(4):
        mask = np.ones(y.shape[1], dtype=bool)
        mask[:s_lo] = False
        mask[s_hi:] = False
        for p in ploc[bases == c]:
            mask[max(int(p - 1.5 * spacing0), 0) : int(p + 1.5 * spacing0) + 1] = False
        values = y[c, mask]
        steps = np.flatnonzero(mask[1:] & mask[:-1])
        white.append(
            mad_sd(np.diff(y[c])[steps]) / math.sqrt(2)
            if len(steps) > 50
            else float("nan")
        )
        floors.append(float(np.median(values)) if len(values) else float("nan"))
        sds.append(mad_sd(values) if len(values) else float("nan"))
        runs = np.flatnonzero(mask[1:] & mask[:-1])
        if len(runs) > 50:
            v = y[c] - (floors[-1] if math.isfinite(floors[-1]) else 0)
            lags.append(float(np.corrcoef(v[runs], v[runs + 1])[0, 1]))
    out["floor"] = floors
    out["floor_sd"] = sds
    out["floor_white_sd"] = white
    out["floor_lag1"] = float(np.nanmedian(lags)) if lags else float("nan")
    out["floor_zero_fraction"] = float(np.mean(y[:, s_lo:s_hi] <= 0))
    # Cross-talk at isolated apexes.
    xt = np.full((4, 4), np.nan)
    for b in range(4):
        for c in range(4):
            if b == c:
                continue
            vals = []
            for k in np.flatnonzero(base == b):
                i = int(idx[k])
                if bases[i - 1] == c or bases[i + 1] == c:
                    continue
                vals.append(y[c, int(apexes[k])] / max(y[b, int(apexes[k])], 1))
            if len(vals) > 10:
                xt[b, c] = float(np.median(vals))
    out["crosstalk"] = [[None if math.isnan(v) else v for v in row] for row in xt]
    # Realism statistic: second share at clean call positions.
    srt = np.sort(at[clean], axis=1)
    ratio = srt[:, 2] / np.maximum(srt[:, 3], 1e-9)
    out["second_ratio_q"] = np.percentile(ratio, [10, 50, 90]).tolist()
    out["primary_share_q"] = np.percentile(share[clean], [10, 50, 90]).tolist()
    return out


def select(manifest: Path, samples: int, seed: int) -> list[tuple[str, dict]]:
    rows = []
    with manifest.open() as handle:
        lines = (line for line in handle if not line.startswith("#"))
        for row in csv.DictReader(lines, delimiter="\t"):
            if row.get("retired", "0") != "0":
                continue
            rows.append(row)
    rng = random.Random(seed)
    recent = [r for r in rows if r["batch_date"] >= "2026-01"]
    older = [r for r in rows if r["batch_date"] < "2026-01"]
    picked = rng.sample(recent, min(len(recent), samples * 2 // 3)) + rng.sample(
        older, min(len(older), samples - samples * 2 // 3)
    )
    jobs = []
    for row in picked:
        paths = row["trace_paths"].split(";")
        for primer in PRIMERS:
            match = sorted(p for p in paths if primer in Path(p).name)
            if match:
                info = {
                    k: row[k] for k in ("sample_id", "batch", "batch_date", "stratum")
                }
                jobs.append((match[-1], info))
    return jobs


def run(job: tuple[str, dict]) -> dict:
    path, info = job
    try:
        out = fit_trace(path)
    except Exception as error:  # noqa: BLE001 - one failed trace must not stop the fit
        out = {"path": path, "status": f"error: {error}"}
    out.update(info)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--samples", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20261021)
    parser.add_argument("--workers", type=int, default=32)
    args = parser.parse_args()
    jobs = select(args.manifest, args.samples, args.seed)
    with ProcessPoolExecutor(args.workers) as pool, args.out.open("w") as handle:
        for result in pool.map(run, jobs, chunksize=4):
            handle.write(json.dumps(result) + "\n")
    print(f"{len(jobs)} traces -> {args.out}")


if __name__ == "__main__":
    main()
