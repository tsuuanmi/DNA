"""Physical forward model of a Sanger electropherogram and its renderer (research only).

A read is a mixture of templates. Every template's fragment that ends `n` bases after
the primer migrates to the same place, so a length mixture shares one peak grid and
differs only in which base each template shows at each grid point. Per peak:

    centre(k, n) = grid(n) + dye_offset[b] + position_context[ctx] + det(ctx7) + iid
    log height   = envelope(n) + dye_gain[b] + height_context[ctx] + det(ctx7) + iid + log f_k
    shape        = Gaussian, sigma(n) = width line (EMG tails were fitted and found absent)

`ctx` is the template's previous, own and next base; `ctx7` its six bases up to and
including the next one. The `det` terms are deterministic hashes of that context: the
real residuals recur at one reference position across samples (intraclass correlation
0.77 for height, 0.85 for position), so shifted copies of one downstream sequence
carry one shifted pattern. The analysed channels add a background of minor peaks on
the grid (lognormal, calibrated to the clean-call second/primary distribution),
spectral cross-talk, a slowly varying floor and white noise, and are clipped at zero
as the instrument's analysed data are.

Calls are produced by an oracle-grid basecaller: one call per grid point (no merged
or split calls), at the apex of the strongest channel within a third of a spacing of
the grid point, its base the strongest channel there, its quality from that channel's
share. This is what `dna` reads as PLOC/PBAS/PCON.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

BASES = "ACGT"
# Real reads' median vendor quality by the primary's share at the call (80 traces).
QV_SHARE = (0.25, 0.6, 0.75, 0.85, 0.925, 0.96)
QV_VALUE = (1.0, 9.0, 22.0, 33.0, 51.0, 55.0)
# The background terms calibrated against real reads (realism.py): minor peaks on the
# grid (log level and spread, relative to the envelope), and the floor's scale, shift
# (RFU, so about a fifth of the analysed samples sit at zero) and smoothing (samples).
CALIBRATION = {
    "minor_mu": -5.3,
    "minor_sd": 1.4,
    "floor_scale": 0.6,
    "floor_shift": -3.0,
    "floor_len": 4.0,
    "position_smooth": 2.5,
}


@dataclass
class Params:
    """One read's physical parameters (samples, RFU, natural-log heights)."""

    first_sample: float = 30.0
    spacing0: float = 12.1  # samples per call at call 0
    spacing_slope: float = -0.0006  # per call
    sigma0: float = 2.95  # Gaussian sigma (samples) at call 0
    sigma_slope: float = 0.0  # per call
    dye_offset: tuple[float, ...] = (-0.23, -0.05, 0.24, 0.05)
    envelope: tuple[tuple[int, float], ...] = (
        (0, 6.6),
        (80, 6.8),
        (300, 6.1),
        (700, 4.9),
    )
    dye_gain: tuple[float, ...] = (-0.06, -0.006, 0.047, 0.017)
    height_context: dict[str, float] = field(default_factory=dict)
    position_context: dict[str, float] = field(default_factory=dict)
    height_det_sd: float = 0.10
    height_iid_sd: float = 0.055
    position_det_sd: float = 0.77
    position_iid_sd: float = 0.32
    position_smooth: float = 2.5  # calls; the mobility pattern's correlation length
    minor_mu: float = math.log(0.006)  # log of a minor peak's height / envelope
    minor_sd: float = 1.0
    crosstalk: tuple[tuple[float, ...], ...] = ((0,) * 4,) * 4
    floor: tuple[float, ...] = (3.0, 4.0, 2.0, 2.0)
    floor_sd: float = 4.0  # slowly varying floor (RFU)
    floor_len: float = 3.0  # its Gaussian smoothing (samples)
    white_sd: float = 1.2  # RFU
    seed: int = 0
    pattern_seed: int = (
        7  # the sequence-determined pattern's seed (shared across reads)
    )
    # The read start: real reads open with unresolved signal (primer, unincorporated
    # dye) before the first good call. A noise burst over the first calls stands in
    # for it; without it every penalty of a clean read is relative to a clean worst.
    start_calls: int = 15
    start_noise: float = 0.4  # relative to the envelope at call 0

    def scaled_noise(self, factor: float) -> Params:
        """The same read with its noise (white, floor, minor peaks) times `factor`."""
        return replace(
            self,
            white_sd=self.white_sd * factor,
            floor_sd=self.floor_sd * factor,
            minor_mu=self.minor_mu + math.log(factor) if factor > 0 else -50.0,
        )

    def ablated(self, parts: list[str]) -> Params:
        """The read with model parts switched off, to attribute an effect to them.

        `noise`: no white noise, floor, minor peaks or cross-talk; `heights`: every peak
        at its envelope (no dye gain, context or residual); `position`: every peak on
        the grid (no dye offset, context or residual).
        """
        out = self
        if "noise" in parts:
            out = replace(out.scaled_noise(0.0), crosstalk=((0.0,) * 4,) * 4)
        if "heights" in parts:
            out = replace(
                out,
                height_context={},
                dye_gain=(0.0,) * 4,
                height_det_sd=0.0,
                height_iid_sd=0.0,
            )
        if "position" in parts:
            out = replace(
                out,
                position_context={},
                dye_offset=(0.0,) * 4,
                position_det_sd=0.0,
                position_iid_sd=0.0,
            )
        return out

    def grid(self, n: np.ndarray) -> np.ndarray:
        return self.first_sample + self.spacing0 * n + 0.5 * self.spacing_slope * n * n

    def spacing(self, n: np.ndarray) -> np.ndarray:
        return self.spacing0 + self.spacing_slope * n

    def sigma(self, n: np.ndarray) -> np.ndarray:
        return np.maximum(self.sigma0 + self.sigma_slope * n, 0.8)

    def log_envelope(self, n: np.ndarray) -> np.ndarray:
        knots = np.array([k for k, _ in self.envelope], dtype=float)
        values = np.array([v for _, v in self.envelope], dtype=float)
        out = np.interp(n, knots, values)
        if len(knots) >= 2:
            slope = min((values[-1] - values[-2]) / (knots[-1] - knots[-2]), 0.0)
            beyond = n > knots[-1]
            out[beyond] = values[-1] + slope * (n[beyond] - knots[-1])
        return out


def _hash_normal(text: str, seed: int) -> float:
    """A standard normal value determined by `text` and `seed`."""
    digest = hashlib.blake2b(f"{seed}:{text}".encode(), digest_size=16).digest()
    u1 = (int.from_bytes(digest[:8], "big") + 0.5) / 2**64
    u2 = (int.from_bytes(digest[8:], "big") + 0.5) / 2**64
    return math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)


def _smooth(values: np.ndarray, width: float) -> np.ndarray:
    """`values` smoothed by a Gaussian of `width` calls, rescaled to unit variance."""
    if width <= 0:
        return values
    k = np.arange(-int(3 * width) - 1, int(3 * width) + 2)
    kernel = np.exp(-0.5 * (k / width) ** 2)
    kernel /= math.sqrt(np.sum(kernel**2))
    return np.convolve(values, kernel, mode="same")


@dataclass
class Artefact:
    """An optional artefact: `spike`, `blob` or `noise` (see `render`)."""

    kind: str
    call: float = 0.0  # where (grid call index)
    channel: int = 0  # blob channel
    height: float = 1.0  # relative to the local envelope
    width: float = 3.0  # blob width (calls, full width at half maximum)


@dataclass
class Rendered:
    channels: np.ndarray  # (4, samples), A/C/G/T
    ploc: np.ndarray
    pbas: str
    pcon: np.ndarray
    truth_grid: np.ndarray  # the grid position of each call


def render(
    templates: list[tuple[str, float]],
    params: Params,
    calls: int | None = None,
    artefacts: list[Artefact] | None = None,
) -> Rendered:
    """Renders the analysed channels of a mixture of `templates` (sequence, fraction)."""
    rng = np.random.default_rng(params.seed)
    total = sum(f for _, f in templates)
    count = calls or min(len(t) for t, _ in templates)
    n = np.arange(count, dtype=float)
    grid = params.grid(n)
    sigma = params.sigma(n)
    env = params.log_envelope(n)
    samples = int(grid[-1] + 6 * params.spacing0)
    t = np.arange(samples, dtype=float)
    y = np.zeros((4, samples))

    def add(channel: int, centre: float, height: float, width: float) -> None:
        lo = max(int(centre - 5 * width), 0)
        hi = min(int(centre + 5 * width) + 2, samples)
        if hi <= lo:
            return
        y[channel, lo:hi] += height * np.exp(-0.5 * ((t[lo:hi] - centre) / width) ** 2)

    for template, fraction in templates:
        if fraction <= 0:
            continue
        share = fraction / total
        padded = "N" + template + "N" * 8
        # Sequence-determined mobility: a hashed value per context, smoothed along the
        # template so neighbouring peaks drift together (real lag-1 about 0.47).
        hashed = np.array(
            [
                _hash_normal("p" + padded[max(i - 4, 0) : i + 3], params.pattern_seed)
                for i in range(count)
            ]
        )
        drift = _smooth(hashed, params.position_smooth)
        for i in range(count):
            base = template[i]
            b = BASES.find(base)
            if b < 0:
                continue
            tri = padded[i : i + 3]
            ctx7 = padded[max(i - 4, 0) : i + 3]
            centre = (
                grid[i]
                + params.dye_offset[b]
                + params.position_context.get(tri, 0.0)
                + params.position_det_sd * drift[i]
                + params.position_iid_sd * rng.standard_normal()
            )
            log_h = (
                env[i]
                + params.dye_gain[b]
                + params.height_context.get(tri, 0.0)
                + params.height_det_sd * _hash_normal("h" + ctx7, params.pattern_seed)
                + params.height_iid_sd * rng.standard_normal()
            )
            height = share * math.exp(log_h)
            add(b, centre, height, sigma[i])
            for c in range(4):
                leak = params.crosstalk[b][c]
                if c != b and leak > 0:
                    add(c, centre, height * leak, sigma[i])
    # Background of minor peaks on the grid, every channel.
    minor = np.exp(params.minor_mu + params.minor_sd * rng.standard_normal((count, 4)))
    for i in range(count):
        level = math.exp(env[i])
        for c in range(4):
            add(c, grid[i] + params.dye_offset[c], level * minor[i, c], sigma[i])
    artefacts = list(artefacts or [])
    if params.start_calls > 0 and params.start_noise > 0:
        artefacts.append(
            Artefact(
                "noise",
                call=-0.5,
                height=params.start_noise,
                width=params.start_calls + 0.5,
            )
        )
    for artefact in artefacts:
        i = int(min(max(artefact.call, 0), count - 1))
        level = math.exp(env[i])
        centre = float(params.grid(np.array([artefact.call]))[0])
        if artefact.kind == "spike":
            for c in range(4):
                add(c, centre, level * artefact.height, 0.9)
        elif artefact.kind == "blob":
            width = artefact.width * params.spacing0 / 2.355
            add(artefact.channel, centre, level * artefact.height, width)
        elif artefact.kind == "noise":
            lo = int(params.grid(np.array([artefact.call]))[0])
            hi = int(params.grid(np.array([artefact.call + artefact.width]))[0])
            y[:, lo:hi] += (
                level * artefact.height * np.abs(rng.standard_normal((4, hi - lo)))
            )
        else:
            raise ValueError(f"unknown artefact {artefact.kind}")
    # Floor: a slowly varying component plus white noise, clipped at zero.
    if params.floor_sd > 0:
        k = np.arange(-int(4 * params.floor_len), int(4 * params.floor_len) + 1)
        kernel = np.exp(-0.5 * (k / params.floor_len) ** 2)
        kernel /= math.sqrt(np.sum(kernel**2))
        for c in range(4):
            slow = np.convolve(rng.standard_normal(samples), kernel, mode="same")
            y[c] += params.floor[c] + params.floor_sd * slow
    y += params.white_sd * rng.standard_normal(y.shape)
    y = np.clip(y, 0, 32767)
    ploc, pbas, pcon = basecall(y, grid, params.spacing(n))
    return Rendered(y, ploc, pbas, pcon, grid)


def basecall(
    y: np.ndarray, grid: np.ndarray, spacing: np.ndarray
) -> tuple[np.ndarray, str, np.ndarray]:
    """Oracle-grid basecaller: one call per grid point at the strongest apex near it."""
    top = y.max(axis=0)
    ploc = np.empty(len(grid), dtype=np.int64)
    for i, (g, s) in enumerate(zip(grid, spacing, strict=True)):
        lo = max(int(math.floor(g - s / 3)), 0)
        hi = min(int(math.ceil(g + s / 3)) + 1, y.shape[1])
        ploc[i] = lo + int(np.argmax(top[lo:hi]))
    ploc = np.maximum.accumulate(ploc)
    at = y[:, ploc]
    best = at.argmax(axis=0)
    share = at.max(axis=0) / np.maximum(at.sum(axis=0), 1e-9)
    # The vendor's median quality per primary share, measured on real reads.
    pcon = np.rint(np.interp(share, QV_SHARE, QV_VALUE)).astype(int)
    pcon[at.max(axis=0) <= 0] = 1
    pbas = "".join(BASES[b] for b in best)
    return ploc, pbas, pcon


def centred(table: dict[str, list[float]]) -> dict[str, float]:
    """A context table (mean, sd, n per trinucleotide) centred per own (middle) base."""
    out: dict[str, float] = {}
    for base in BASES:
        cells = {k: v for k, v in table.items() if k[1] == base}
        weight = sum(v[2] for v in cells.values())
        mean = sum(v[0] * v[2] for v in cells.values()) / weight if weight else 0.0
        out.update({k: v[0] - mean for k, v in cells.items()})
    return out


def from_fit(
    trace: dict, pooled: dict, seed: int = 0, calibration: dict | None = None
) -> Params:
    """A read's parameters from one fitted trace and the pooled tables."""
    spacing = {int(k): v for k, v in trace["spacing_at"].items()}
    keys = sorted(spacing)
    if len(keys) >= 2:
        slope, intercept = np.polyfit(keys, [spacing[k] for k in keys], 1)
    else:
        slope, intercept = 0.0, spacing[keys[0]]
    width = trace["width_line"]  # sigma vs call/1000
    env = trace["envelope_poly"]
    lo, hi = trace["span"]
    knots = [(int(k), float(np.polyval(env, k / 1000))) for k in range(lo, hi + 1, 25)]
    knots = [(0, knots[0][1])] + knots
    # The trace's dye offsets and gains already hold each base's mean context
    # effect, so the tables are centred per own base.
    height_ctx = centred(pooled["height_context"])
    position_ctx = centred(pooled["position_context"])
    rep = pooled.get("reproducibility", {})
    h_icc = rep.get("height", {}).get("icc", 0.77)
    p_icc = rep.get("position", {}).get("icc", 0.85)
    h_rem = trace["height_resid_sd"] * math.sqrt(
        1 - pooled.get("height_explained", 0.565)
    )
    p_rem = trace["position_jitter_sd"] * math.sqrt(
        1 - pooled.get("position_explained", 0.312)
    )
    white = [w for w in trace.get("floor_white_sd") or [] if w and math.isfinite(w)]
    floor_sd = [w for w in trace.get("floor_sd") or [] if w and math.isfinite(w)]
    xt = [
        [(v or 0.0) if i != j else 0.0 for j, v in enumerate(row)]
        for i, row in enumerate(trace["crosstalk"])
    ]
    cal = {**CALIBRATION, **(calibration or {})}
    return Params(
        first_sample=30.0,
        spacing0=float(intercept),
        spacing_slope=float(slope),
        sigma0=float(width[1]),
        sigma_slope=float(width[0]) / 1000.0,
        dye_offset=tuple(trace["dye_offset"]),
        envelope=tuple(knots),
        dye_gain=tuple(trace["dye_gain"]),
        height_context=height_ctx,
        position_context=position_ctx,
        height_det_sd=h_rem * math.sqrt(h_icc),
        height_iid_sd=h_rem * math.sqrt(1 - h_icc),
        position_det_sd=p_rem * math.sqrt(p_icc),
        position_iid_sd=p_rem * math.sqrt(1 - p_icc),
        minor_mu=cal.get("minor_mu", math.log(0.006)),
        minor_sd=cal.get("minor_sd", 1.0),
        crosstalk=tuple(tuple(row) for row in xt)
        if cal.get("crosstalk", True)
        else ((0.0,) * 4,) * 4,
        floor=tuple(float(v) + cal.get("floor_shift", -2.0) for v in trace["floor"]),
        floor_sd=cal.get("floor_scale", 0.8)
        * (float(np.median(floor_sd)) if floor_sd else 3.0),
        floor_len=cal["floor_len"],
        position_smooth=cal["position_smooth"],
        white_sd=float(np.median(white)) if white else 1.5,
        seed=seed,
    )


def load_pooled(path: Path) -> dict:
    return json.loads(Path(path).read_text())
