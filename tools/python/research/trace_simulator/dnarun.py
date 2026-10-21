"""Runs `dna evidence` on traces and keeps what the benchmark reads (research only).

The binary comes from `$DNA_BIN` (default: `target/release/dna` of the repository);
the configuration is the repository's `config/`, copied into a scratch working
directory so `dna` writes its results and logs there, not into the checkout.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from reference import REPO


def binary() -> str:
    return os.environ.get("DNA_BIN", str(REPO / "target" / "release" / "dna"))


def workdir(root: Path) -> Path:
    """A working directory holding a copy of the repository configuration."""
    root.mkdir(parents=True, exist_ok=True)
    if not (root / "config").exists():
        shutil.copytree(REPO / "config", root / "config")
    return root


def summarise(document: dict) -> dict:
    """The parts of a `dna.read_evidence/v1` document the benchmark reads."""
    read = document["read"]
    calls = document["calls"]
    return {
        "call_count": read["call_count"],
        "trim": read["trim"],
        "informative": read.get("informative"),
        "segments": read["callability"]["segments"],
        "tracts": read.get("tracts", []),
        "mixtures": [
            {
                k: m.get(k)
                for k in ("calls", "run_calls", "explained", "residual", "components")
            }
            for m in read.get("mixtures", [])
        ],
        "bases": "".join(call["base"] for call in calls),
        "profiles": [call.get("profile") or [0.0, 0.0, 0.0, 0.0] for call in calls],
        "second": [
            (call["second"]["share"] if call.get("second") else 0.0) for call in calls
        ],
        "stands_out": [
            bool((call.get("second") or {}).get("stands_out")) for call in calls
        ],
        "withheld": [
            (call.get("identity") or {}).get("withheld", []) for call in calls
        ],
        "cause": [(call.get("signal") or {}).get("cause") for call in calls],
    }


def evidence(trace: Path, root: Path) -> dict:
    """`dna evidence` of one trace, summarised; `{"error": ...}` when it fails."""
    root = root.resolve()
    with tempfile.TemporaryDirectory(dir=root) as scratch:
        run_dir = Path(scratch)
        os.symlink(root / "config", run_dir / "config")
        result = subprocess.run(
            [binary(), "evidence", str(trace)],
            cwd=run_dir,
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            return {"error": result.stderr.strip()[-400:]}
        out = run_dir / "results" / f"{Path(trace).stem}.evidence.json"
        document = json.loads(out.read_text())
        if "read" not in document:  # a rejected read: no per-read evidence
            return {"error": "no read evidence: " + ", ".join(sorted(document))[:200]}
        return summarise(document)


def evidence_many(traces: list[Path], root: Path, workers: int = 32) -> list[dict]:
    workdir(root)
    with ThreadPoolExecutor(workers) as pool:
        return list(pool.map(lambda t: evidence(t, root), traces))
