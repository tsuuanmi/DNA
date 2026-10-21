"""rCRS access, read placement by seed votes, and read templates (research only)."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import cache
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
RCRS = REPO / "references" / "rCRS.fasta"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


@cache
def rcrs() -> str:
    """The rCRS sequence, upper case."""
    lines = RCRS.read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">")).upper()


def revcomp(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


@dataclass(frozen=True)
class Placement:
    """Where a read's call 0 lies: strand and 1-based rCRS position it reads."""

    strand: str  # "F" or "R"
    start: int  # 1-based rCRS position of call 0 (the first base read)
    votes: int


def place(calls: str, seed: int = 12, window: int = 300) -> Placement | None:
    """Places `calls` on rCRS by the modal diagonal of exact seed matches."""
    reference = rcrs()
    best: Placement | None = None
    for strand, text in (("F", reference), ("R", revcomp(reference))):
        index: dict[str, list[int]] = {}
        for offset in range(len(text) - seed + 1):
            index.setdefault(text[offset : offset + seed], []).append(offset)
        votes: Counter[int] = Counter()
        head = calls[:window]
        for offset in range(len(head) - seed + 1):
            for hit in index.get(head[offset : offset + seed], ()):
                votes[hit - offset] += 1
        if not votes:
            continue
        diagonal, count = votes.most_common(1)[0]
        if best is None or count > best.votes:
            position = diagonal + 1 if strand == "F" else len(text) - diagonal
            best = Placement(strand, position, count)
    return best


def template_map(
    strand: str, start: int, length: int, edits: dict[int, str] | None = None
) -> list[int]:
    """The rCRS position each base of `template(...)` comes from, in read order."""
    reference = rcrs()
    edits = edits or {}
    out: list[int] = []
    positions = (
        range(start, len(reference) + 1) if strand == "F" else range(start, 0, -1)
    )
    for position in positions:
        out.extend([position] * len(edits.get(position, reference[position - 1])))
        if len(out) >= length:
            break
    return out[:length]


def template(
    strand: str, start: int, length: int, edits: dict[int, str] | None = None
) -> str:
    """The read-order template from 1-based `start` on `strand`, `length` bases.

    `edits` maps a 1-based reference position to the reference-strand text that
    replaces its base (an empty string deletes it, several bases insert after
    it), so a C-tract's length is set by editing one of its positions.
    """
    reference = rcrs()
    edits = edits or {}
    pieces = [
        edits.get(position, base) for position, base in enumerate(reference, start=1)
    ]
    if strand == "F":
        out: list[str] = []
        for position in range(start, len(reference) + 1):
            out.append(pieces[position - 1])
            if sum(map(len, out)) >= length:
                break
        return "".join(out)[:length]
    out = []
    for position in range(start, 0, -1):
        out.append(revcomp(pieces[position - 1]))
        if sum(map(len, out)) >= length:
            break
    return "".join(out)[:length]
