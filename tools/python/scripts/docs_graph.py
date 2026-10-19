"""Documentation knowledge graph: check, sync, trace, impact, candidates, retire.

Implements docs/governance/documentation-architecture.md for this repository. The
profile is the front matter of docs/governance/documentation.md. Structure rules
(DOC-STR-*, DOC-LNK-001) live in validate_docs_structure.py.
"""

from __future__ import annotations

import argparse
import fnmatch
import os
import re
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import UTC, date, datetime, timedelta
from pathlib import Path

PROFILE_DOC = "governance/documentation.md"

SCALAR_KEYS = ("type", "status", "id", "date", "created", "summary", "reviewed")
RELATION_KEYS = (
    "implements",
    "decided-by",
    "covers",
    "verified-by",
    "supersedes",
    "evidence",
    "related",
)
LIST_KEYS = (
    "owners",
    "implementation",
    "family",
    "family-labels",
    "unverified",
    *RELATION_KEYS,
)
PROFILE_KEYS = (
    "level",
    "areas",
    "id-prefixes",
    "source-roots",
    "forbidden-names",
    "ledger",
    "reachability-depth",
    "review-days",
    "max-lines",
)
KEY_ORDER = (
    "id",
    "type",
    "status",
    "date",
    "created",
    "summary",
    "reviewed",
    "owners",
    "family",
    *RELATION_KEYS,
)

LIVING_STATUS = ("draft", "active", "deprecated")
STATUS_BY_TYPE = {
    "decision": (
        "proposed",
        "accepted",
        "rejected",
        "superseded",
        "superseded-in-part",
        "deprecated",
    ),
    "proposal": (
        "draft",
        "proposed",
        "accepted",
        "rejected",
        "implementing",
        "implemented",
        "superseded",
        "withdrawn",
    ),
    "research": ("active", "concluded", "promoted"),
}
NON_LIVING = ("decision", "proposal", "research")
RETIRABLE_STATUS = {
    "decision": ("superseded", "rejected", "deprecated"),
    "proposal": ("implemented", "rejected", "superseded", "withdrawn"),
    "research": ("concluded", "promoted"),
}
GROUPS = {
    "proposal": {
        "open": ("draft", "proposed", "accepted", "implementing"),
        "done": ("implemented", "rejected", "superseded", "withdrawn"),
    },
    "research": {"open": ("active",), "done": ("concluded", "promoted")},
}
STATUS_LABEL = {
    "accepted": "Accepted",
    "proposed": "Proposed",
    "rejected": "Rejected",
    "superseded": "Superseded",
    "superseded-in-part": "Superseded in part",
    "deprecated": "Deprecated",
}

ITEM_DEF = re.compile(r"^\s*[-*]\s+\*\*((?:[A-Z][A-Z0-9]*-)+\d+)[:*]", re.MULTILINE)
WIKILINK = re.compile(r"\[\[([^\]|#]+)(?:#[^\]|]*)?(?:\|[^\]]*)?\]\]")
MD_LINK = re.compile(r"!?\[[^\]]*\]\(([^)\s]+)")
FENCE = re.compile(r"^[ \t]*(```|~~~).*?^[ \t]*\1[ \t]*$", re.MULTILINE | re.DOTALL)
INLINE_CODE = re.compile(r"`[^`\n]+`")
GENERATED = re.compile(
    r"(<!-- docs:begin (\S+) -->\n)(.*?)(<!-- docs:end \2 -->)", re.DOTALL
)
COMMENT = re.compile(r"<!--.*?-->", re.DOTALL)
CODE_PATH = re.compile(r"^(?:[\w.\-]+/)+[\w.\-]*$")
FULL_TEXT = re.compile(r"^[0-9a-f]{7,40}:\S.*$")
WAIVER = re.compile(r"docs-unchanged:\s*(\S+)\s*[—–-]+\s*(\S.*)")


def today() -> date:
    return datetime.now(UTC).date()


@dataclass
class Finding:
    rule: str
    severity: str
    file: str
    line: int
    message: str


@dataclass(eq=False)
class Doc:
    path: Path
    rel: str
    ident: str
    fm: dict[str, list[str]]
    fm_errors: list[str]
    text: str
    body: str
    body_offset: int
    router: bool
    type: str
    status: str
    id: str | None
    area: str
    items: dict[str, int] = field(default_factory=dict)

    def first(self, key: str) -> str:
        values = self.fm.get(key, [])
        return values[0] if values else ""


@dataclass
class LedgerRow:
    id: str
    title: str
    outcome: str
    superseded_by: list[str]
    retired: str
    full_text: str


@dataclass
class Profile:
    level: int
    areas: dict[str, str]
    id_prefixes: dict[str, str]
    source_roots: list[str]
    forbidden_names: list[str]
    ledger: str
    reachability_depth: int
    review_days: int
    max_lines: int


def unquote(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
        return value[1:-1]
    return value


def split_front_matter(text: str) -> tuple[dict[str, list[str]], list[str], str, int]:
    """Return fields, errors, body, and the number of lines before the body."""
    if not text.startswith("---\n"):
        return {}, [], text, 0
    end = text.find("\n---\n", 3)
    if end < 0:
        return {}, ["front matter is not closed"], text, 0
    block = text[4:end]
    body = text[end + 5 :]
    offset = text[: end + 5].count("\n")
    fields: dict[str, list[str]] = {}
    errors: list[str] = []
    key = ""
    for raw in block.split("\n"):
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith("- ") and key:
            fields[key].append(unquote(stripped[2:]))
            continue
        name, sep, value = raw.partition(":")
        if not sep or raw[:1] in " \t":
            errors.append(f"cannot parse front matter line: {raw!r}")
            continue
        key = name.strip()
        value = value.strip()
        if key in fields:
            errors.append(f"duplicate front matter key: {key}")
        if value == "":
            fields[key] = []
        elif value.startswith("["):
            if not value.endswith("]"):
                errors.append(f"unterminated list for {key}")
                fields[key] = []
                continue
            inner = value[1:-1]
            fields[key] = [unquote(part) for part in split_list(inner)]
        else:
            fields[key] = [unquote(value)]
    return fields, errors, body, offset


def split_list(inner: str) -> list[str]:
    """Split an inline list on commas that are not inside quotes."""
    parts: list[str] = []
    current = ""
    quote = ""
    for char in inner:
        if quote:
            quote = "" if char == quote else quote
        elif char in "\"'":
            quote = char
        elif char == ",":
            parts.append(current)
            current = ""
            continue
        current += char
    parts.append(current)
    return [part.strip() for part in parts if part.strip()]


def quote_if_needed(value: str) -> str:
    if (
        re.search(r"(^[\[\]{}&*!|>'\"%@`#-])|(: )|(\s#)|,", value)
        or value != value.strip()
    ):
        return '"' + value.replace('"', "'") + '"'
    return value


def render_front_matter(fields: dict[str, list[str]]) -> str:
    keys = [k for k in KEY_ORDER if k in fields]
    keys += [k for k in fields if k not in KEY_ORDER]
    lines = ["---"]
    for key in keys:
        values = fields[key]
        if key in SCALAR_KEYS:
            lines.append(f"{key}: {quote_if_needed(values[0]) if values else ''}")
            continue
        inline = f"{key}: [{', '.join(quote_if_needed(v) for v in values)}]"
        if len(inline) <= 100:
            lines.append(inline)
        else:
            lines.append(f"{key}:")
            lines.extend(f"  - {quote_if_needed(v)}" for v in values)
    lines.append("---")
    return "\n".join(lines) + "\n"


def mask(text: str, pattern: re.Pattern[str]) -> str:
    """Blank matches of pattern, keeping newlines so line numbers survive."""
    return pattern.sub(lambda m: re.sub(r"[^\n]", " ", m.group(0)), text)


def prose(text: str) -> str:
    """Text a reader treats as prose: no fences, generated blocks, comments, or code."""
    text = mask(text, FENCE)
    text = mask(text, GENERATED)
    text = mask(text, COMMENT)
    return mask(text, INLINE_CODE)


def line_of(text: str, index: int) -> int:
    return text.count("\n", 0, index) + 1


def split_qualified(entry: str) -> tuple[str, str]:
    head, _, rest = entry.strip().partition(" ")
    return head, rest.strip()


def strip_wiki(entry: str) -> str:
    match = re.fullmatch(r"\[\[([^\]|]+)(?:\|[^\]]*)?\]\]", entry.strip())
    return match.group(1).strip() if match else entry.strip()


def number_of(identifier: str) -> int:
    match = re.search(r"-(\d+)$", identifier)
    return int(match.group(1)) if match else -1


def load_profile(root: Path) -> Profile:
    """Read the repository profile from the front matter of the policy document."""
    path = root / "docs" / PROFILE_DOC
    fields, errors, _body, _offset = split_front_matter(path.read_text("utf-8"))
    if errors or "level" not in fields:
        raise SystemExit(f"docs/{PROFILE_DOC}: repository profile front matter missing")

    def pairs(key: str) -> dict[str, str]:
        return dict(p.split("=", 1) for p in fields.get(key, []))

    def number(key: str, default: int) -> int:
        return int(fields[key][0]) if fields.get(key) else default

    return Profile(
        level=number("level", 1),
        areas=pairs("areas"),
        id_prefixes=pairs("id-prefixes"),
        source_roots=fields.get("source-roots", []),
        forbidden_names=fields.get("forbidden-names", []),
        ledger=fields.get("ledger", ["governance/retired.md"])[0],
        reachability_depth=number("reachability-depth", 3),
        review_days=number("review-days", 365),
        max_lines=number("max-lines", 500),
    )


class Graph:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.docs_dir = root / "docs"
        self.findings: list[Finding] = []
        self.profile = self.load_profile()
        self.docs: list[Doc] = []
        self.by_ident: dict[str, Doc] = {}
        self.by_id: dict[str, Doc] = {}
        self.by_path: dict[Path, Doc] = {}
        self.items: dict[str, Doc] = {}
        self.namespaces: dict[str, list[Doc]] = defaultdict(list)
        self.ledger_doc: Doc | None = None
        self._tracked: list[str] | None = None
        self._sources: set[str] | None = None
        self.ledger: dict[str, LedgerRow] = {}
        self.successors: dict[str, list[tuple[Doc, str]]] = defaultdict(list)
        self.load_docs()
        self.index()

    def add(
        self, rule: str, severity: str, doc: Doc | str, line: int, message: str
    ) -> None:
        name = doc if isinstance(doc, str) else f"docs/{doc.rel}"
        self.findings.append(Finding(rule, severity, name, line, message))

    # ------------------------------------------------------------------ loading

    def load_profile(self) -> Profile:
        return load_profile(self.root)

    def ident_of(self, rel: str) -> str:
        ident = rel[:-3]
        if ident == "README":
            return ""
        return ident.removesuffix("/README")

    def load_docs(self) -> None:
        for path in sorted(self.docs_dir.rglob("*.md")):
            rel = path.relative_to(self.docs_dir).as_posix()
            text = path.read_text("utf-8")
            fields, errors, body, offset = split_front_matter(text)
            first = rel.split("/", 1)[0]
            area_type = self.profile.areas.get(first, "")
            router = path.name == "README.md" and not fields
            doc_type = fields.get("type", [area_type or "router"])[0]
            if router:
                doc_type = "router"
            status = fields.get("status", ["active"])[0]
            doc = Doc(
                path=path,
                rel=rel,
                ident=self.ident_of(rel),
                fm=fields,
                fm_errors=errors,
                text=text,
                body=body,
                body_offset=offset,
                router=router,
                type=doc_type,
                status=status,
                id=fields["id"][0] if fields.get("id") else None,
                area=first if "/" in rel else "",
            )
            for match in ITEM_DEF.finditer(prose(body)):
                doc.items.setdefault(
                    match.group(1), line_of(body, match.start()) + offset
                )
            self.docs.append(doc)

    def index(self) -> None:
        for doc in self.docs:
            self.by_path[doc.path] = doc
            self.by_ident[doc.ident] = doc
            if doc.id:
                if doc.id in self.by_id:
                    self.add("DOC-ID-002", "error", doc, 1, f"duplicate id {doc.id}")
                self.by_id[doc.id] = doc
            for item in doc.items:
                if item in self.items and self.items[item] is not doc:
                    self.add(
                        "DOC-ID-002",
                        "error",
                        doc,
                        doc.items[item],
                        f"{item} is also defined in docs/{self.items[item].rel}",
                    )
                self.items.setdefault(item, doc)
                namespace = re.sub(r"-\d+$", "", item)
                if doc not in self.namespaces[namespace]:
                    self.namespaces[namespace].append(doc)
            for entry in doc.fm.get("supersedes", []):
                target, qualifier = split_qualified(strip_wiki(entry))
                self.successors[target].append((doc, qualifier))
        self.load_ledger()

    def load_ledger(self) -> None:
        self.ledger_doc = self.by_ident.get(self.profile.ledger[:-3])
        if self.ledger_doc is None:
            return
        for number, line in enumerate(self.ledger_doc.body.split("\n")):
            if not line.startswith("|"):
                continue
            cells = [c.strip().strip("`") for c in line.strip().strip("|").split("|")]
            if len(cells) != 6 or cells[0] in ("id", "") or set(cells[0]) <= set("-: "):
                continue
            row = LedgerRow(
                id=cells[0],
                title=cells[1],
                outcome=cells[2],
                superseded_by=[
                    x.strip() for x in cells[3].split(",") if x.strip() not in ("", "-")
                ],
                retired=cells[4],
                full_text=cells[5],
            )
            if row.id in self.ledger:
                self.add(
                    "DOC-RET-002",
                    "error",
                    self.ledger_doc,
                    number + 1 + self.ledger_doc.body_offset,
                    f"duplicate ledger row {row.id}",
                )
            self.ledger[row.id] = row

    # ---------------------------------------------------------------- resolving

    def prefixes(self, kind: str) -> str:
        return self.profile.id_prefixes.get(kind, "")

    def resolve(self, target: str) -> Doc | str | None:
        """A document, a marker string for a ledger/namespace entry, or None."""
        target = strip_wiki(target)
        if target in self.by_id:
            return self.by_id[target]
        if target in self.items:
            return self.items[target]
        if target in self.ledger:
            return "ledger"
        if target in self.namespaces:
            return self.namespaces[target][0]
        ident = target.removesuffix(".md").strip("/")
        ident = ident.removesuffix("/README")
        if ident in self.by_ident:
            return self.by_ident[ident]
        tails = [
            d for d in self.docs if d.ident.rsplit("/", 1)[-1] == ident and d.ident
        ]
        if len(tails) == 1:
            return tails[0]
        return None

    def id_regex(self) -> re.Pattern[str]:
        names = {self.prefixes("decision"), self.prefixes("proposal"), *self.namespaces}
        names.discard("")
        if not names:
            return re.compile(r"(?!x)x")
        alternation = "|".join(
            sorted((re.escape(n) for n in names), key=len, reverse=True)
        )
        return re.compile(rf"\b((?:{alternation})-\d+)\b")

    def covered_paths(self, entry: str) -> list[str]:
        entry = entry.strip()
        if any(ch in entry for ch in "*?["):
            return [p for p in self.tracked_paths() if fnmatch.fnmatchcase(p, entry)]
        target = self.root / entry
        if target.is_dir():
            prefix = entry.rstrip("/") + "/"
            return [p for p in self.tracked_paths() if p.startswith(prefix)] or [entry]
        return [entry] if target.exists() else []

    def tracked_paths(self) -> list[str]:
        if self._tracked is None:
            result = subprocess.run(
                ["git", "ls-files"],
                cwd=self.root,
                capture_output=True,
                text=True,
                check=False,
            )
            paths = result.stdout.split("\n") if result.returncode == 0 else []
            if not paths or paths == [""]:
                paths = [
                    p.relative_to(self.root).as_posix()
                    for p in self.root.rglob("*")
                    if p.is_file() and ".git" not in p.parts and "target" not in p.parts
                ]
            self._tracked = [p for p in paths if p]
        return self._tracked

    def matches(self, entry: str, path: str) -> bool:
        entry = entry.strip()
        if any(ch in entry for ch in "*?["):
            return fnmatch.fnmatchcase(path, entry)
        base = entry.rstrip("/")
        return path == base or path.startswith(base + "/")

    # ------------------------------------------------------------------- checks

    def living(self, doc: Doc) -> bool:
        """Current-state documents, including standing ones in a historical area."""
        if doc.type == "router":
            return False
        standing = doc.id is None and doc.path.name != "README.md"
        return doc.type not in NON_LIVING or standing

    def check(self) -> list[Finding]:
        for doc in self.docs:
            self.check_identity(doc)
            self.check_edges(doc)
            self.check_citations(doc)
            self.check_freshness(doc)
        self.check_supersession()
        self.check_reachability()
        self.check_ledger()
        self.check_generated()
        self.check_requirements()
        self.check_governance()
        return [f for f in self.findings if self.active(f.rule, f.severity)]

    def active(self, rule: str, severity: str) -> bool:
        level = self.profile.level
        if rule.startswith("DOC-ID") or rule == "DOC-LNK-001":
            return True
        if rule.startswith(("DOC-RET", "DOC-IMP", "DOC-FRS")) or rule in (
            "DOC-GRF-004",
            "DOC-GRF-005",
        ):
            return level >= 3
        return level >= 2

    def check_identity(self, doc: Doc) -> None:
        for message in doc.fm_errors:
            self.add("DOC-ID-003", "error", doc, 1, message)
        known = set(SCALAR_KEYS) | set(LIST_KEYS)
        if doc.rel == PROFILE_DOC:
            known |= set(PROFILE_KEYS)
        for key in doc.fm:
            if key not in known:
                self.add(
                    "DOC-LNK-005",
                    "error",
                    doc,
                    1,
                    f"unknown front matter field {key!r}",
                )
        if doc.router:
            return
        types = set(self.profile.areas.values())
        if doc.type not in types:
            self.add(
                "DOC-ID-001", "error", doc, 1, f"type {doc.type!r} is not an area type"
            )
            return
        prefix = self.prefixes(doc.type)
        numbered = re.match(r"\d+-", doc.path.stem) is not None
        typed_id = bool(prefix) and (numbered or doc.id is not None)
        allowed = (
            STATUS_BY_TYPE.get(doc.type, LIVING_STATUS) if typed_id else LIVING_STATUS
        )
        if doc.status not in allowed:
            self.add(
                "DOC-ID-003",
                "error",
                doc,
                1,
                f"status {doc.status!r} is not allowed for {doc.type}",
            )
        if not typed_id:
            return
        if doc.id is None or not doc.fm.get("status"):
            self.add(
                "DOC-ID-003",
                "error",
                doc,
                1,
                f"{doc.type} documents require id and status",
            )
            return
        digits = re.match(r"(\d+)-", doc.path.stem)
        if digits and doc.id != f"{prefix}-{digits.group(1)}":
            self.add(
                "DOC-ID-002",
                "error",
                doc,
                1,
                f"id {doc.id} does not match file name prefix {digits.group(1)}",
            )
        if not doc.id.startswith(prefix + "-"):
            self.add(
                "DOC-ID-002",
                "error",
                doc,
                1,
                f"id {doc.id} does not start with {prefix}-",
            )

    def check_edges(self, doc: Doc) -> None:
        for key in RELATION_KEYS:
            for entry in doc.fm.get(key, []):
                value = strip_wiki(entry)
                if key == "covers":
                    if not self.covered_paths(value):
                        self.add(
                            "DOC-GRF-003",
                            "error",
                            doc,
                            1,
                            f"covers {value!r} matches no path",
                        )
                    continue
                if key == "verified-by":
                    if self.covered_paths(value) or self.resolve(value) is not None:
                        continue
                    self.add(
                        "DOC-LNK-005",
                        "error",
                        doc,
                        1,
                        f"verified-by {value!r} resolves to nothing",
                    )
                    continue
                head = split_qualified(value)[0] if key == "supersedes" else value
                if self.resolve(head) is None:
                    self.add(
                        "DOC-LNK-005",
                        "error",
                        doc,
                        1,
                        f"{key} {head!r} resolves to nothing",
                    )

    def citations(self, doc: Doc, regex: re.Pattern[str]) -> list[tuple[str, int]]:
        text = prose(doc.body)
        return [
            (m.group(1), line_of(text, m.start()) + doc.body_offset)
            for m in regex.finditer(text)
        ]

    def check_citations(self, doc: Doc) -> None:
        if doc is self.ledger_doc:
            return
        text = prose(doc.body)
        for match in WIKILINK.finditer(text):
            target = match.group(1).strip()
            if self.resolve(target) is None:
                line = line_of(text, match.start()) + doc.body_offset
                self.add(
                    "DOC-LNK-002",
                    "error",
                    doc,
                    line,
                    f"wikilink [[{target}]] resolves to nothing",
                )
        regex = self.id_regex()
        for ident, line in self.citations(doc, regex):
            if self.resolve(ident) is None and not (
                doc.type == "proposal" and doc.path.stem == "0000-template"
            ):
                self.add(
                    "DOC-LNK-003", "error", doc, line, f"{ident} resolves to nothing"
                )
        if self.living(doc):
            decisions, items = self.missing_declarations(doc)
            for ident in decisions:
                self.add(
                    "DOC-LNK-006",
                    "error",
                    doc,
                    1,
                    f"cites {ident} but does not declare decided-by (run sync)",
                )
            for ident in items:
                self.add(
                    "DOC-LNK-006",
                    "error",
                    doc,
                    1,
                    f"cites {ident} but does not declare implements (run sync)",
                )
            self.check_code_paths(doc)

    def check_code_paths(self, doc: Doc) -> None:
        text = mask(mask(mask(doc.body, FENCE), GENERATED), COMMENT)
        roots = {r.split("/", 1)[0].split("*", 1)[0] for r in self.profile.source_roots}
        roots |= {"tests", "tools", "config"}
        seen: set[str] = set()
        for match in INLINE_CODE.finditer(text):
            value = match.group(0).strip("`")
            first = value.split("/", 1)[0]
            if (
                first not in roots
                or not CODE_PATH.match(value)
                or "*" in value
                or value in seen
            ):
                continue
            seen.add(value)
            if not (self.root / value).exists():
                line = line_of(text, match.start()) + doc.body_offset
                self.add(
                    "DOC-LNK-004", "warning", doc, line, f"path {value} does not exist"
                )

    def missing_declarations(self, doc: Doc) -> tuple[list[str], list[str]]:
        """Cited decisions not in decided-by; cited requirement items not covered."""
        decisions: set[str] = set()
        items: set[str] = set()
        decision_prefix = self.prefixes("decision")
        for ident, _line in self.citations(doc, self.id_regex()):
            if decision_prefix and ident.startswith(decision_prefix + "-"):
                if ident in self.by_id or ident in self.ledger:
                    decisions.add(ident)
            elif ident in self.items and self.items[ident] is not doc:
                items.add(ident)
        declared = {
            split_qualified(strip_wiki(e))[0]
            for k in RELATION_KEYS
            for e in doc.fm.get(k, [])
        }
        implements = {strip_wiki(e) for e in doc.fm.get("implements", [])}
        missing_decisions = sorted(
            d for d in decisions if d not in declared and d != doc.id
        )
        missing_items = sorted(
            i
            for i in items
            if i not in implements and re.sub(r"-\d+$", "", i) not in implements
        )
        return missing_decisions, missing_items

    def check_freshness(self, doc: Doc) -> None:
        lines = doc.text.count("\n")
        if lines > self.profile.max_lines and self.living(doc):
            self.add(
                "DOC-FRS-002",
                "warning",
                doc,
                1,
                f"{lines} lines exceeds the {self.profile.max_lines}-line budget",
            )
        reviewed = doc.first("reviewed")
        if reviewed and self.living(doc):
            try:
                age = today() - date.fromisoformat(reviewed)
            except ValueError:
                self.add(
                    "DOC-ID-003",
                    "error",
                    doc,
                    1,
                    f"reviewed {reviewed!r} is not a date",
                )
                return
            if age > timedelta(days=self.profile.review_days):
                self.add("DOC-FRS-001", "warning", doc, 1, f"last reviewed {reviewed}")

    def check_supersession(self) -> None:
        for target, entries in self.successors.items():
            doc = self.by_id.get(target)
            if doc is None:
                if target not in self.ledger:
                    for successor, _q in entries:
                        self.add(
                            "DOC-GRF-001",
                            "error",
                            successor,
                            1,
                            f"supersedes unknown {target}",
                        )
                continue
            if doc.status not in ("superseded", "superseded-in-part"):
                self.add(
                    "DOC-GRF-001",
                    "error",
                    doc,
                    1,
                    f"{target} is superseded by {entries[0][0].id} but status is {doc.status}",
                )
        for doc in self.docs:
            if (
                doc.status in ("superseded", "superseded-in-part")
                and doc.id
                and doc.id not in self.successors
            ):
                self.add(
                    "DOC-GRF-001",
                    "error",
                    doc,
                    1,
                    f"status {doc.status} but no record supersedes {doc.id}",
                )

    def links_of(self, doc: Doc) -> set[Doc]:
        found: set[Doc] = set()
        text = mask(mask(doc.body, FENCE), COMMENT)
        for match in MD_LINK.finditer(text):
            target = match.group(1).split("#", 1)[0]
            if not target or target.startswith(("http://", "https://", "mailto:")):
                continue
            resolved = (doc.path.parent / target).resolve()
            if resolved.is_dir():
                resolved = resolved / "README.md"
            if resolved in self.by_path:
                found.add(self.by_path[resolved])
        for match in WIKILINK.finditer(text):
            target = self.resolve(match.group(1))
            if isinstance(target, Doc):
                found.add(target)
        for ident, _line in self.citations(doc, self.id_regex()):
            target = self.resolve(ident)
            if isinstance(target, Doc):
                found.add(target)
        return found

    def check_reachability(self) -> None:
        root = self.by_ident.get("")
        if root is None:
            return
        seen = {root: 0}
        queue = [root]
        for current in queue:
            if seen[current] >= self.profile.reachability_depth:
                continue
            for linked in self.links_of(current):
                if linked not in seen:
                    seen[linked] = seen[current] + 1
                    queue.append(linked)
        for doc in self.docs:
            if doc not in seen:
                self.add(
                    "DOC-GRF-002",
                    "error",
                    doc,
                    1,
                    f"not reachable from docs/README.md within {self.profile.reachability_depth} links",
                )

    def check_ledger(self) -> None:
        prefixes = {p: t for t, p in self.profile.id_prefixes.items()}
        if self.ledger_doc is None:
            self.add(
                "DOC-RET-002",
                "error",
                f"docs/{self.profile.ledger}",
                1,
                "ledger document is missing",
            )
            return
        for row in self.ledger.values():
            if row.id in self.by_id or row.id in self.by_ident:
                self.add(
                    "DOC-RET-001",
                    "error",
                    self.ledger_doc,
                    1,
                    f"{row.id} is retired but still live",
                )
            if not (
                FULL_TEXT.match(row.full_text)
                or row.full_text in ("never-committed", "-")
            ):
                self.add(
                    "DOC-RET-002",
                    "error",
                    self.ledger_doc,
                    1,
                    f"{row.id} full-text must be <commit>:<path>",
                )
            for successor in row.superseded_by:
                if self.resolve(successor) is None:
                    self.add(
                        "DOC-RET-002",
                        "error",
                        self.ledger_doc,
                        1,
                        f"{row.id} superseded-by {successor} resolves to nothing",
                    )
        for prefix in prefixes:
            numbers = {
                number_of(i)
                for i in [*self.by_id, *self.ledger]
                if i.startswith(prefix + "-") and number_of(i) >= 0
            }
            if not numbers:
                continue
            gaps = [n for n in range(1, max(numbers) + 1) if n not in numbers]
            if gaps:
                listed = ", ".join(f"{prefix}-{n:04d}" for n in gaps)
                self.add(
                    "DOC-RET-003",
                    "warning",
                    self.ledger_doc,
                    1,
                    f"numbering gaps not in the ledger: {listed}",
                )

    def check_requirements(self) -> None:
        sources = self.source_text()
        implemented: set[str] = set()
        verified: set[str] = set()
        waived: set[str] = set()
        for doc in self.docs:
            implemented |= {strip_wiki(e) for e in doc.fm.get("implements", [])}
            verified |= {strip_wiki(e) for e in doc.fm.get("verified-by", [])}
            waived |= {strip_wiki(e) for e in doc.fm.get("unverified", [])}
        for doc in self.docs:
            if not doc.items or doc.type != "requirement":
                continue
            unimplemented = [
                i
                for i in doc.items
                if i not in implemented and re.sub(r"-\d+$", "", i) not in implemented
            ]
            unverified = [
                i
                for i in doc.items
                if i not in sources
                and i not in verified
                and i not in waived
                and re.sub(r"-\d+$", "", i) not in waived
            ]
            if unimplemented:
                self.add(
                    "DOC-GRF-004",
                    "warning",
                    doc,
                    1,
                    f"{len(unimplemented)} of {len(doc.items)} items are implemented by no document",
                )
            if unverified:
                self.add(
                    "DOC-GRF-004",
                    "warning",
                    doc,
                    1,
                    f"{len(unverified)} of {len(doc.items)} items are cited by no code, test, or verified-by",
                )

    def check_governance(self) -> None:
        for doc in self.docs:
            if (
                doc.type == "decision"
                and doc.status == "accepted"
                and not self.inverse(doc)["Governs"]
            ):
                self.add(
                    "DOC-GRF-005",
                    "warning",
                    doc,
                    1,
                    f"{doc.id} governs no current document",
                )

    def source_text(self) -> set[str]:
        """Requirement items cited anywhere in code, tests, or tooling."""
        if self._sources is not None:
            return self._sources
        cited: set[str] = set()
        pattern = re.compile(r"\b(?:[A-Z][A-Z0-9]*-)+\d+\b")
        for path in self.tracked_paths():
            if path.startswith("docs/") or not path.endswith(
                (".rs", ".py", ".toml", ".yml")
            ):
                continue
            try:
                text = (self.root / path).read_text("utf-8")
            except (OSError, UnicodeDecodeError):
                continue
            cited |= {m for m in pattern.findall(text) if m in self.items}
        self._sources = cited
        return cited

    # ----------------------------------------------------------------- generate

    def title(self, doc: Doc) -> str:
        match = re.search(r"^# (.+)$", mask(doc.body, FENCE), re.MULTILINE)
        title = match.group(1) if match else doc.ident
        title = re.sub(r"^[A-Z]+-\d+:\s*", "", title)
        return re.sub(r"^(Proposal|Research):\s*", "", title)

    def link(self, source: Doc, target: Doc) -> str:
        return Path(os.path.relpath(target.path, source.path.parent)).as_posix()

    def ref(self, source: Doc, identifier: str) -> str:
        target = self.by_id.get(identifier)
        if target is not None:
            return f"[{identifier}]({self.link(source, target)})"
        if identifier in self.ledger and self.ledger_doc is not None:
            return f"[{identifier}]({self.link(source, self.ledger_doc)})"
        return identifier

    def successor_text(self, source: Doc, doc: Doc) -> str:
        label = STATUS_LABEL.get(doc.status, doc.status)
        if doc.status not in ("superseded", "superseded-in-part"):
            return label
        entries = [
            self.ref(source, s.id or s.ident) + (f" ({q})" if q else "")
            for s, q in self.successors.get(doc.id or "", [])
        ]
        return f"{label} by {', '.join(entries)}" if entries else label

    def block(self, doc: Doc, name: str) -> list[str] | None:
        kind, _, rest = name.partition(":")
        if kind == "index":
            return self.gen_index(doc, rest)
        if kind == "families":
            return self.gen_families(doc)
        if kind == "coverage":
            return self.gen_coverage(doc)
        if name == "lifecycle":
            return self.gen_lifecycle(doc)
        if name == "backlinks":
            return self.gen_backlinks(doc)
        return None

    def gen_index(self, doc: Doc, spec: str) -> list[str]:
        doc_type, _, group = spec.partition(":")
        members = sorted(
            (
                d
                for d in self.docs
                if d.type == doc_type and d.id and d.path.stem[:4] != "0000"
            ),
            key=lambda d: (number_of(d.id or ""), d.id),
        )
        if group:
            members = [
                d
                for d in members
                if d.status in GROUPS.get(doc_type, {}).get(group, ())
            ]
        if doc_type == "decision":
            entries: list[tuple[int, str]] = []
            for d in members:
                entries.append(
                    (
                        number_of(d.id or ""),
                        f"| [{d.id}]({self.link(doc, d)}) | {self.title(d)} | {self.successor_text(doc, d)} |",
                    )
                )
            if self.ledger_doc is not None:
                for row in self.ledger.values():
                    if row.id.startswith(self.prefixes("decision") + "-"):
                        by = (
                            f", superseded by {', '.join(row.superseded_by)}"
                            if row.superseded_by
                            else ""
                        )
                        entries.append(
                            (
                                number_of(row.id),
                                f"| [{row.id}]({self.link(doc, self.ledger_doc)}) | {row.title} | Retired{by} |",
                            )
                        )
            return [
                "| ADR | Decision | Status |",
                "|---|---|---|",
                *(row for _n, row in sorted(entries)),
            ]
        lines = []
        for d in members:
            label = (
                f"{d.id} — {self.title(d)}"
                if d.id and doc_type == "proposal"
                else self.title(d)
            )
            summary = d.first("summary")
            lines.append(
                f"- [{label}]({self.link(doc, d)}) — **{STATUS_LABEL.get(d.status, d.status.capitalize())}**"
                + (f"; {summary}" if summary else "")
            )
        return lines or ["_None._"]

    def gen_families(self, doc: Doc) -> list[str]:
        labels = dict(p.split("=", 1) for p in doc.fm.get("family-labels", []))
        members: dict[str, list[Doc]] = defaultdict(list)
        for d in self.docs:
            for family in d.fm.get("family", []):
                members[family].append(d)
        rows = ["| Concern | Decisions († superseded) | Governs |", "|---|---|---|"]
        for slug, label in labels.items():
            decisions = sorted(
                members.get(slug, []), key=lambda d: number_of(d.id or "")
            )
            if not decisions:
                continue
            ids = {d.id for d in decisions}
            cited = ", ".join(
                f"[{d.id}]({self.link(doc, d)})"
                + ("†" if d.status.startswith("superseded") else "")
                for d in decisions
            )
            governed = sorted(
                {
                    g.ident or "docs"
                    for g in self.docs
                    if self.living(g)
                    and ids & {strip_wiki(e) for e in g.fm.get("decided-by", [])}
                }
            )
            links = [f"[{g}]({self.link(doc, self.by_ident[g])})" for g in governed[:6]]
            more = f" +{len(governed) - 6} more" if len(governed) > 6 else ""
            rows.append(f"| {label} | {cited} | {', '.join(links) or '—'}{more} |")
        return rows

    def gen_lifecycle(self, doc: Doc) -> list[str]:
        lines = [f"- **Status:** {self.successor_text(doc, doc)}"]
        if doc.first("date"):
            lines.append(f"- **Date:** {doc.first('date')}")
        supersedes = [
            self.ref(doc, head) + (f" ({q})" if q else "")
            for head, q in (
                split_qualified(strip_wiki(e)) for e in doc.fm.get("supersedes", [])
            )
        ]
        if supersedes:
            lines.append(f"- **Supersedes:** {', '.join(supersedes)}")
        superseded = {
            split_qualified(strip_wiki(e))[0] for e in doc.fm.get("supersedes", [])
        }
        related = [
            self.ref(doc, strip_wiki(e))
            for e in doc.fm.get("related", [])
            if strip_wiki(e) not in superseded
        ]
        if related:
            lines.append(f"- **Related:** {', '.join(related)}")
        return lines

    def inverse(self, doc: Doc) -> dict[str, list[Doc]]:
        result: dict[str, list[Doc]] = {"Governs": [], "Implemented in": []}
        if doc.id:
            result["Governs"] = [
                d
                for d in self.docs
                if d is not doc
                and self.living(d)
                and doc.id in {strip_wiki(e) for e in d.fm.get("decided-by", [])}
            ]
        if doc.items:
            names = set(doc.items) | {re.sub(r"-\d+$", "", i) for i in doc.items}
            result["Implemented in"] = [
                d
                for d in self.docs
                if d is not doc
                and self.living(d)
                and names & {strip_wiki(e) for e in d.fm.get("implements", [])}
            ]
        return result

    def gen_backlinks(self, doc: Doc) -> list[str]:
        lines = []
        for label, docs in self.inverse(doc).items():
            if docs:
                links = ", ".join(
                    f"[{d.ident or 'docs'}]({self.link(doc, d)})"
                    for d in sorted(docs, key=lambda d: d.ident)
                )
                lines.append(f"- **{label}:** {links}")
        return lines or ["- _Not yet referenced by a current document._"]

    def gen_coverage(self, doc: Doc) -> list[str]:
        cited = self.source_text()
        rows = [
            "| Namespace | Defined in | Items | Implemented in | Decisions | Source covered | Cited in code | Unverified |",
            "|---|---|---|---|---|---|---|---|",
        ]
        for namespace in sorted(self.namespaces):
            definers = self.namespaces[namespace]
            if not any(d.type == "requirement" for d in definers):
                continue
            items = sorted(
                i
                for d in definers
                for i in d.items
                if re.sub(r"-\d+$", "", i) == namespace
            )
            names = {namespace, *items}
            implementers = [
                d
                for d in self.docs
                if self.living(d)
                and d not in definers
                and names & {strip_wiki(e) for e in d.fm.get("implements", [])}
            ]
            decisions = sorted(
                {
                    strip_wiki(e)
                    for d in [*definers, *implementers]
                    for e in d.fm.get("decided-by", [])
                },
                key=number_of,
            )
            covered = sorted(
                {strip_wiki(e) for d in implementers for e in d.fm.get("covers", [])}
            )
            cites = len([i for i in items if i in cited])
            waived = {
                strip_wiki(e) for d in definers for e in d.fm.get("unverified", [])
            }
            unverified = len([i for i in items if i in waived and i not in cited])

            def link(d: Doc) -> str:
                return f"[{d.ident or 'docs'}]({self.link(doc, d)})"

            shown = ", ".join(f"`{c}`" for c in covered[:5]) + (
                f" +{len(covered) - 5}" if len(covered) > 5 else ""
            )
            rows.append(
                f"| `{namespace}` | {', '.join(link(d) for d in definers)} | {len(items)} | "
                f"{', '.join(link(d) for d in implementers) or '—'} | "
                f"{', '.join(self.ref(doc, d) for d in decisions) or '—'} | {shown or '—'} | {cites}/{len(items)} | {unverified} |"
            )
        return rows

    def expected_text(self, doc: Doc) -> str:
        text = doc.text
        fences = [(m.start(), m.end()) for m in FENCE.finditer(text)]

        def replace(match: re.Match[str]) -> str:
            if any(a <= match.start() < b for a, b in fences):
                return match.group(0)
            generated = self.block(doc, match.group(2))
            if generated is None:
                self.add(
                    "DOC-GEN-001",
                    "error",
                    doc,
                    line_of(text, match.start()),
                    f"unknown generator {match.group(2)!r}",
                )
                return match.group(0)
            return match.group(1) + "\n".join(generated) + "\n" + match.group(4)

        return GENERATED.sub(replace, text)

    def check_generated(self) -> None:
        for doc in self.docs:
            if self.expected_text(doc) != doc.text:
                self.add(
                    "DOC-GEN-001",
                    "error",
                    doc,
                    1,
                    "generated block is stale (run sync)",
                )

    # --------------------------------------------------------------------- sync

    def sync(self) -> list[str]:
        """Rewrite declarations and generated blocks; return changed files."""
        changed: list[str] = []
        for doc in self.docs:
            if not self.living(doc):
                continue
            decisions, items = self.missing_declarations(doc)
            if not decisions and not items:
                continue
            fields = {k: list(v) for k, v in doc.fm.items()}
            fields.setdefault("decided-by", [])
            fields["decided-by"] = sorted(
                {*fields["decided-by"], *decisions}, key=lambda i: (i[:3], number_of(i))
            )
            if not fields["decided-by"]:
                del fields["decided-by"]
            if items:
                namespaces = sorted({re.sub(r"-\d+$", "", i) for i in items})
                fields["implements"] = sorted(
                    {*fields.get("implements", []), *namespaces}
                )
            doc.path.write_text(
                render_front_matter(fields)
                + ("\n" if doc.body[:1] != "\n" and not doc.fm else "")
                + doc.body,
                "utf-8",
            )
            changed.append(f"docs/{doc.rel}")
        if changed:
            return changed + self.rewrite_after(changed)
        return self.rewrite_after([])

    def rewrite_after(self, already: list[str]) -> list[str]:
        refreshed = Graph(self.root)
        changed: list[str] = []
        for doc in refreshed.docs:
            expected = refreshed.expected_text(doc)
            if expected != doc.text:
                doc.path.write_text(expected, "utf-8")
                changed.append(f"docs/{doc.rel}")
        return changed

    # -------------------------------------------------------------------- trace

    def seeds(self, target: str) -> tuple[list[Doc], list[str]]:
        notes: list[str] = []
        if (self.root / target).exists() and not target.startswith("docs/"):
            path = target.rstrip("/")
            seeds = [
                d
                for d in self.docs
                if any(self.matches(e, path) for e in d.fm.get("covers", []))
            ]
            directory = (
                (self.root / path)
                if (self.root / path).is_dir()
                else (self.root / path).parent
            )
            while directory != self.root.parent and directory != self.root:
                readme = directory / "README.md"
                if readme.exists():
                    notes.append(
                        f"nearest source README: {readme.relative_to(self.root).as_posix()}"
                    )
                    break
                directory = directory.parent
            return seeds, notes
        target = target.removeprefix("docs/")
        resolved = self.resolve(target)
        if isinstance(resolved, Doc):
            seeds = [resolved]
            if target in self.namespaces:
                seeds = list(self.namespaces[target])
            return seeds, notes
        if resolved == "ledger":
            row = self.ledger[strip_wiki(target)]
            notes.append(
                f"retired {row.retired}: {row.outcome}; full text {row.full_text}"
            )
        return [], notes

    def trace(self, target: str) -> list[str]:
        seeds, notes = self.seeds(target)
        if not seeds and not notes:
            return [f"trace {target}", "  nothing resolves to this target"]
        groups: dict[str, dict[str, tuple[Doc, str]]] = defaultdict(dict)

        def note(doc: Doc, why: str) -> None:
            groups[doc.type].setdefault(doc.rel, (doc, why))

        for seed in seeds:
            note(seed, "seed")
            for key in (
                "implements",
                "decided-by",
                "evidence",
                "related",
                "supersedes",
            ):
                for entry in seed.fm.get(key, []):
                    head = split_qualified(strip_wiki(entry))[0]
                    if head in self.ledger:
                        row = self.ledger[head]
                        notes.append(
                            f"{head} is retired: {row.outcome} ({row.full_text})"
                        )
                        continue
                    for definer in self.namespaces.get(head) or [self.resolve(head)]:
                        if isinstance(definer, Doc):
                            note(definer, f"{key} of {seed.ident or 'docs'}")
            for label, docs in self.inverse(seed).items():
                for doc in docs:
                    note(doc, f"{label.lower()} {seed.id or seed.ident}")
            for successor, _q in self.successors.get(seed.id or "", []):
                note(successor, f"supersedes {seed.id}")
        for doc, _why in list(groups.get("requirement", {}).values()):
            for implementer in self.inverse(doc)["Implemented in"]:
                note(implementer, f"implements {doc.ident}")
        out = [f"trace {target}", *(f"  {n}" for n in notes)]
        for kind in sorted(groups):
            out.append(f"  {kind}:")
            ordered = sorted(groups[kind].values(), key=lambda pair: pair[0].rel)
            out += [f"    docs/{d.rel}  ({why})" for d, why in ordered[:25]]
            if len(ordered) > 25:
                out.append(f"    … +{len(ordered) - 25} more")
        for label, key in (
            ("source (covers)", "covers"),
            ("verified-by", "verified-by"),
        ):
            values = sorted({e for s in seeds for e in s.fm.get(key, [])})
            if values:
                out.append(f"  {label}:")
                out += [f"    {v}" for v in values[:25]]
        return out

    # ------------------------------------------------------------------- impact

    def impact(self, base: str, waivers: str) -> list[Finding]:
        diff = subprocess.run(
            ["git", "diff", "--name-only", base],
            cwd=self.root,
            capture_output=True,
            text=True,
            check=False,
        )
        if diff.returncode != 0:
            raise SystemExit(f"git diff {base} failed: {diff.stderr.strip()}")
        changed = {p for p in diff.stdout.split("\n") if p}
        waived = {
            m.group(1).removesuffix(".md").removeprefix("docs/")
            for m in WAIVER.finditer(waivers)
        }
        found: list[Finding] = []
        for doc in self.docs:
            if (
                f"docs/{doc.rel}" in changed
                or doc.ident in waived
                or f"docs/{doc.rel}" in waived
            ):
                continue
            hits = sorted(
                {
                    p
                    for p in changed
                    if not p.startswith("docs/")
                    and any(self.matches(e, p) for e in doc.fm.get("covers", []))
                }
            )
            if hits:
                found.append(
                    Finding(
                        "DOC-IMP-001",
                        "error",
                        f"docs/{doc.rel}",
                        1,
                        f"covers changed {', '.join(hits[:3])}{' …' if len(hits) > 3 else ''}; update it or add 'docs-unchanged: {doc.ident or 'docs'} — reason'",
                    )
                )
        return found

    # ------------------------------------------------------------ retire / list

    def candidates(self) -> list[str]:
        out = []
        for doc in self.docs:
            if doc.router:
                continue
            retirable = RETIRABLE_STATUS.get(doc.type, ("deprecated",))
            if doc.status in retirable or doc.status == "deprecated":
                out.append(f"{doc.id or doc.ident}\t{doc.status}\tdocs/{doc.rel}")
            elif (
                doc.type == "decision"
                and doc.status == "accepted"
                and not self.inverse(doc)["Governs"]
            ):
                out.append(
                    f"{doc.id}\taccepted, governs no current document\tdocs/{doc.rel}"
                )
        return out

    def retire(self, target: str, outcome: str, superseded_by: list[str]) -> list[str]:
        doc = self.resolve(target)
        if not isinstance(doc, Doc) or doc.router:
            raise SystemExit(f"{target}: not a retirable document")
        identity = doc.id or doc.ident
        allowed = (*RETIRABLE_STATUS.get(doc.type, ()), "deprecated")
        if doc.status not in allowed:
            raise SystemExit(
                f"{identity}: status {doc.status} is not retirable ({', '.join(allowed)})"
            )
        if self.ledger_doc is None:
            raise SystemExit("ledger document is missing")
        head = subprocess.run(
            ["git", "rev-parse", "--short=12", "HEAD"],
            cwd=self.root,
            capture_output=True,
            text=True,
            check=False,
        )
        relative = doc.path.relative_to(self.root).as_posix()
        committed = subprocess.run(
            ["git", "cat-file", "-e", f"HEAD:{relative}"],
            cwd=self.root,
            capture_output=True,
            check=False,
        )
        if head.returncode != 0 or committed.returncode != 0:
            raise SystemExit(
                f"{relative} is not committed in HEAD; commit it before retiring so its full text is recoverable"
            )
        successors = superseded_by or [
            s.id for s, _q in self.successors.get(doc.id or "", []) if s.id
        ]
        # Inbound current-state links go to the ledger; wikilinks and bare IDs keep resolving.
        touched: list[str] = []
        ledger_path = self.ledger_doc.path
        for path in [
            *self.root.glob("README.md"),
            *self.root.glob("AGENTS.md"),
            *self.docs_dir.rglob("*.md"),
            *self.root.glob("crates/**/README.md"),
            *self.root.glob("src/**/README.md"),
        ]:
            if path == doc.path:
                continue
            text = path.read_text("utf-8")

            def swap(match: re.Match[str], path: Path = path) -> str:
                target_path = match.group(2).split("#", 1)[0]
                if (
                    target_path
                    and (path.parent / target_path).resolve() == doc.path.resolve()
                ):
                    return f"{match.group(1)}({Path(os.path.relpath(ledger_path, path.parent)).as_posix()}"
                return match.group(0)

            updated = re.sub(r"(\[[^\]]*\])\(([^)\s]+)", swap, text)
            if updated != text:
                path.write_text(updated, "utf-8")
                touched.append(path.relative_to(self.root).as_posix())
        title = self.title(doc).replace("|", "/")
        row = f"| {identity} | {title} | {outcome.replace('|', '/')} | {', '.join(successors) or '-'} | {today().isoformat()} | `{head.stdout.strip()}:{relative}` |\n"
        text = self.ledger_doc.text
        marker = "<!-- ledger:end -->"
        text = (
            text.replace(marker, row + marker)
            if marker in text
            else text.rstrip("\n") + "\n" + row
        )
        ledger_path.write_text(text, "utf-8")
        doc.path.unlink()
        return [
            f"retired {identity}: {relative}",
            f"links rewritten in: {', '.join(touched) or 'none'}",
            f"record the commit with: git tag docs-retired/{today().isoformat()} {head.stdout.strip()}",
        ]


def print_findings(findings: list[Finding], strict: bool, cap: int = 10) -> int:
    errors = [f for f in findings if f.severity == "error"]
    warnings = [f for f in findings if f.severity == "warning"]
    for finding in sorted(errors, key=lambda f: (f.file, f.line)):
        print(f"{finding.file}:{finding.line}: error {finding.rule} {finding.message}")
    shown: dict[str, int] = defaultdict(int)
    hidden: dict[str, int] = defaultdict(int)
    for finding in sorted(warnings, key=lambda f: (f.rule, f.file, f.line)):
        if shown[finding.rule] < cap:
            shown[finding.rule] += 1
            print(
                f"{finding.file}:{finding.line}: warning {finding.rule} {finding.message}"
            )
        else:
            hidden[finding.rule] += 1
    for rule, count in hidden.items():
        print(f"warning {rule}: {count} more not shown")
    return 1 if errors or (strict and warnings) else 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[3]
    )
    sub = parser.add_subparsers(dest="command", required=True)
    check = sub.add_parser("check", help="run the documentation checks")
    check.add_argument("--strict", action="store_true", help="fail on warnings too")
    sync = sub.add_parser("sync", help="rewrite declarations and generated blocks")
    sync.add_argument("--check", action="store_true", help="only report staleness")
    trace = sub.add_parser(
        "trace", help="print the trace closure of a path, ID, or document"
    )
    trace.add_argument("target")
    impact = sub.add_parser(
        "impact", help="documents whose covered paths changed against a base"
    )
    impact.add_argument("base")
    impact.add_argument("--waivers", default=os.environ.get("DOCS_WAIVERS", ""))
    sub.add_parser("candidates", help="list retirement candidates")
    retire = sub.add_parser("retire", help="move a document into the ledger")
    retire.add_argument("target")
    retire.add_argument("--outcome", required=True)
    retire.add_argument("--superseded-by", action="append", default=[])
    args = parser.parse_args(argv)

    graph = Graph(args.root.resolve())
    if args.command == "check":
        return print_findings(graph.check(), args.strict)
    if args.command == "sync":
        if args.check:
            stale = [
                f for f in graph.check() if f.rule in ("DOC-GEN-001", "DOC-LNK-006")
            ]
            return print_findings(stale, False)
        for name in graph.sync():
            print(f"updated {name}")
        return 0
    if args.command == "trace":
        print("\n".join(graph.trace(args.target)))
        return 0
    if args.command == "impact":
        return print_findings(graph.impact(args.base, args.waivers), False)
    if args.command == "candidates":
        print("\n".join(graph.candidates()) or "none")
        return 0
    print("\n".join(graph.retire(args.target, args.outcome, args.superseded_by)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
