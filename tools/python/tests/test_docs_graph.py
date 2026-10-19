from __future__ import annotations

import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import ModuleType

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "tools" / "python" / "scripts" / "docs_graph.py"


def load_module() -> ModuleType:
    spec = importlib.util.spec_from_file_location("docs_graph", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["docs_graph"] = module
    spec.loader.exec_module(module)
    return module


docs_graph = load_module()

PROFILE = """---
level: 3
areas: [requirements=requirement, design=design, decisions=decision, governance=governance]
id-prefixes: [decision=ADR]
source-roots: [src]
forbidden-names: [legacy]
ledger: governance/retired.md
reachability-depth: 4
review-days: 365
max-lines: 500
---

# Policy
"""
LEDGER = """# Ledger

| id | title | outcome | superseded-by | retired | full-text |
|---|---|---|---|---|---|
<!-- ledger:end -->
"""
BASE = {
    "docs/README.md": (
        "# Docs\n\n[policy](governance/documentation.md) [ledger](governance/retired.md) "
        "[req](requirements/a.md) [design](design/a.md) [adr](decisions/adr/0001-x.md) "
        "[adr2](decisions/adr/0002-y.md)\n"
    ),
    "docs/governance/documentation.md": PROFILE,
    "docs/governance/retired.md": LEDGER,
    "docs/requirements/a.md": "# A\n\n- **SRS-A-001:** DNA MUST do it.\n",
    "docs/design/a.md": (
        "---\nimplements: [SRS-A]\ndecided-by: [ADR-0001]\ncovers: [src/a.rs]\n---\n\n"
        "# Design\n\nUses ADR-0001 for SRS-A-001.\n"
    ),
    "docs/decisions/adr/0001-x.md": "---\nid: ADR-0001\nstatus: accepted\n---\n\n# ADR-0001: X\n",
    "docs/decisions/adr/0002-y.md": "---\nid: ADR-0002\nstatus: accepted\n---\n\n# ADR-0002: Y\n",
    "src/a.rs": "// SRS-A-001\n",
}


def build(root: Path, files: dict[str, str]) -> None:
    for name, text in {**BASE, **files}.items():
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")


def run(root: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--root", str(root), *args],
        check=False,
        capture_output=True,
        text=True,
    )


class GraphTests(unittest.TestCase):
    def check(self, files: dict[str, str]) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            build(root, files)
            return run(root, "check")

    def test_accepts_a_consistent_tree(self) -> None:
        result = self.check({})
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_rejects_unresolved_wikilink_and_id(self) -> None:
        result = self.check(
            {
                "docs/design/a.md": "---\ndecided-by: [ADR-0001]\n---\n\n[[nowhere]] and ADR-0099.\n"
            }
        )
        self.assertEqual(result.returncode, 1)
        self.assertIn("DOC-LNK-002", result.stdout)
        self.assertIn("DOC-LNK-003", result.stdout)

    def test_prose_citation_must_be_declared_and_sync_declares_it(self) -> None:
        files = {"docs/design/a.md": "# Design\n\nSee ADR-0002 and SRS-A-001.\n"}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            build(root, files)
            first = run(root, "check")
            self.assertEqual(first.returncode, 1)
            self.assertIn("DOC-LNK-006", first.stdout)
            self.assertEqual(run(root, "sync").returncode, 0)
            fixed = run(root, "check")
            self.assertEqual(fixed.returncode, 0, fixed.stdout)
            text = (root / "docs/design/a.md").read_text(encoding="utf-8")
            self.assertIn("decided-by: [ADR-0002]", text)
            self.assertIn("implements: [SRS-A]", text)

    def test_supersession_must_be_symmetric(self) -> None:
        result = self.check(
            {
                "docs/decisions/adr/0002-y.md": (
                    "---\nid: ADR-0002\nstatus: accepted\nsupersedes: [ADR-0001]\n---\n\n# ADR-0002: Y\n"
                )
            }
        )
        self.assertEqual(result.returncode, 1)
        self.assertIn("DOC-GRF-001", result.stdout)

    def test_stale_generated_block_is_reported_then_synced(self) -> None:
        index = "# Decisions\n\n<!-- docs:begin index:decision -->\nold\n<!-- docs:end index:decision -->\n"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            build(root, {"docs/decisions/adr/README.md": index})
            stale = run(root, "sync", "--check")
            self.assertEqual(stale.returncode, 1)
            self.assertIn("DOC-GEN-001", stale.stdout)
            run(root, "sync")
            self.assertEqual(run(root, "sync", "--check").returncode, 0)
            self.assertIn(
                "| [ADR-0001](0001-x.md) | X | Accepted |",
                (root / "docs/decisions/adr/README.md").read_text(encoding="utf-8"),
            )

    def test_retired_id_resolves_through_the_ledger_but_not_while_live(self) -> None:
        ledger = LEDGER.replace(
            "<!-- ledger:end -->",
            "| ADR-0003 | Gone | replaced | ADR-0001 | 2026-10-19 | `abc1234:docs/x.md` |\n<!-- ledger:end -->",
        )
        ok = self.check(
            {
                "docs/governance/retired.md": ledger,
                "docs/design/a.md": (
                    "---\ndecided-by: [ADR-0001, ADR-0003]\n---\n\nSee [[ADR-0003]].\n"
                ),
            }
        )
        self.assertEqual(ok.returncode, 0, ok.stdout)
        live = self.check(
            {
                "docs/governance/retired.md": ledger.replace("ADR-0003", "ADR-0001"),
            }
        )
        self.assertIn("DOC-RET-001", live.stdout)

    def test_trace_of_a_source_path_names_the_design_and_requirement(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            build(root, {})
            result = run(root, "trace", "src/a.rs")
            self.assertIn("docs/design/a.md", result.stdout)
            self.assertIn("docs/requirements/a.md", result.stdout)
            self.assertIn("docs/decisions/adr/0001-x.md", result.stdout)

    def test_impact_requires_the_covering_document_or_a_waiver(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            build(root, {})
            for command in (
                ["git", "init", "-q"],
                ["git", "add", "-A"],
                [
                    "git",
                    "-c",
                    "user.email=a@b",
                    "-c",
                    "user.name=n",
                    "commit",
                    "-qm",
                    "base",
                ],
            ):
                subprocess.run(command, cwd=root, check=True, capture_output=True)
            (root / "src/a.rs").write_text("// changed\n", encoding="utf-8")
            failing = run(root, "impact", "HEAD")
            self.assertEqual(failing.returncode, 1)
            self.assertIn("DOC-IMP-001", failing.stdout)
            waived = run(
                root,
                "impact",
                "HEAD",
                "--waivers",
                "docs-unchanged: design/a — comment only",
            )
            self.assertEqual(waived.returncode, 0, waived.stdout)


class FrontMatterTests(unittest.TestCase):
    def test_quoted_list_entries_may_contain_commas(self) -> None:
        text = '---\nsupersedes: ["ADR-0060 §1, §6", ADR-0058]\n---\nbody\n'
        fields, errors, body, _offset = docs_graph.split_front_matter(text)
        self.assertEqual(errors, [])
        self.assertEqual(fields["supersedes"], ["ADR-0060 §1, §6", "ADR-0058"])
        self.assertEqual(body, "body\n")

    def test_render_round_trips(self) -> None:
        fields = {"id": ["ADR-0001"], "status": ["accepted"], "related": ["PROP-0001"]}
        parsed, errors, _body, _offset = docs_graph.split_front_matter(
            docs_graph.render_front_matter(fields) + "x"
        )
        self.assertEqual(errors, [])
        self.assertEqual(parsed, fields)


if __name__ == "__main__":
    unittest.main()
