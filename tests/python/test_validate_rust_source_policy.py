from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "validate_rust_source_policy.py"


class RustSourcePolicyTests(unittest.TestCase):
    def run_policy(self, source: str) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "lib.rs").write_text(source, encoding="utf-8")
            return subprocess.run(
                [sys.executable, str(SCRIPT), "--root", str(root)],
                check=False,
                capture_output=True,
                text=True,
            )

    def test_accepts_clean_production_rust(self) -> None:
        result = self.run_policy(
            "#![forbid(unsafe_code)]\n"
            "#![deny(deprecated)]\n"
            "pub(crate) fn score(value: i32) -> i32 { value + 1 }\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_rejects_deprecated_api_declarations(self) -> None:
        result = self.run_policy(
            '#[deprecated(note = "old API")]\npub fn old_api() {}\n'
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("deprecated-api", result.stderr)

    def test_rejects_warning_suppression_for_stale_code(self) -> None:
        result = self.run_policy("#[allow(dead_code)]\nfn unused_path() {}\n")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("lint-suppression", result.stderr)

    def test_rejects_explicit_compatibility_identifiers(self) -> None:
        result = self.run_policy("fn legacy_adapter() {}\n")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("legacy-identifier", result.stderr)

    def test_rejects_compatibility_feature_gates(self) -> None:
        result = self.run_policy('#[cfg(feature = "legacy-output")]\nfn output() {}\n')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("compat-feature", result.stderr)


if __name__ == "__main__":
    unittest.main()
