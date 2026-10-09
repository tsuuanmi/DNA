from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "tools" / "python" / "scripts" / "validate_module_layers.py"


class CrateMapTests(unittest.TestCase):
    def run_map(self, files: dict[str, str]) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, source in files.items():
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(source, encoding="utf-8")
            return subprocess.run(
                [sys.executable, str(SCRIPT), "--root", str(root)],
                check=False,
                capture_output=True,
                text=True,
            )

    def test_accepts_dependencies_on_allowed_crates(self) -> None:
        result = self.run_map(
            {
                "lib.rs": "mod model;\nmod pipeline;\nuse crate::pipeline::run;\n",
                "model/mod.rs": "pub(crate) mod nucleotide;\npub(crate) mod alignment;\n",
                "model/nucleotide.rs": "pub(crate) struct Nucleotide;\n",
                "model/alignment.rs": "use crate::model::nucleotide::Nucleotide;\n",
                "alignment/mod.rs": (
                    "use crate::{error::Error, model::{alignment::A, nucleotide::N}};\n"
                ),
                "error.rs": "pub(crate) struct Error;\n",
                "pipeline/mod.rs": "use crate::alignment::align;\nuse crate::basecalling::call;\n",
                "basecalling/mod.rs": "use crate::error::Error;\n",
            }
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_rejects_a_plugin_depending_on_another_plugin(self) -> None:
        result = self.run_map(
            {
                "alignment/mod.rs": "use crate::model::{basecalls::BaseCalls, nucleotide::N};\n",
                "model/basecalls.rs": "pub(crate) struct BaseCalls;\n",
                "model/nucleotide.rs": "pub(crate) struct N;\n",
                "conformance.rs": "use crate::read_call::call_read;\n",
                "read_call.rs": "pub(crate) fn call_read() {}\n",
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "core module alignment must not depend on sanger module model::basecalls",
            result.stderr,
        )
        self.assertIn(
            "post module conformance must not depend on core module read_call",
            result.stderr,
        )
        self.assertNotIn("model::nucleotide", result.stderr)

    def test_rejects_the_kernel_or_a_plugin_depending_on_the_facade(self) -> None:
        result = self.run_map(
            {
                "profile/mod.rs": "use crate::model::alignment::Orientation;\n",
                "model/alignment.rs": "pub(crate) enum Orientation {}\n",
                "read_processing/mod.rs": "use crate::config::Config;\n",
                "config/mod.rs": "pub(crate) struct Config;\n",
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "kernel module profile must not depend on core module model::alignment",
            result.stderr,
        )
        self.assertIn(
            "sanger module read_processing must not depend on dna module config",
            result.stderr,
        )

    def test_assigns_input_children_individually(self) -> None:
        result = self.run_map(
            {
                "input/sanger/abif/decode.rs": "use crate::input::sequence::load;\n",
                "input/sequence.rs": "use crate::input::sanger::abif::load;\n",
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "sanger module input::sanger::abif must not depend on dna module input",
            result.stderr,
        )
        self.assertNotIn("dna module input must not", result.stderr)

    def test_rejects_a_cycle_within_a_crate(self) -> None:
        result = self.run_map(
            {
                "callability/mod.rs": "use crate::read_processing::SangerConfig;\n",
                "read_processing/mod.rs": "use crate::callability::analyze;\n",
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "module dependency cycle: callability -> read_processing -> callability",
            result.stderr,
        )

    def test_ignores_tests_and_comments(self) -> None:
        result = self.run_map(
            {
                "model/nucleotide.rs": (
                    "//! See [`run`](crate::pipeline::run).\n"
                    "pub(crate) struct N; // not crate::pipeline\n"
                    "#[cfg(test)]\nmod tests {\n    use crate::pipeline::run;\n}\n"
                ),
                "alignment/tests.rs": "use crate::pipeline::run;\n",
                "pipeline.rs": "pub(crate) fn run() {}\n",
            }
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_rejects_unmapped_modules_and_dependencies(self) -> None:
        result = self.run_map(
            {
                "scheduler.rs": "pub(crate) fn hook() {}\n",
                "pipeline.rs": "use crate::model::fresh::Thing;\n",
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("scheduler: unmapped module", result.stderr)
        self.assertIn("unmapped dependency crate::model::fresh::Thing", result.stderr)


if __name__ == "__main__":
    unittest.main()
