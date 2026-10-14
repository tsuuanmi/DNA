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

    def test_accepts_the_allowed_crate_graph(self) -> None:
        result = self.run_map(
            {
                "Cargo.toml": '[dependencies]\ndna-core = "0"\ndna-sanger = "0"\n',
                "src/lib.rs": "pub use dna_kernel::error;\n",
                "src/pipeline/mod.rs": (
                    "use dna_core::read_call::call_read;\n"
                    "use dna_sanger::read_processing::process;\n"
                    "use crate::{report::build, model::result::Result};\n"
                ),
                "src/report.rs": "use crate::model::result::Result;\n",
                "crates/dna-kernel/Cargo.toml": "[dependencies]\n",
                "crates/dna-kernel/src/error.rs": "pub struct Error;\n",
                "crates/dna-core/Cargo.toml": (
                    '[dependencies]\ndna-kernel = "0"\n'
                    '[dev-dependencies]\ndna-kernel = { version = "0" }\n'
                ),
                "crates/dna-core/src/read_call.rs": (
                    "use dna_kernel::{error::Error, model::reference::Reference};\n"
                    "use crate::model::alignment::Alignment;\n"
                ),
                "crates/dna-core/src/model/alignment.rs": "use crate::model::variant::V;\n",
            }
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_rejects_a_plugin_depending_on_another_plugin(self) -> None:
        result = self.run_map(
            {
                "crates/dna-core/Cargo.toml": '[dependencies]\ndna-sanger = "0"\n',
                "crates/dna-core/src/alignment/mod.rs": (
                    "use dna_sanger::model::basecalls::BaseCalls;\n"
                ),
                "crates/dna-post/src/conformance.rs": "use dna_core::read_call::call_read;\n",
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("core must not depend on dna-sanger", result.stderr)
        self.assertIn(
            "core module alignment must not depend on sanger module model::basecalls",
            result.stderr,
        )
        self.assertIn(
            "post module conformance must not depend on core module read_call",
            result.stderr,
        )

    def test_rejects_the_kernel_depending_on_any_crate(self) -> None:
        result = self.run_map(
            {
                "crates/dna-kernel/src/profile/mod.rs": "use dna_core::model::alignment::A;\n",
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "kernel module profile must not depend on core module model::alignment",
            result.stderr,
        )

    def test_rejects_a_cycle_within_a_crate(self) -> None:
        result = self.run_map(
            {
                "crates/dna-sanger/src/callability/mod.rs": (
                    "use crate::read_processing::SangerConfig;\n"
                ),
                "crates/dna-sanger/src/read_processing/mod.rs": (
                    "use crate::callability::analyze;\n"
                ),
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "sanger: module dependency cycle: callability -> read_processing -> callability",
            result.stderr,
        )

    def test_ignores_tests_and_comments(self) -> None:
        result = self.run_map(
            {
                "crates/dna-core/src/model/variant.rs": (
                    "//! See [`process`](dna_sanger::read_processing::process).\n"
                    "pub struct V; // not dna_sanger::x\n"
                    "#[cfg(test)]\nmod tests {\n    use dna_sanger::x;\n}\n"
                ),
                "crates/dna-core/src/sample/tests.rs": "use dna_sanger::x;\n",
            }
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_accepts_the_harness_on_the_facade_and_plugin_crates(self) -> None:
        result = self.run_map(
            {
                "crates/dna-eval/Cargo.toml": (
                    '[dependencies]\ndna = "0"\ndna-kernel = "0"\ndna-core = "0"\n'
                ),
                "crates/dna-eval/src/main.rs": "mod runner;\n",
                "crates/dna-eval/src/runner.rs": (
                    "use dna::sample_call;\n"
                    "use dna_kernel::checksum::hex_sha256;\n"
                    "use dna_core::consensus::build;\n"
                    "use crate::compare::Failure;\n"
                ),
            }
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_rejects_any_crate_depending_on_the_harness(self) -> None:
        result = self.run_map(
            {
                "Cargo.toml": '[dev-dependencies]\ndna-eval = "0"\n',
                "src/pipeline/mod.rs": "use dna_eval::compare::compare;\n",
                "crates/dna-core/Cargo.toml": '[dependencies]\ndna-eval = "0"\n',
                "crates/dna-eval/Cargo.toml": '[dependencies]\ndna = "0"\n',
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("dna must not depend on dna-eval", result.stderr)
        self.assertIn("core must not depend on dna-eval", result.stderr)
        self.assertIn(
            "dna module pipeline must not depend on eval module compare", result.stderr
        )

    def test_rejects_a_plugin_depending_on_the_facade(self) -> None:
        result = self.run_map(
            {
                "crates/dna-post/Cargo.toml": '[dependencies]\ndna = "0"\n',
                "crates/dna-post/src/conformance.rs": "use dna::sample_call;\n",
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("post must not depend on dna\n", result.stderr)
        self.assertIn(
            "post module conformance must not depend on dna module sample_call",
            result.stderr,
        )

    def test_rejects_an_unknown_crate(self) -> None:
        result = self.run_map({"crates/dna-ngs/src/lib.rs": "pub fn read() {}\n"})
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("crates/dna-ngs: unknown crate", result.stderr)


if __name__ == "__main__":
    unittest.main()
