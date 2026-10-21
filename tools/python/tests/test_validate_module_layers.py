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
                "Cargo.toml": '[dependencies]\ndna-call = "0"\ndna-sanger = "0"\n',
                "src/lib.rs": "pub use dna_kernel::error;\n",
                "src/pipeline/mod.rs": (
                    "use dna_call::read_call::call_read;\n"
                    "use dna_sanger::read_processing::process;\n"
                    "use crate::{report::build, model::result::Result};\n"
                ),
                "src/report.rs": "use crate::model::result::Result;\n",
                "crates/dna-kernel/Cargo.toml": "[dependencies]\n",
                "crates/dna-kernel/src/error.rs": "pub struct Error;\n",
                "crates/dna-call/Cargo.toml": (
                    '[dependencies]\ndna-kernel = "0"\ndna-align = "0"\n'
                    '[dev-dependencies]\ndna-align = { version = "0" }\n'
                ),
                "crates/dna-call/src/read_call.rs": (
                    "use dna_kernel::{error::Error, model::reference::Reference};\n"
                    "use dna_align::model::alignment::Alignment;\n"
                    "use crate::model::variant::V;\n"
                ),
                "crates/dna-call/src/model/called_read.rs": "use crate::model::variant::V;\n",
                "crates/dna-consensus/Cargo.toml": '[dependencies]\ndna-call = "0"\n',
                "crates/dna-consensus/src/lib.rs": "use dna_call::Eligibility;\n",
            }
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_rejects_a_plugin_depending_on_another_plugin(self) -> None:
        result = self.run_map(
            {
                "crates/dna-align/Cargo.toml": '[dependencies]\ndna-sanger = "0"\n',
                "crates/dna-align/src/orient.rs": (
                    "use dna_sanger::model::basecalls::BaseCalls;\n"
                ),
                "crates/dna-empop/src/conformance.rs": "use dna_call::read_call::call_read;\n",
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("align must not depend on dna-sanger", result.stderr)
        self.assertIn(
            "align module orient must not depend on sanger module model::basecalls",
            result.stderr,
        )
        self.assertIn(
            "empop module conformance must not depend on call module read_call",
            result.stderr,
        )

    def test_rejects_a_core_step_depending_on_a_later_one(self) -> None:
        result = self.run_map(
            {
                "crates/dna-align/Cargo.toml": '[dependencies]\ndna-call = "0"\n',
                "crates/dna-align/src/orient.rs": "use dna_call::Eligibility;\n",
                "crates/dna-consensus/Cargo.toml": '[dependencies]\ndna-sample = "0"\n',
                "crates/dna-consensus/src/observe.rs": "use dna_sample::aggregate;\n",
                "crates/dna-eval/Cargo.toml": '[dependencies]\ndna-call = "0"\n',
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("align must not depend on dna-call", result.stderr)
        self.assertIn("consensus must not depend on dna-sample", result.stderr)
        self.assertIn("eval must not depend on dna-call", result.stderr)
        self.assertIn(
            "align module orient must not depend on call module Eligibility",
            result.stderr,
        )

    def test_rejects_the_kernel_depending_on_any_crate(self) -> None:
        result = self.run_map(
            {
                "crates/dna-kernel/src/profile/mod.rs": "use dna_align::model::alignment::A;\n",
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "kernel module profile must not depend on align module model::alignment",
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
                "crates/dna-call/src/model/variant.rs": (
                    "//! See [`process`](dna_sanger::read_processing::process).\n"
                    "pub struct V; // not dna_sanger::x\n"
                    "#[cfg(test)]\nmod tests {\n    use dna_sanger::x;\n}\n"
                ),
                "crates/dna-sample/src/aggregate/tests.rs": "use dna_sanger::x;\n",
            }
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_accepts_the_harness_on_the_facade_and_plugin_crates(self) -> None:
        result = self.run_map(
            {
                "crates/dna-eval/Cargo.toml": (
                    '[dependencies]\ndna = "0"\ndna-kernel = "0"\ndna-consensus = "0"\n'
                ),
                "crates/dna-eval/src/main.rs": "mod runner;\n",
                "crates/dna-eval/src/runner.rs": (
                    "use dna::sample_call;\n"
                    "use dna_kernel::checksum::hex_sha256;\n"
                    "use dna_consensus::build;\n"
                    "use crate::compare::Failure;\n"
                ),
            }
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    # SRS-COMPAT-003
    def test_rejects_any_crate_depending_on_the_harness(self) -> None:
        result = self.run_map(
            {
                "Cargo.toml": '[dev-dependencies]\ndna-eval = "0"\n',
                "src/pipeline/mod.rs": "use dna_eval::compare::compare;\n",
                "crates/dna-call/Cargo.toml": '[dependencies]\ndna-eval = "0"\n',
                "crates/dna-eval/Cargo.toml": '[dependencies]\ndna = "0"\n',
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("dna must not depend on dna-eval", result.stderr)
        self.assertIn("call must not depend on dna-eval", result.stderr)
        self.assertIn(
            "dna module pipeline must not depend on eval module compare", result.stderr
        )

    def test_rejects_a_plugin_depending_on_the_facade(self) -> None:
        result = self.run_map(
            {
                "crates/dna-empop/Cargo.toml": '[dependencies]\ndna = "0"\n',
                "crates/dna-empop/src/conformance.rs": "use dna::sample_call;\n",
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("empop must not depend on dna\n", result.stderr)
        self.assertIn(
            "empop module conformance must not depend on dna module sample_call",
            result.stderr,
        )

    def test_rejects_an_unknown_crate(self) -> None:
        result = self.run_map({"crates/dna-ngs/src/lib.rs": "pub fn read() {}\n"})
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("crates/dna-ngs: unknown crate", result.stderr)

    def test_rejects_a_neutral_facade_module_naming_a_modality(self) -> None:
        result = self.run_map(
            {
                "Cargo.toml": '[dependencies]\ndna-sanger = "0"\n',
                "src/input/consensus.rs": (
                    "let reason = dna_sanger::read_processing::reason(label);\n"
                ),
                "src/model/review_result.rs": ("use dna_trace::tract::Tract;\n"),
                "src/report/review/second.rs": "use dna_sanger::run_signal::Key;\n",
                "src/input/evidence.rs": "use dna_sanger::read_processing::reason;\n",
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "src/input/consensus.rs:1: a modality-neutral facade module must not "
            "name dna-sanger",
            result.stderr,
        )
        self.assertIn("src/model/review_result.rs:1", result.stderr)
        # The review's Sanger part and the Sanger evidence reader may.
        self.assertNotIn("second.rs", result.stderr)
        self.assertNotIn("evidence.rs", result.stderr)


if __name__ == "__main__":
    unittest.main()
