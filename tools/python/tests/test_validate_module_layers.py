from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import ClassVar

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
                "crates/dna-eval-lib/Cargo.toml": '[dependencies]\ndna-call = "0"\n',
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("align must not depend on dna-call", result.stderr)
        self.assertIn("consensus must not depend on dna-sample", result.stderr)
        self.assertIn("eval-lib must not depend on dna-call", result.stderr)
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
                "crates/dna-eval/Cargo.toml": '[dependencies]\ndna-eval-lib = "0"\n',
                "crates/dna-eval/src/main.rs": "use dna_eval_lib::Identity;\n",
                "crates/dna-eval-lib/Cargo.toml": (
                    '[dependencies]\ndna = "0"\ndna-kernel = "0"\ndna-consensus = "0"\n'
                ),
                "crates/dna-eval-lib/src/lib.rs": "mod runner;\n",
                "crates/dna-eval-lib/src/runner.rs": (
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
                "src/report.rs": "use dna_eval_lib::compare::compare;\n",
                "crates/dna-call/Cargo.toml": '[dependencies]\ndna-eval-lib = "0"\n',
                "crates/dna-eval/Cargo.toml": '[dependencies]\ndna = "0"\n',
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("dna must not depend on dna-eval", result.stderr)
        self.assertIn("call must not depend on dna-eval-lib", result.stderr)
        self.assertIn("eval must not depend on dna\n", result.stderr)
        self.assertIn(
            "dna module pipeline must not depend on eval module compare", result.stderr
        )
        self.assertIn(
            "dna module report must not depend on eval-lib module compare",
            result.stderr,
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

    STAGE_TREE: ClassVar[dict[str, str]] = {
        "Cargo.lock": "",
        "Cargo.toml": '[dependencies]\ndna-kernel = "0"\ndna-call = "0"\n',
        "rust-toolchain.toml": "",
        "src/lib.rs": "mod config;\nmod pipeline;\nmod report;\n",
        "src/config.rs": "pub(crate) struct Config(dna_call::CallConfig);\n",
        "src/pipeline/mod.rs": (
            "mod evidence;\nmod review;\n"
            "pub(crate) use evidence::evidence_document;\n"
            "use crate::report::review;\n"
        ),
        "src/pipeline/evidence.rs": (
            "use crate::report::{self, build};\n"
            "use super::super::config::Config;\n"
            "pub(crate) fn evidence_document() { report::publish(); build(); }\n"
            "fn kernel() { dna_kernel::checksum::hex(); }\n"
        ),
        "src/pipeline/review.rs": "pub(crate) fn review() {}\n",
        "src/report/mod.rs": (
            "mod atomic;\nmod evidence;\nmod review;\n"
            "pub(crate) use atomic::publish;\n"
            "pub(crate) use evidence::{build as build, Built};\n"
            "pub(crate) use review::review;\n"
        ),
        "src/report/atomic.rs": "pub(crate) fn publish() {}\n",
        "src/report/evidence.rs": "pub(crate) fn build() {}\n",
        "src/report/review.rs": "pub(crate) fn review() {}\n",
        "crates/dna-kernel/Cargo.toml": "[dependencies]\n",
        "crates/dna-kernel/src/lib.rs": "pub mod checksum;\n",
        "crates/dna-call/Cargo.toml": '[dependencies]\ndna-kernel = "0"\n',
        "crates/dna-call/src/lib.rs": "pub struct CallConfig;\n",
        "crates/dna-eval/Cargo.toml": '[dependencies]\ndna-eval-lib = "0"\n',
        "crates/dna-eval-lib/Cargo.toml": '[dependencies]\ndna = "0"\n',
    }
    STAGES = (
        "# the read stage\n"
        "[read]\n"
        "Cargo.lock\nCargo.toml\nrust-toolchain.toml\n"
        "crates/dna-kernel\n"
        "src/lib.rs: config, pipeline, report\n"
        "src/config.rs\n"
        "src/pipeline/mod.rs: evidence_document\n"
        "src/pipeline/evidence.rs\n"
        "src/report/mod.rs: publish, build\n"
        "src/report/atomic.rs\n"
        "src/report/evidence.rs\n"
        "!crates/dna-call\n"
    )

    def run_stages(self, stages: str, **files: str) -> subprocess.CompletedProcess[str]:
        return self.run_map(
            {**self.STAGE_TREE, "crates/dna-eval/stages.txt": stages, **files}
        )

    def test_accepts_a_closed_stage_list(self) -> None:
        result = self.run_stages(self.STAGES)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_rejects_a_stage_list_that_misses_a_used_source(self) -> None:
        stages = (
            self.STAGES.replace("src/report/atomic.rs\n", "")
            .replace("src/config.rs\n", "")
            .replace("crates/dna-kernel\n", "")
            .replace("rust-toolchain.toml\n", "")
            .replace("src/report/mod.rs: publish, build\n", "")
        )
        result = self.run_stages(stages)
        self.assertNotEqual(result.returncode, 0)
        for message in (
            "[read] src/pipeline/evidence.rs uses src/report/atomic.rs, which is not listed",
            "[read] src/pipeline/evidence.rs uses src/config.rs, which is not listed",
            "[read] src/pipeline/evidence.rs uses dna-kernel, which is not listed",
            "[read] src/pipeline/evidence.rs uses src/report/mod.rs, which is not listed",
            "[read] src/report/evidence.rs is in src/report/mod.rs, which is not listed",
            "[read] lists no rust-toolchain.toml",
        ):
            self.assertIn(message, result.stderr)

    def test_rejects_unused_exclusions_and_missing_items(self) -> None:
        stages = (
            self.STAGES.replace(
                "src/pipeline/mod.rs: evidence_document\n",
                "src/pipeline/mod.rs: evidence_document, analyze\n",
            )
            + "!crates/dna-sanger\n"
        )
        result = self.run_stages(stages)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("[read] src/pipeline/mod.rs has no item analyze", result.stderr)
        self.assertIn("[read] crates/dna-sanger does not exist", result.stderr)
        self.assertIn(
            "[read] !crates/dna-sanger is used by no listed source", result.stderr
        )

    def test_rejects_a_listed_crate_whose_dependency_is_not_listed(self) -> None:
        stages = self.STAGES.replace("!crates/dna-call\n", "crates/dna-call\n").replace(
            "crates/dna-kernel\n", ""
        )
        result = self.run_stages(stages)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "[read] crates/dna-call depends on dna-kernel, which is not listed",
            result.stderr,
        )

    def test_rejects_an_unknown_crate(self) -> None:
        result = self.run_map({"crates/dna-ngs/src/lib.rs": "pub fn read() {}\n"})
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("crates/dna-ngs: unknown crate", result.stderr)


if __name__ == "__main__":
    unittest.main()
