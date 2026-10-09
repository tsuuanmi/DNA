from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "tools" / "python" / "scripts" / "validate_module_layers.py"


class ModuleLayerTests(unittest.TestCase):
    def run_layers(self, files: dict[str, str]) -> subprocess.CompletedProcess[str]:
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

    def test_accepts_downward_dependencies(self) -> None:
        result = self.run_layers(
            {
                "lib.rs": "mod model;\nmod pipeline;\nuse crate::pipeline::run;\n",
                "model.rs": "pub(crate) struct Call;\n",
                "pipeline/mod.rs": "use crate::model::Call;\nuse crate::{error::Error, model};\n",
                "error.rs": "pub(crate) struct Error;\n",
            }
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_rejects_an_upward_dependency(self) -> None:
        result = self.run_layers(
            {
                "model.rs": "use crate::pipeline::run;\n",
                "pipeline.rs": "pub(crate) fn run() {}\n",
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "layer 0 module model must not depend on layer 4 module pipeline",
            result.stderr,
        )

    def test_rejects_a_same_layer_cycle(self) -> None:
        result = self.run_layers(
            {
                "profile.rs": "use crate::reference::Reference;\n",
                "reference.rs": "use crate::profile::Profile;\n",
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "module dependency cycle: profile -> reference -> profile", result.stderr
        )

    def test_ignores_tests_and_comments(self) -> None:
        result = self.run_layers(
            {
                "model.rs": (
                    "//! See [`run`](crate::pipeline::run).\n"
                    "pub(crate) struct Call; // not crate::pipeline\n"
                    "#[cfg(test)]\nmod tests {\n    use crate::pipeline::run;\n}\n"
                ),
                "model/tests.rs": "use crate::pipeline::run;\n",
                "pipeline.rs": "pub(crate) fn run() {}\n",
            }
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_rejects_an_unmapped_module(self) -> None:
        result = self.run_layers({"plugin.rs": "pub(crate) fn hook() {}\n"})
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("plugin: unmapped module", result.stderr)

    def test_rejects_a_neutral_module_depending_on_sanger(self) -> None:
        result = self.run_layers(
            {
                "lib.rs": "mod alignment;\nmod model;\nmod basecalling;\n",
                "model/mod.rs": "pub(crate) mod basecalls;\npub(crate) mod nucleotide;\n",
                "model/basecalls.rs": "pub(crate) struct BaseCalls;\n",
                "model/nucleotide.rs": "pub(crate) struct Nucleotide;\n",
                "basecalling/mod.rs": "pub(crate) fn call() {}\n",
                "alignment/mod.rs": (
                    "use crate::model::{basecalls::BaseCalls, nucleotide::Nucleotide};\n"
                    "use crate::basecalling::call;\n"
                ),
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("must not depend on model::basecalls", result.stderr)
        self.assertIn("must not depend on Sanger module basecalling", result.stderr)
        self.assertNotIn("model::nucleotide", result.stderr)


if __name__ == "__main__":
    unittest.main()
