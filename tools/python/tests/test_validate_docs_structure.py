from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from scripts import validate_docs_structure as docs

SCHEMA = "../../../docs/reference/schemas/basecalls-v3.schema.json"


class OutputContractTests(unittest.TestCase):
    def validate(self, files: dict[str, str]) -> list[str]:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            schema = (
                root / "docs" / "reference" / "schemas" / "basecalls-v3.schema.json"
            )
            schema.parent.mkdir(parents=True)
            schema.write_text("{}\n", encoding="utf-8")
            for name, text in files.items():
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text, encoding="utf-8")
            return docs.output_contract_violations(root)

    def test_accepts_a_named_schema_and_none_yet(self) -> None:
        errors = self.validate(
            {
                "crates/dna-sanger/src/README.md": (
                    "# dna-sanger\n\nOutput contract:\n"
                    f"[`dna.basecalls/v3`]({SCHEMA}) (`dna basecall`).\n"
                ),
                "crates/dna-eval/src/README.md": (
                    "# dna-eval\n\nOutput contract: none yet (PROP-0011 increment 3).\n"
                ),
                "crates/dna-sanger/src/abif/README.md": "# ABIF\n",
                "src/README.md": "# Rust Source\n",
            }
        )
        self.assertEqual(errors, [])

    def test_rejects_a_crate_readme_without_output_contract(self) -> None:
        errors = self.validate({"crates/dna-call/src/README.md": "# dna-call\n"})
        self.assertEqual(len(errors), 1)
        self.assertIn("names no output contract", errors[0])
        self.assertIn("crates/dna-call/src/README.md", errors[0])

    def test_rejects_a_missing_schema(self) -> None:
        errors = self.validate(
            {
                "crates/dna-call/src/README.md": (
                    "# dna-call\n\nOutput contract: `dna.analysis/v10`\n"
                    "(analysis-v10.schema.json).\n"
                ),
                "src/report/README.md": (
                    "# Report\n\nOutput contract: review-v9.schema.json.\n"
                ),
            }
        )
        self.assertEqual(len(errors), 2)
        self.assertIn("missing schema", errors[0])
        self.assertIn("analysis-v10.schema.json", errors[0])
        self.assertIn("review-v9.schema.json", errors[1])

    def test_rejects_a_contract_without_schema(self) -> None:
        errors = self.validate(
            {
                "crates/dna-call/src/README.md": (
                    "# dna-call\n\nOutput contract: `dna.analysis/v9`.\n"
                ),
            }
        )
        self.assertEqual(len(errors), 1)
        self.assertIn("names no schema", errors[0])


if __name__ == "__main__":
    unittest.main()
