from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from scripts import validate_configuration_reference as reference

CONFIG = """schema_version = 3
profile = "profiles/target.toml"

[alpha]
ratio = 0.50
count = 4
edges = [-1.0, 1.0]

[beta]
flag = 2
"""

TABLE = """# Configuration

## Schema version 3

| Section | Key | Value | Validation |
|---|---|---:|---|
| root | `schema_version` | `3` | exactly 3 |
| | `profile` | `profiles/target.toml` | a path |
| `alpha` | `ratio` | `0.5` | fraction |
| | `count` | `4` | positive |
| | `edges` | `[-1.0, 1.0]` | ascending |
| `beta` | `flag` | `2` | positive |

Prose after the table.
"""

TYPES = "const SCHEMA_VERSION: u32 = 3;\n"


class ConfigurationReferenceTests(unittest.TestCase):
    def validate(
        self, config: str = CONFIG, table: str = TABLE, types: str = TYPES
    ) -> list[str]:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for path, text in (
                (reference.CONFIG, config),
                (reference.REFERENCE, table),
                (reference.TYPES, types),
            ):
                (root / path).parent.mkdir(parents=True, exist_ok=True)
                (root / path).write_text(text, encoding="utf-8")
            return reference.violations(root)

    def test_accepts_a_matching_table(self) -> None:
        self.assertEqual(self.validate(), [])

    def test_rejects_a_key_missing_from_the_table(self) -> None:
        errors = self.validate(
            table=TABLE.replace("| `beta` | `flag` | `2` | positive |\n", "")
        )
        self.assertTrue(any("beta.flag is not in" in error for error in errors), errors)

    def test_rejects_a_key_missing_from_the_configuration(self) -> None:
        errors = self.validate(config=CONFIG.replace("count = 4\n", ""))
        self.assertTrue(
            any("alpha.count is not in" in error for error in errors), errors
        )

    def test_rejects_a_different_value_or_type(self) -> None:
        errors = self.validate(
            table=TABLE.replace("`0.5`", "`0.6`").replace("`4`", "`4.0`")
        )
        self.assertEqual(len(errors), 2, errors)

    def test_an_empty_section_cell_continues_the_section_above(self) -> None:
        errors = self.validate(table=TABLE.replace("| `beta` | `flag`", "| | `flag`"))
        self.assertTrue(any("alpha.flag" in error for error in errors), errors)

    def test_rejects_differing_schema_versions(self) -> None:
        for errors in (
            self.validate(
                table=TABLE.replace("## Schema version 3", "## Schema version 2")
            ),
            self.validate(types=TYPES.replace("= 3;", "= 4;")),
        ):
            self.assertTrue(
                any("schema versions differ" in error for error in errors), errors
            )


if __name__ == "__main__":
    unittest.main()
