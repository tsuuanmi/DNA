from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

from scripts import measure_command as measurement


class MeasureCommandTests(unittest.TestCase):
    def test_requires_command(self) -> None:
        with self.assertRaisesRegex(ValueError, "command is required"):
            measurement.normalized_command(["--"])

    @unittest.skipUnless(sys.platform == "linux", "Linux evidence contract")
    def test_records_successful_command_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "performance.json"
            status = measurement.main(
                [
                    "--output",
                    str(output),
                    "--",
                    sys.executable,
                    "-c",
                    "sum(range(1000))",
                ]
            )

            self.assertEqual(status, 0)
            document = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(document["schema_version"], "dna.performance/v1")
            self.assertEqual(document["exit_code"], 0)
            self.assertEqual(document["platform"]["system"], "Linux")
            self.assertGreaterEqual(document["wall_seconds"], 0)
            self.assertGreater(document["peak_rss_kib"], 0)

    @unittest.skipUnless(sys.platform == "linux", "Linux evidence contract")
    def test_preserves_failure_status_and_writes_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "performance.json"
            status = measurement.main(
                [
                    "--output",
                    str(output),
                    "--",
                    sys.executable,
                    "-c",
                    "raise SystemExit(7)",
                ]
            )

            self.assertEqual(status, 7)
            document = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(document["exit_code"], 7)


if __name__ == "__main__":
    unittest.main()
