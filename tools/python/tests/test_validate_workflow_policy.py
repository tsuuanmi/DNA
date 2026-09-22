from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from scripts import validate_workflow_policy as policy


class WorkflowPolicyTests(unittest.TestCase):
    def validate(self, text: str) -> list[str]:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "ci.yml").write_text(text, encoding="utf-8")
            return policy.violations(root)

    def test_accepts_full_sha_action_pin(self) -> None:
        self.assertEqual(
            self.validate(
                "steps:\n"
                "  - uses: actions/checkout@"
                "3d3c42e5aac5ba805825da76410c181273ba90b1 # v7.0.1\n"
            ),
            [],
        )

    def test_rejects_mutable_action_tag(self) -> None:
        errors = self.validate("steps:\n  - uses: actions/checkout@v7\n")
        self.assertEqual(len(errors), 1)
        self.assertIn("full 40-character commit SHA", errors[0])

    def test_rejects_privileged_untrusted_trigger(self) -> None:
        errors = self.validate("on:\n  pull_request_target:\n")
        self.assertEqual(len(errors), 1)
        self.assertIn("pull_request_target", errors[0])

    def test_rejects_write_all_permissions(self) -> None:
        errors = self.validate("permissions: write-all\n")
        self.assertEqual(len(errors), 1)
        self.assertIn("write-all", errors[0])

    def test_accepts_complete_ci_aggregate(self) -> None:
        workflow = (
            "jobs:\n"
            "  rust-quality:\n"
            "    runs-on: ubuntu-latest\n"
            "  dependency-policy:\n"
            "    runs-on: ubuntu-latest\n"
            "  ci-success:\n"
            "    needs:\n"
            "      - rust-quality\n"
            "      - dependency-policy\n"
            "    runs-on: ubuntu-latest\n"
        )
        self.assertEqual(self.validate(workflow), [])

    def test_rejects_incomplete_ci_aggregate(self) -> None:
        workflow = (
            "jobs:\n"
            "  rust-quality:\n"
            "    runs-on: ubuntu-latest\n"
            "  dependency-policy:\n"
            "    runs-on: ubuntu-latest\n"
            "  ci-success:\n"
            "    needs:\n"
            "      - rust-quality\n"
            "    runs-on: ubuntu-latest\n"
        )
        errors = self.validate(workflow)
        self.assertEqual(len(errors), 1)
        self.assertIn("dependency-policy", errors[0])


if __name__ == "__main__":
    unittest.main()
