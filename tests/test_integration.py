"""
Integration smoke test: run MultiQC against example data.

This test actually invokes multiqc as a subprocess against the example
CSV files to verify the full plugin lifecycle works end-to-end.

Marked with @pytest.mark.integration — skip with: pytest -m "not integration"
"""

import subprocess
import tempfile
from pathlib import Path

import pytest

from tests.conftest import EXAMPLES_DIR


@pytest.mark.integration
class TestMultiqcIntegration:
    def test_multiqc_runs_with_example_data(self):
        """multiqc should complete successfully against the examples directory."""
        with tempfile.TemporaryDirectory() as outdir:
            result = subprocess.run(
                [
                    "multiqc",
                    str(EXAMPLES_DIR),
                    "--outdir", outdir,
                    "--no-data-dir",
                    "--force",
                ],
                capture_output=True,
                text=True,
                timeout=120,
            )
            assert result.returncode == 0, (
                f"multiqc failed with return code {result.returncode}\n"
                f"stdout: {result.stdout}\n"
                f"stderr: {result.stderr}"
            )
            # Check that the report was generated
            html_files = list(Path(outdir).glob("*.html"))
            assert len(html_files) >= 1, f"No HTML report found in {outdir}"

    @pytest.mark.integration
    def test_fusilli_modules_appear_in_report(self):
        """At least one FUSILLI section should appear in the generated report."""
        with tempfile.TemporaryDirectory() as outdir:
            result = subprocess.run(
                [
                    "multiqc",
                    str(EXAMPLES_DIR),
                    "--outdir", outdir,
                    "--no-data-dir",
                    "--force",
                ],
                capture_output=True,
                text=True,
                timeout=120,
            )
            if result.returncode != 0:
                pytest.skip(f"multiqc failed: {result.stderr[:200]}")

            html_files = list(Path(outdir).glob("*.html"))
            assert html_files, "No HTML report generated"

            report_text = html_files[0].read_text()
            fusilli_markers = ["fusilli_detection", "fusilli_diversity",
                               "fusilli_preprocessing", "fusilli_partners"]
            found = [m for m in fusilli_markers if m in report_text]
            assert found, (
                "No FUSILLI module anchors found in report. "
                "Plugin may not be installed or modules failed silently."
            )
