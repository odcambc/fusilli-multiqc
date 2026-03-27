"""Tests for fusilli_multiqc.modules.detection."""

from pathlib import Path
from unittest.mock import patch, MagicMock

import pandas as pd
import pytest

from tests.conftest import EXAMPLES_DIR, make_find_log_files


class TestDetectionModule:
    """Test the detection module using example data."""

    def _make_module(self, file_map=None):
        """Instantiate detection module with mocked MultiQC internals."""
        if file_map is None:
            file_map = {
                "fusilli_detection/fusion_qc_metrics": str(EXAMPLES_DIR / "fusion_qc_metrics.csv"),
                "fusilli_detection/sensitivity_metrics": str(EXAMPLES_DIR / "sensitivity_metrics.csv"),
                "fusilli_detection/fusion_counts": str(EXAMPLES_DIR / "fusion_counts_summary.csv"),
                "fusilli_preprocessing/decay_metrics": str(EXAMPLES_DIR / "decay_metrics.csv"),
            }

        from fusilli_multiqc.modules.detection import MultiqcModule

        with patch("multiqc.base_module.BaseMultiqcModule.__init__", return_value=None):
            with patch.object(MultiqcModule, "find_log_files", side_effect=make_find_log_files(file_map)):
                with patch.object(MultiqcModule, "add_data_source"):
                    with patch.object(MultiqcModule, "add_section") as mock_section:
                        with patch.object(MultiqcModule, "general_stats_addcols") as mock_stats:
                            module = MultiqcModule()

        return module, mock_section, mock_stats

    def test_loads_data(self):
        module, _, _ = self._make_module()
        assert module.fusion_qc_data is not None
        assert module.sensitivity_data is not None

    def test_detection_efficiency_calculated(self):
        module, _, _ = self._make_module()
        df = module.fusion_qc_data
        assert "detection_efficiency" in df.columns
        assert (df["detection_efficiency"] >= 0).all()
        assert (df["detection_efficiency"] <= 1).all()

    def test_coverage_metrics_calculated(self):
        module, _, _ = self._make_module()
        df = module.fusion_qc_data
        assert "variant_coverage" in df.columns
        assert "partner_coverage" in df.columns

    def test_partner_coverage_uses_max_detected_not_100(self):
        """Regression: partner_coverage should not divide by hardcoded 100.

        When the CSV doesn't pre-compute partner_coverage, the module should
        use max(unique_partners_detected) as the divisor, so the max sample = 1.0.
        """
        # Build a minimal DataFrame WITHOUT pre-computed partner_coverage
        df = pd.DataFrame({
            "sample": ["A", "B"],
            "unique_partners_detected": [80, 50],
            "reads_processed": [1000, 1000],
            "matched_reads": [100, 50],
            "partner_end_reads": [200, 100],
            "unique_fusions": [10, 5],
            "expected_fusions": [100, 100],
        })
        csv_path = str(EXAMPLES_DIR.parent / "tests" / "_tmp_qc.csv")
        df.to_csv(csv_path, index=False)
        try:
            module, _, _ = self._make_module(file_map={
                "fusilli_detection/fusion_qc_metrics": csv_path,
            })
            result = module.fusion_qc_data
            assert "partner_coverage" in result.columns
            assert result["partner_coverage"].max() == pytest.approx(1.0)
            # Sample B should be 50/80 = 0.625
            assert result.loc[result["sample"] == "B", "partner_coverage"].iloc[0] == pytest.approx(0.625)
        finally:
            Path(csv_path).unlink(missing_ok=True)

    def test_breakpoint_coverage_fallback(self):
        """breakpoint_coverage falls back to variant_coverage when dedicated columns missing."""
        # Build a minimal DataFrame WITHOUT pre-computed breakpoint_coverage
        df = pd.DataFrame({
            "sample": ["A", "B"],
            "reads_processed": [1000, 1000],
            "matched_reads": [100, 50],
            "partner_end_reads": [200, 100],
            "unique_fusions": [10, 5],
            "expected_fusions": [100, 100],
        })
        csv_path = str(EXAMPLES_DIR.parent / "tests" / "_tmp_bp.csv")
        df.to_csv(csv_path, index=False)
        try:
            module, _, _ = self._make_module(file_map={
                "fusilli_detection/fusion_qc_metrics": csv_path,
            })
            result = module.fusion_qc_data
            assert "breakpoint_coverage" in result.columns
            pd.testing.assert_series_equal(
                result["breakpoint_coverage"], result["variant_coverage"], check_names=False
            )
        finally:
            Path(csv_path).unlink(missing_ok=True)

    def test_variant_called_retention(self):
        """variant_called_retention should be calculated when decay data is available."""
        module, _, _ = self._make_module()
        df = module.fusion_qc_data
        assert "variant_called_retention" in df.columns
        assert (df["variant_called_retention"] >= 0).all()
        assert (df["variant_called_retention"] <= 1).all()

    def test_sections_added(self):
        _, mock_section, _ = self._make_module()
        section_names = [call.kwargs.get("name", call.args[0] if call.args else None)
                         for call in mock_section.call_args_list]
        assert "Detection Efficiency" in section_names
        assert "Library Coverage" in section_names

    def test_general_stats_added(self):
        _, _, mock_stats = self._make_module()
        assert mock_stats.called

    def test_raises_user_warning_no_data(self):
        """Module should raise UserWarning when no data files are found."""
        from fusilli_multiqc.modules.detection import MultiqcModule

        with patch("multiqc.base_module.BaseMultiqcModule.__init__", return_value=None):
            with patch.object(MultiqcModule, "find_log_files", return_value=iter([])):
                with pytest.raises(UserWarning):
                    MultiqcModule()
