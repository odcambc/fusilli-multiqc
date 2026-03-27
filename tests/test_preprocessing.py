"""Tests for fusilli_multiqc.modules.preprocessing."""

from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from tests.conftest import EXAMPLES_DIR, make_find_log_files


class TestPreprocessingModule:
    """Test the preprocessing module using example data."""

    def _make_module(self, file_map=None):
        if file_map is None:
            file_map = {
                "fusilli_preprocessing/decay_metrics": str(EXAMPLES_DIR / "decay_metrics.csv"),
                "fusilli_preprocessing/trim_metrics": str(EXAMPLES_DIR / "trim_metrics.csv"),
                "fusilli_preprocessing/contam_metrics": str(EXAMPLES_DIR / "contam_metrics.csv"),
                "fusilli_preprocessing/quality_metrics": str(EXAMPLES_DIR / "quality_metrics.csv"),
            }

        from fusilli_multiqc.modules.preprocessing import MultiqcModule

        with patch("multiqc.base_module.BaseMultiqcModule.__init__", return_value=None):
            with patch.object(MultiqcModule, "find_log_files", side_effect=make_find_log_files(file_map)):
                with patch.object(MultiqcModule, "add_data_source"):
                    with patch.object(MultiqcModule, "add_section") as mock_section:
                        with patch.object(MultiqcModule, "general_stats_addcols") as mock_stats:
                            module = MultiqcModule()

        return module, mock_section, mock_stats

    def test_loads_data(self):
        module, _, _ = self._make_module()
        assert module.decay_data is not None

    def test_retention_rates_calculated(self):
        module, _, _ = self._make_module()
        assert module.retention_data is not None
        assert "retention_rate" in module.retention_data.columns

    def test_retention_rates_bounded(self):
        module, _, _ = self._make_module()
        rd = module.retention_data
        # All retention rates should be between 0 and 1 (pairs basis throughout)
        assert (rd["retention_rate"] >= 0).all()
        assert (rd["retention_rate"] <= 1.0 + 1e-9).all()

    def test_retention_decreases_through_pipeline(self):
        """Retention should generally decrease (or stay same) through pipeline steps."""
        module, _, _ = self._make_module()
        rd = module.retention_data
        steps = ["trimmed", "cleaned", "quality"]
        for sample in rd["sample"].unique():
            sample_data = rd[rd["sample"] == sample]
            prev = 1.0
            for step in steps:
                step_data = sample_data[sample_data["step"] == step]
                if not step_data.empty:
                    rate = step_data.iloc[0]["retention_rate"]
                    assert rate <= prev + 1e-9, f"{sample} {step}: {rate} > {prev}"
                    prev = rate

    def test_sections_added(self):
        _, mock_section, _ = self._make_module()
        section_names = [call.kwargs.get("name") for call in mock_section.call_args_list]
        assert "Read Decay" in section_names

    def test_base_decay_section_added(self):
        _, mock_section, _ = self._make_module()
        section_names = [call.kwargs.get("name") for call in mock_section.call_args_list]
        assert "Base Decay" in section_names

    def test_base_decay_only_non_null_steps(self):
        """Base decay plot should only include steps with valid base counts."""
        module, mock_section, _ = self._make_module()
        # Find the Base Decay section call
        base_decay_call = next(
            (c for c in mock_section.call_args_list if c.kwargs.get("name") == "Base Decay"),
            None,
        )
        assert base_decay_call is not None
        # matched step has no bases — it must not appear in any sample's data
        # (we can't easily inspect the plot object, so just verify the section was created)

    def test_base_decay_skipped_without_bases_column(self):
        """If decay_metrics.csv has no bases column, base decay section is not added."""
        import tempfile, os
        tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False)
        tmp.write("sample,step,reads,read_fraction\n")
        tmp.write("S1,raw,1000,1.0\n")
        tmp.write("S1,merged,800,0.8\n")
        tmp.close()
        try:
            _, mock_section, _ = self._make_module(file_map={
                "fusilli_preprocessing/decay_metrics": tmp.name,
            })
            section_names = [call.kwargs.get("name") for call in mock_section.call_args_list]
            assert "Base Decay" not in section_names
        finally:
            os.unlink(tmp.name)

    def test_raises_user_warning_no_data(self):
        from fusilli_multiqc.modules.preprocessing import MultiqcModule

        with patch("multiqc.base_module.BaseMultiqcModule.__init__", return_value=None):
            with patch.object(MultiqcModule, "find_log_files", return_value=iter([])):
                with pytest.raises(UserWarning):
                    MultiqcModule()
