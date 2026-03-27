"""Tests for fusilli_multiqc.modules.partners."""

from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from tests.conftest import EXAMPLES_DIR, make_find_log_files


class TestPartnersModule:
    """Test the partners module using example data."""

    def _make_module(self, file_map=None):
        if file_map is None:
            file_map = {
                "fusilli_partners/partner_counts": str(EXAMPLES_DIR / "partner_counts_summary.csv"),
            }

        from fusilli_multiqc.modules.partners import MultiqcModule

        with patch("multiqc.base_module.BaseMultiqcModule.__init__", return_value=None):
            with patch.object(MultiqcModule, "find_log_files", side_effect=make_find_log_files(file_map)):
                with patch.object(MultiqcModule, "add_data_source"):
                    with patch.object(MultiqcModule, "add_section") as mock_section:
                        with patch.object(MultiqcModule, "general_stats_addcols"):
                            module = MultiqcModule()

        return module, mock_section

    def test_loads_data(self):
        module, _ = self._make_module()
        assert module.partner_data is not None
        assert "partner_name" in module.partner_data.columns

    def test_processed_data_structure(self):
        module, _ = self._make_module()
        pd_data = module.processed_data
        assert pd_data is not None
        assert "detection_matrix" in pd_data
        assert "count_matrix" in pd_data
        assert "partner_coverage" in pd_data
        assert "end_linker_data" in pd_data
        assert "samples" in pd_data

    def test_detection_matrix_binary(self):
        """Detection matrix values should be 0 or 1."""
        module, _ = self._make_module()
        dm = module.processed_data["detection_matrix"]
        for partner, samples in dm.items():
            for sample, val in samples.items():
                assert val in (0, 1), f"{partner}/{sample} = {val}"

    def test_no_top50_cutoff(self):
        """Regression: all partners should be included, not just top 50."""
        module, _ = self._make_module()
        pc = module.processed_data["partner_coverage"]
        num_partners = len(module.partner_data)
        assert len(pc) == num_partners

    def test_end_linker_data(self):
        module, _ = self._make_module()
        eld = module.processed_data["end_linker_data"]
        assert not eld.empty
        assert "end_counts" in eld.columns
        assert "linker_counts" in eld.columns

    def test_sections_added(self):
        _, mock_section = self._make_module()
        section_names = [call.kwargs.get("name") for call in mock_section.call_args_list]
        assert "Partner Detection Heatmap" in section_names
        assert "Partner Coverage" in section_names
        assert "Partner End vs Linker" in section_names

    def test_raises_user_warning_no_data(self):
        from fusilli_multiqc.modules.partners import MultiqcModule

        with patch("multiqc.base_module.BaseMultiqcModule.__init__", return_value=None):
            with patch.object(MultiqcModule, "find_log_files", return_value=iter([])):
                with pytest.raises(UserWarning):
                    MultiqcModule()
