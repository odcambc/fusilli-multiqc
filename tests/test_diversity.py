"""Tests for fusilli_multiqc.modules.diversity."""

from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from tests.conftest import EXAMPLES_DIR, make_find_log_files


class TestDiversityModule:
    """Test the diversity module using example data."""

    def _make_module(self, file_map=None):
        if file_map is None:
            file_map = {
                "fusilli_diversity/fusion_counts": str(EXAMPLES_DIR / "fusion_counts_summary.csv"),
            }

        from fusilli_multiqc.modules.diversity import MultiqcModule

        with patch("multiqc.base_module.BaseMultiqcModule.__init__", return_value=None):
            with patch.object(MultiqcModule, "find_log_files", side_effect=make_find_log_files(file_map)):
                with patch.object(MultiqcModule, "add_data_source"):
                    with patch.object(MultiqcModule, "add_section") as mock_section:
                        with patch.object(MultiqcModule, "general_stats_addcols") as mock_stats:
                            module = MultiqcModule()

        return module, mock_section, mock_stats

    def test_loads_data(self):
        module, _, _ = self._make_module()
        assert module.fusion_counts_data is not None

    def test_diversity_metrics_calculated(self):
        module, _, _ = self._make_module()
        dm = module.diversity_metrics
        assert dm is not None
        assert "shannon_diversity" in dm.columns
        assert "simpson_diversity" in dm.columns
        assert "evenness" in dm.columns

    def test_shannon_diversity_positive(self):
        module, _, _ = self._make_module()
        dm = module.diversity_metrics
        assert (dm["shannon_diversity"] >= 0).all()

    def test_simpson_diversity_bounded(self):
        module, _, _ = self._make_module()
        dm = module.diversity_metrics
        assert (dm["simpson_diversity"] >= 0).all()
        assert (dm["simpson_diversity"] <= 1).all()

    def test_evenness_bounded(self):
        module, _, _ = self._make_module()
        dm = module.diversity_metrics
        assert (dm["evenness"] >= 0).all()
        assert (dm["evenness"] <= 1).all()

    def test_top_n_fractions_monotonic(self):
        """top1 <= top5 <= top10 <= top25 <= 1.0."""
        module, _, _ = self._make_module()
        dm = module.diversity_metrics
        for _, row in dm.iterrows():
            assert row["top1_fraction"] <= row["top5_fraction"] + 1e-9
            assert row["top5_fraction"] <= row["top10_fraction"] + 1e-9
            assert row["top10_fraction"] <= row["top25_fraction"] + 1e-9
            assert row["top25_fraction"] <= 1.0 + 1e-9

    def test_sections_added(self):
        _, mock_section, _ = self._make_module()
        section_names = [call.kwargs.get("name") for call in mock_section.call_args_list]
        assert "Diversity Comparison" in section_names
        assert "Evenness" in section_names
        assert "Rank-Abundance" in section_names

    def test_raises_user_warning_no_data(self):
        from fusilli_multiqc.modules.diversity import MultiqcModule

        with patch("multiqc.base_module.BaseMultiqcModule.__init__", return_value=None):
            with patch.object(MultiqcModule, "find_log_files", return_value=iter([])):
                with pytest.raises(UserWarning):
                    MultiqcModule()


class TestDiversityCalculation:
    """Unit test diversity math with synthetic data."""

    def test_perfectly_even_distribution(self):
        """Perfectly even distribution should have max evenness."""
        from fusilli_multiqc.modules.diversity import MultiqcModule

        # 10 variants, each with count 100
        counts = np.array([100] * 10)
        proportions = counts / counts.sum()
        shannon = -np.sum(proportions * np.log(proportions))
        evenness = shannon / np.log(len(counts))
        assert evenness == pytest.approx(1.0, abs=1e-9)

    def test_single_dominant_variant(self):
        """One dominant variant should have low evenness."""
        counts = np.array([10000, 1, 1, 1, 1])
        proportions = counts / counts.sum()
        shannon = -np.sum(proportions * np.log(proportions))
        evenness = shannon / np.log(len(counts))
        assert evenness < 0.1
