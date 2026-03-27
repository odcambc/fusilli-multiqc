"""Tests for fusilli_multiqc.hooks."""

from unittest.mock import MagicMock, patch


class TestExecutionStart:
    def test_registers_all_search_patterns(self):
        """execution_start should register all expected search patterns on config.sp."""
        mock_config = MagicMock()
        mock_config.sp = {}

        with patch("fusilli_multiqc.hooks.config", mock_config):
            from fusilli_multiqc.hooks import execution_start
            execution_start()

        expected_patterns = [
            "fusilli_detection/fusion_qc_metrics",
            "fusilli_detection/sensitivity_metrics",
            "fusilli_detection/fusion_counts",
            "fusilli_preprocessing/decay_metrics",
            "fusilli_preprocessing/trim_metrics",
            "fusilli_preprocessing/contam_metrics",
            "fusilli_preprocessing/quality_metrics",
            "fusilli_diversity/fusion_counts",
            "fusilli_partners/partner_counts",
            "fusilli_partners/unmerged_counts",
            "fusilli_partners/unmerged_partner_counts",
        ]

        for pattern in expected_patterns:
            assert pattern in mock_config.sp, f"Missing search pattern: {pattern}"

    def test_patterns_have_fn_key(self):
        """Each registered pattern must have an 'fn' key with a glob."""
        mock_config = MagicMock()
        mock_config.sp = {}

        with patch("fusilli_multiqc.hooks.config", mock_config):
            from fusilli_multiqc.hooks import execution_start
            execution_start()

        for key, pattern in mock_config.sp.items():
            assert "fn" in pattern, f"Pattern {key} missing 'fn' key"
            assert pattern["fn"].endswith(".csv"), f"Pattern {key} fn doesn't match CSV"

    def test_shared_patterns(self):
        """fusion_counts_summary is shared between detection and diversity modules."""
        mock_config = MagicMock()
        mock_config.sp = {}

        with patch("fusilli_multiqc.hooks.config", mock_config):
            from fusilli_multiqc.hooks import execution_start
            execution_start()

        assert mock_config.sp["fusilli_detection/fusion_counts"].get("shared") is True
        assert mock_config.sp["fusilli_diversity/fusion_counts"].get("shared") is True

    def test_creates_sp_if_missing(self):
        """If config has no sp attribute, execution_start should create it."""
        mock_config = MagicMock(spec=[])  # no attributes at all

        with patch("fusilli_multiqc.hooks.config", mock_config):
            from fusilli_multiqc.hooks import execution_start
            execution_start()

        assert hasattr(mock_config, "sp")
