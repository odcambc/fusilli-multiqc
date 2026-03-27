"""Shared fixtures for fusilli-multiqc tests."""

import os
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest


EXAMPLES_DIR = Path(__file__).parent.parent / "examples"


@pytest.fixture
def examples_dir():
    """Return the path to the examples directory."""
    return EXAMPLES_DIR


@pytest.fixture
def fusion_qc_metrics_df():
    """Load fusion_qc_metrics.csv from examples."""
    return pd.read_csv(EXAMPLES_DIR / "fusion_qc_metrics.csv")


@pytest.fixture
def sensitivity_metrics_df():
    """Load sensitivity_metrics.csv from examples."""
    return pd.read_csv(EXAMPLES_DIR / "sensitivity_metrics.csv")


@pytest.fixture
def decay_metrics_df():
    """Load decay_metrics.csv from examples."""
    return pd.read_csv(EXAMPLES_DIR / "decay_metrics.csv")


@pytest.fixture
def fusion_counts_summary_df():
    """Load fusion_counts_summary.csv from examples."""
    return pd.read_csv(EXAMPLES_DIR / "fusion_counts_summary.csv")


@pytest.fixture
def partner_counts_summary_df():
    """Load partner_counts_summary.csv from examples."""
    return pd.read_csv(EXAMPLES_DIR / "partner_counts_summary.csv")


def make_find_log_files(file_map):
    """
    Create a mock find_log_files function that yields file dicts
    based on a pattern-to-filepath mapping.

    Args:
        file_map: dict mapping pattern strings to file paths, e.g.
            {"fusilli_detection/fusion_qc_metrics": "/path/to/fusion_qc_metrics.csv"}
    """
    def find_log_files(pattern, **kwargs):
        if pattern in file_map:
            path = Path(file_map[pattern])
            yield {"root": str(path.parent), "fn": path.name}
    return find_log_files


def make_module_instance(module_class, file_map):
    """
    Instantiate a MultiQC module class with mocked BaseMultiqcModule methods.

    This patches find_log_files, add_data_source, add_section,
    general_stats_addcols, and the parent __init__ so that the module
    can be tested without a running MultiQC instance.

    Args:
        module_class: The MultiqcModule class to instantiate
        file_map: dict mapping search pattern to CSV file path

    Returns:
        The instantiated module
    """
    with patch.object(module_class, "__bases__", (object,)):
        # We need to patch at the class level to bypass super().__init__
        with patch("multiqc.base_module.BaseMultiqcModule.__init__", return_value=None):
            with patch.object(module_class, "find_log_files", side_effect=make_find_log_files(file_map)):
                with patch.object(module_class, "add_data_source"):
                    with patch.object(module_class, "add_section"):
                        with patch.object(module_class, "general_stats_addcols"):
                            instance = module_class()
    return instance
