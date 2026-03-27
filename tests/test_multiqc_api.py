"""
MultiQC API guard tests.

These tests verify that the MultiQC APIs we depend on haven't changed
in a way that would break our plugin. If any of these fail after a
MultiQC upgrade, the plugin likely needs updating.
"""

import inspect

import pytest
from multiqc.base_module import BaseMultiqcModule
from multiqc import config


class TestBaseMultiqcModuleAPI:
    """Verify BaseMultiqcModule has the interface we rely on."""

    def test_init_signature(self):
        """BaseMultiqcModule.__init__ must accept name, anchor, info kwargs."""
        sig = inspect.signature(BaseMultiqcModule.__init__)
        params = list(sig.parameters.keys())
        assert "name" in params, "BaseMultiqcModule.__init__ missing 'name' parameter"
        assert "anchor" in params, "BaseMultiqcModule.__init__ missing 'anchor' parameter"
        assert "info" in params, "BaseMultiqcModule.__init__ missing 'info' parameter"

    def test_has_find_log_files(self):
        assert hasattr(BaseMultiqcModule, "find_log_files"), "Missing find_log_files method"
        assert callable(getattr(BaseMultiqcModule, "find_log_files"))

    def test_has_add_data_source(self):
        assert hasattr(BaseMultiqcModule, "add_data_source"), "Missing add_data_source method"

    def test_has_add_section(self):
        assert hasattr(BaseMultiqcModule, "add_section"), "Missing add_section method"

    def test_has_general_stats_addcols(self):
        assert hasattr(BaseMultiqcModule, "general_stats_addcols"), "Missing general_stats_addcols method"


class TestMultiqcConfigAPI:
    """Verify multiqc.config has the attributes we depend on."""

    def test_config_has_sp(self):
        """config.sp must exist (or be settable) for search pattern registration."""
        # Our hooks set config.sp; it should either already exist or be settable
        assert hasattr(config, "sp") or not hasattr(config, "__slots__"), \
            "Cannot set config.sp — MultiQC config API may have changed"


class TestMultiqcPlotsAPI:
    """Verify multiqc.plots modules we import still exist."""

    def test_bargraph_importable(self):
        from multiqc.plots import bargraph
        assert hasattr(bargraph, "plot")

    def test_linegraph_importable(self):
        from multiqc.plots import linegraph
        assert hasattr(linegraph, "plot")

    def test_scatter_importable(self):
        from multiqc.plots import scatter
        assert hasattr(scatter, "plot")

    def test_heatmap_importable(self):
        from multiqc.plots import heatmap
        assert hasattr(heatmap, "plot")


class TestEntryPoints:
    """Verify our entry points are properly configured."""

    def test_module_entry_points_importable(self):
        """All four module classes should be importable."""
        from fusilli_multiqc.modules.detection import MultiqcModule as Detection
        from fusilli_multiqc.modules.diversity import MultiqcModule as Diversity
        from fusilli_multiqc.modules.preprocessing import MultiqcModule as Preprocessing
        from fusilli_multiqc.modules.partners import MultiqcModule as Partners

        for cls in [Detection, Diversity, Preprocessing, Partners]:
            assert issubclass(cls, BaseMultiqcModule), \
                f"{cls.__name__} is not a subclass of BaseMultiqcModule"

    def test_hooks_entry_point_importable(self):
        from fusilli_multiqc.hooks import execution_start
        assert callable(execution_start)
