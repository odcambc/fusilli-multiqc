"""
MultiQC hooks for FUSILLI plugin.

This module registers search patterns for FUSILLI files BEFORE MultiQC
scans for files. This is critical - patterns must be registered via hooks,
not at module import time.

Pattern keys must begin with the module name (e.g., fusilli_detection/pattern_name).
"""

from multiqc import config


def execution_start():
    """
    Register search patterns at execution start, before file scanning begins.

    This hook runs early in MultiQC's initialization, allowing us to
    register file patterns that MultiQC will use during its file search phase.
    Pattern keys must begin with the module name.
    """
    # Initialize sp dict if needed
    if not hasattr(config, "sp"):
        config.sp = {}

    # Detection module patterns (module name: fusilli_detection)
    config.sp["fusilli_detection/fusion_qc_metrics"] = {
        "fn": "*fusion_qc_metrics.csv",
    }
    config.sp["fusilli_detection/sensitivity_metrics"] = {
        "fn": "*sensitivity_metrics.csv",
    }
    config.sp["fusilli_detection/fusion_counts"] = {
        "fn": "*fusion_counts_summary.csv",
        "shared": True,
    }

    # Preprocessing module patterns (module name: fusilli_preprocessing)
    config.sp["fusilli_preprocessing/decay_metrics"] = {
        "fn": "*decay_metrics.csv",
    }
    config.sp["fusilli_preprocessing/trim_metrics"] = {
        "fn": "*trim_metrics.csv",
    }
    config.sp["fusilli_preprocessing/contam_metrics"] = {
        "fn": "*contam_metrics.csv",
    }
    config.sp["fusilli_preprocessing/quality_metrics"] = {
        "fn": "*quality_metrics.csv",
    }

    # Diversity module patterns (module name: fusilli_diversity)
    config.sp["fusilli_diversity/fusion_counts"] = {
        "fn": "*fusion_counts_summary.csv",
        "shared": True,
    }

    # Partners module patterns (module name: fusilli_partners)
    config.sp["fusilli_partners/partner_counts"] = {
        "fn": "*partner_counts_summary.csv",
    }
    config.sp["fusilli_partners/unmerged_counts"] = {
        "fn": "*unmerged_counts_summary.csv",
    }
    config.sp["fusilli_partners/unmerged_partner_counts"] = {
        "fn": "*unmerged_partner_counts_summary.csv",
    }
