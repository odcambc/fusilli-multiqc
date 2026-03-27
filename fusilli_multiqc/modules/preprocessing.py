"""
MultiQC module for FUSILLI preprocessing metrics.

This module visualizes:
- Read decay through preprocessing steps
- Retention rates per step
- Step loss breakdown
"""

import logging
from collections import OrderedDict
from pathlib import Path
from typing import Optional

import pandas as pd

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, bargraph

from fusilli_multiqc.utils import parse_csv_file

logger = logging.getLogger(__name__)

import os


class MultiqcModule(BaseMultiqcModule):
    """
    FUSILLI Preprocessing Metrics Module

    Parses decay_metrics.csv and step-specific metrics to visualize read retention
    through the preprocessing pipeline.
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="FUSILLI Preprocessing Metrics",
            anchor="fusilli_preprocessing",
            info="Read retention and loss metrics through preprocessing steps.",
        )

        # Find and parse input files
        self.decay_data = None
        self.trim_data = None
        self.contam_data = None
        self.quality_data = None

        # Search for decay_metrics.csv using registered pattern
        for f in self.find_log_files("fusilli_preprocessing/decay_metrics"):
            file_path = os.path.join(f["root"], f["fn"])
            self.decay_data = parse_csv_file(file_path)
            if self.decay_data is not None:
                self.add_data_source(f)
                break

        # Search for step-specific metrics (optional)
        for f in self.find_log_files("fusilli_preprocessing/trim_metrics"):
            file_path = os.path.join(f["root"], f["fn"])
            self.trim_data = parse_csv_file(file_path)
            if self.trim_data is not None:
                break

        for f in self.find_log_files("fusilli_preprocessing/contam_metrics"):
            file_path = os.path.join(f["root"], f["fn"])
            self.contam_data = parse_csv_file(file_path)
            if self.contam_data is not None:
                break

        for f in self.find_log_files("fusilli_preprocessing/quality_metrics"):
            file_path = os.path.join(f["root"], f["fn"])
            self.quality_data = parse_csv_file(file_path)
            if self.quality_data is not None:
                break

        # If no data found, exit
        if self.decay_data is None:
            raise UserWarning

        # Calculate retention rates
        self.retention_data = self.calculate_retention_rates()

        # Generate plots
        if self.decay_data is not None:
            self.read_decay_plot()
            self.base_decay_plot()

        if self.retention_data is not None and not self.retention_data.empty:
            # Only create plots if we have valid retention data
            try:
                self.retention_rate_plot()
            except (ValueError, KeyError) as e:
                logger.warning(f"Could not create retention rate plot: {e}")

            try:
                self.step_loss_plot()
            except (ValueError, KeyError) as e:
                logger.warning(f"Could not create step loss plot: {e}")

            try:
                self.add_summary_table()
            except (ValueError, KeyError) as e:
                logger.warning(f"Could not add summary table: {e}")

    def calculate_retention_rates(self) -> Optional[pd.DataFrame]:
        """
        Calculate retention rates for each preprocessing step.

        decay_metrics.csv stores all values in read pairs: pre-merge steps divide
        the combined R1+R2 count by 2; the merged step reports joined pairs directly.
        Retention is therefore simply pairs_at_step / raw_pairs for every step.

        Returns:
            DataFrame with columns: sample, step, retention_rate, or None if no data
        """
        if self.decay_data is None:
            return None

        raw_data = self.decay_data[self.decay_data["step"] == "raw"].copy()
        if raw_data.empty:
            return None

        raw_reads = dict(zip(raw_data["sample"], raw_data["reads"]))

        retention_list = []
        steps = ["trimmed", "cleaned", "quality", "merged", "matched"]

        for step in steps:
            step_data = self.decay_data[self.decay_data["step"] == step].copy()
            for _, row in step_data.iterrows():
                sample = row["sample"]
                reads = row["reads"]
                raw = raw_reads.get(sample, 0)
                retention = (reads / raw) if raw > 0 else 0.0
                retention_list.append({
                    "sample": sample,
                    "step": step,
                    "reads": reads,
                    "raw_reads": raw,
                    "retention_rate": retention,
                    "loss_rate": 1.0 - retention,
                })

        return pd.DataFrame(retention_list)

    def read_decay_plot(self) -> None:
        """Create read decay line plot through preprocessing steps."""
        if self.decay_data is None:
            return

        # Prepare data for line plot
        # MultiQC 1.33 expects data as dict of sample -> dict of x -> y
        # Use numeric indices so the x-axis respects pipeline order,
        # then map to labels via categories.
        plot_data = OrderedDict()
        steps = ["raw", "cleaned", "trimmed", "quality", "merged", "matched"]

        for sample in self.decay_data["sample"].unique():
            sample_data = self.decay_data[self.decay_data["sample"] == sample]
            sample_dict = OrderedDict()
            for step in steps:
                step_row = sample_data[sample_data["step"] == step]
                if not step_row.empty:
                    reads = step_row.iloc[0]["reads"]
                    # Convert to native Python type (not numpy)
                    if hasattr(reads, 'item'):
                        reads = reads.item()
                    sample_dict[step] = float(reads)
            if sample_dict:
                plot_data[sample] = sample_dict

        if not plot_data:
            return

        pconfig = {
            "id": "fusilli_read_decay",
            "title": "FUSILLI: Read Decay Through Preprocessing",
            "ylab": "Number of Reads",
            "xlab": "Preprocessing Step",
            "categories": True,
            "tt_label": "<b>{point.x}</b><br>{point.y:,.0f} reads",
            "ylog": True,
        }

        self.add_section(
            name="Read Decay",
            anchor="fusilli_read_decay",
            description="Read pairs at each preprocessing and detection step (log scale).",
            helptext="""
            This plot shows read pair counts through the full pipeline
            (in order: raw → cleaned → trimmed → quality → merged → matched).
            All values are in **read pairs**: pre-merge steps divide the combined
            R1+R2 count by 2; the merged and matched steps count individual reads/detections.
            With a high merge rate the line should remain nearly flat into the merged step.
            The drop to **matched** reflects the fraction of reads assigned to a fusion variant.
            - **raw**: Initial fragment pairs
            - **cleaned**: Pairs after contaminant removal
            - **trimmed**: Pairs after adapter trimming
            - **quality**: Pairs after quality filtering
            - **merged**: Pairs that successfully merged into a single read
            - **matched**: Merged reads matched to a fusion breakpoint
            """,
            plot=linegraph.plot(plot_data, pconfig),
        )

    def base_decay_plot(self) -> None:
        """Create base decay line plot through preprocessing steps."""
        if self.decay_data is None:
            return

        if "bases" not in self.decay_data.columns:
            return

        steps = ["raw", "cleaned", "trimmed", "quality", "merged"]
        plot_data = OrderedDict()

        for sample in self.decay_data["sample"].unique():
            sample_data = self.decay_data[self.decay_data["sample"] == sample]
            sample_dict = OrderedDict()
            for step in steps:
                step_row = sample_data[sample_data["step"] == step]
                if not step_row.empty:
                    bases = step_row.iloc[0]["bases"]
                    if pd.notna(bases):
                        val = bases
                        if hasattr(val, 'item'):
                            val = val.item()
                        sample_dict[step] = float(val)
            if sample_dict:
                plot_data[sample] = sample_dict

        if not plot_data:
            return

        pconfig = {
            "id": "fusilli_base_decay",
            "title": "FUSILLI: Base Decay Through Preprocessing",
            "ylab": "Number of Bases",
            "xlab": "Preprocessing Step",
            "categories": True,
            "tt_label": "<b>{point.x}</b><br>{point.y:,.0f} bases",
            "ylog": True,
        }

        self.add_section(
            name="Base Decay",
            anchor="fusilli_base_decay",
            description="Total bases at each preprocessing step (log scale).",
            helptext="""
            This plot shows total base counts through the preprocessing steps.
            Unlike read pairs, bases from R1 and R2 are tracked independently
            (i.e., not normalized — each read's bases are legitimately distinct).

            A larger drop in bases compared to reads at a given step indicates
            that shorter reads are being preferentially removed or trimmed.
            - **raw**: Total input bases across all reads
            - **cleaned**: Bases after contaminant removal
            - **trimmed**: Bases after adapter trimming (includes length reduction)
            - **quality**: Bases after quality filtering
            - **merged**: Bases in successfully merged reads
            """,
            plot=linegraph.plot(plot_data, pconfig),
        )

    def retention_rate_plot(self) -> None:
        """Create retention rate grouped bar chart."""
        if self.retention_data is None or self.retention_data.empty:
            return

        plot_data = OrderedDict()
        for _, row in self.retention_data.iterrows():
            sample = row["sample"]
            step = row["step"]
            if sample not in plot_data:
                plot_data[sample] = {}
            # Convert to native Python type
            retention_rate = row["retention_rate"]
            if hasattr(retention_rate, 'item'):
                retention_rate = retention_rate.item()
            plot_data[sample][step] = float(retention_rate)

        if not plot_data:
            return

        # Use explicit pipeline order so bars render as grouped (not overlapping)
        categories = ["trimmed", "cleaned", "quality", "merged", "matched"]

        if not categories or not plot_data:
            logger.debug("No retention rate data to plot")
            return

        pconfig = {
            "id": "fusilli_retention_rates",
            "title": "FUSILLI: Retention Rates by Step",
            "ylab": "Retention Rate",
            "xlab": "Sample",
            "stacking": None,
            "ymax": 1.0,
            "ymin": 0.0,
        }

        self.add_section(
            name="Retention Rates",
            anchor="fusilli_retention_rates",
            description="Read retention rates at each preprocessing step.",
            helptext="""
            This plot shows the fraction of reads retained after each preprocessing step:
            - **trimmed**: Retention after adapter trimming
            - **cleaned**: Retention after contaminant removal
            - **quality**: Retention after quality filtering
            - **merged**: Merge rate (fraction of pairs that successfully merged),
              i.e. 2×merged_reads / quality_reads. ~50% read count drop at merge
              is expected; here we show merge success, not read retention.

            Higher retention rates indicate less aggressive filtering.
            """,
            plot=bargraph.plot(plot_data, cats=categories, pconfig=pconfig),
        )

    def step_loss_plot(self) -> None:
        """Create step loss stacked bar chart."""
        if self.retention_data is None or self.retention_data.empty:
            return

        plot_data = OrderedDict()
        for sample in self.retention_data["sample"].unique():
            sample_data = self.retention_data[
                self.retention_data["sample"] == sample
            ]
            sample_dict = {}

            # Calculate loss at each step
            prev_retention = 1.0
            for step in ["trimmed", "cleaned", "quality", "merged", "matched"]:
                step_row = sample_data[sample_data["step"] == step]
                if not step_row.empty:
                    retention = step_row.iloc[0]["retention_rate"]
                    # Convert to native Python type
                    if hasattr(retention, 'item'):
                        retention = retention.item()
                    retention = float(retention)
                    loss = prev_retention - retention
                    sample_dict[step] = float(loss)
                    prev_retention = retention
                else:
                    sample_dict[step] = 0.0

            if sample_dict:
                plot_data[sample] = sample_dict

        if not plot_data:
            return

        # Use explicit pipeline order
        categories = ["trimmed", "cleaned", "quality", "merged"]

        if not categories:
            return

        pconfig = {
            "id": "fusilli_step_loss",
            "title": "FUSILLI: Step Loss Breakdown",
            "ylab": "Fraction of Raw Reads",
            "xlab": "Sample",
            "stacking": "normal",
        }

        self.add_section(
            name="Step Loss Breakdown",
            anchor="fusilli_step_loss",
            description="Fraction of reads lost at each preprocessing step.",
            helptext="""
            This stacked bar chart shows the cumulative loss of reads through preprocessing:
            - **trimmed**: Loss during adapter trimming
            - **cleaned**: Additional loss during contaminant removal
            - **quality**: Additional loss during quality filtering
            - **merged**: Fraction of pairs that did not merge (1 − merge rate).
              The expected ~50% drop in read count at merge is not shown as loss.

            The total height shows the final retention rate.
            """,
            plot=bargraph.plot(plot_data, cats=categories, pconfig=pconfig),
        )

    def add_summary_table(self) -> None:
        """Add retention rate metrics to general stats table."""
        if self.retention_data is None or self.retention_data.empty:
            return

        # Calculate final retention (merged step)
        merged_data = self.retention_data[
            self.retention_data["step"] == "merged"
        ].copy()

        if merged_data.empty:
            return

        headers = OrderedDict()
        headers["retention_rate"] = {
            "title": "Merge Rate",
            "description": "Fraction of pairs that merged (merged step); earlier steps show read retention",
            "format": "{:.3f}",
            "ymax": 1.0,
            "ymin": 0.0,
        }

        # Ensure DataFrame has 'sample' column as index for general_stats_addcols
        if "sample" in merged_data.columns:
            stats_data = merged_data.set_index("sample")
        else:
            stats_data = merged_data

        self.general_stats_addcols(stats_data, headers)
