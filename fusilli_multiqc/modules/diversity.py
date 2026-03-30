"""
MultiQC module for FUSILLI diversity metrics.

This module visualizes:
- Shannon and Simpson diversity indices
- Evenness (Pielou's evenness)
- Top N fractions (top1, top5, top10, top25)
- Variant count distribution
"""

import logging
from collections import OrderedDict
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, bargraph

from fusilli_multiqc.utils import parse_csv_file

logger = logging.getLogger(__name__)

import os


class MultiqcModule(BaseMultiqcModule):
    """
    FUSILLI Diversity Metrics Module

    Parses fusion_counts_summary.csv to calculate and visualize library diversity metrics.
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="FUSILLI Diversity Metrics",
            anchor="fusilli_diversity",
            info="Library diversity metrics including Shannon/Simpson diversity and evenness.",
        )

        # Find and parse input files
        self.fusion_counts_data = None

        # Search for fusion_counts_summary.csv using registered pattern
        for f in self.find_log_files("fusilli_diversity/fusion_counts"):
            file_path = os.path.join(f["root"], f["fn"])
            self.fusion_counts_data = parse_csv_file(file_path)
            if self.fusion_counts_data is not None:
                self.add_data_source(f)
                break

        # If no data found, exit
        if self.fusion_counts_data is None:
            raise UserWarning

        # Calculate diversity metrics
        self.diversity_metrics = self.calculate_diversity_metrics()

        # Generate plots
        if self.diversity_metrics is not None and not self.diversity_metrics.empty:
            self.evenness_plot()
            self.rank_abundance_plot()
            self.variant_count_distribution_plot()
            self.add_summary_table()

    def calculate_diversity_metrics(self) -> Optional[pd.DataFrame]:
        """
        Calculate diversity metrics for each sample.

        Returns:
            DataFrame with columns: sample, shannon_diversity, simpson_diversity, evenness,
            top1_fraction, top5_fraction, top10_fraction, top25_fraction, or None if no data
        """
        if self.fusion_counts_data is None:
            return None

        # Get sample columns (exclude fusion_id, type, total)
        sample_cols = [
            c for c in self.fusion_counts_data.columns
            if c not in ["fusion_id", "type", "total"]
        ]

        if not sample_cols:
            return None

        # Filter to fusion type only
        fusion_data = self.fusion_counts_data[
            self.fusion_counts_data["type"] == "fusion"
        ].copy()

        if fusion_data.empty:
            return None

        metrics_list = []

        for sample in sample_cols:
            # Get counts for this sample
            counts = fusion_data[sample].fillna(0).values
            total_counts = counts.sum()

            if total_counts == 0:
                metrics_list.append({
                    "sample": sample,
                    "evenness": 0.0,
                    "top1_fraction": 0.0,
                    "top5_fraction": 0.0,
                    "top10_fraction": 0.0,
                    "top25_fraction": 0.0,
                    "observed_variants": 0,
                })
                continue

            # Filter out zero counts
            non_zero_counts = counts[counts > 0]
            if len(non_zero_counts) == 0:
                metrics_list.append({
                    "sample": sample,
                    "evenness": 0.0,
                    "top1_fraction": 0.0,
                    "top5_fraction": 0.0,
                    "top10_fraction": 0.0,
                    "top25_fraction": 0.0,
                    "observed_variants": 0,
                })
                continue

            # Calculate proportions
            proportions = non_zero_counts / total_counts

            # Shannon diversity: H' = -Σ(p_i * ln(p_i))
            shannon = -np.sum(proportions * np.log(proportions))

            # Evenness: H' / ln(S) where S is number of observed variants
            observed_variants = len(non_zero_counts)
            if observed_variants > 1:
                evenness = shannon / np.log(observed_variants)
            else:
                evenness = 0.0

            # Top N fractions
            sorted_counts = np.sort(non_zero_counts)[::-1]  # Descending order
            top1 = sorted_counts[0] if len(sorted_counts) > 0 else 0
            top5 = sorted_counts[:5].sum() if len(sorted_counts) >= 5 else sorted_counts.sum()
            top10 = sorted_counts[:10].sum() if len(sorted_counts) >= 10 else sorted_counts.sum()
            top25 = sorted_counts[:25].sum() if len(sorted_counts) >= 25 else sorted_counts.sum()

            metrics_list.append({
                "sample": sample,
                "evenness": float(evenness),
                "top1_fraction": float(top1 / total_counts),
                "top5_fraction": float(top5 / total_counts),
                "top10_fraction": float(top10 / total_counts),
                "top25_fraction": float(top25 / total_counts),
                "observed_variants": int(observed_variants),
            })

        return pd.DataFrame(metrics_list)

    def evenness_plot(self) -> None:
        """Create evenness bar chart."""
        if self.diversity_metrics is None or self.diversity_metrics.empty:
            return

        plot_data = OrderedDict()
        for _, row in self.diversity_metrics.iterrows():
            sample = row["sample"]
            plot_data[sample] = {"Evenness": row["evenness"]}

        if not plot_data:
            return

        categories = list(set().union(*[d.keys() for d in plot_data.values()])) if plot_data else []

        pconfig = {
            "id": "fusilli_evenness",
            "title": "FUSILLI: Library Evenness",
            "ylab": "Evenness (Pielou's)",
            "xlab": "Sample",
            "ymax": 1.0,
            "ymin": 0.0,
        }

        self.add_section(
            name="Evenness",
            anchor="fusilli_evenness",
            description="Pielou's evenness index across samples.",
            helptext="""
            Evenness measures how evenly distributed variant counts are.
            Values range from 0 (highly uneven) to 1 (perfectly even).
            Higher evenness indicates more uniform representation of variants.
            """,
            plot=bargraph.plot(plot_data, cats=categories, pconfig=pconfig),
        )

    def rank_abundance_plot(self) -> None:
        """Create rank-abundance (Whittaker) curve: one line per sample."""
        if self.fusion_counts_data is None:
            return

        fusion_data = self.fusion_counts_data[
            self.fusion_counts_data["type"] == "fusion"
        ].copy()
        if fusion_data.empty:
            return

        sample_cols = [
            c for c in fusion_data.columns
            if c not in ["fusion_id", "type", "total"]
        ]
        if not sample_cols:
            return

        plot_data = OrderedDict()
        for sample in sample_cols:
            counts = fusion_data[sample].fillna(0).values
            counts = counts[counts > 0]
            if len(counts) == 0:
                continue
            total = counts.sum()
            sorted_props = np.sort(counts)[::-1] / total
            plot_data[sample] = {int(rank + 1): float(p) for rank, p in enumerate(sorted_props)}

        if not plot_data:
            return

        pconfig = {
            "id": "fusilli_rank_abundance",
            "title": "FUSILLI: Rank-Abundance Curve",
            "ylab": "Relative Abundance",
            "xlab": "Variant Rank",
            "ylog": True,
            "xlog": True,
            "ymin": 0.0,
        }

        self.add_section(
            name="Rank-Abundance",
            anchor="fusilli_rank_abundance",
            description="Relative abundance of each variant by rank (most to least abundant).",
            helptext="""
            Each line represents one sample. Variants are ranked from most to least abundant
            (rank 1 = highest count). The y-axis shows each variant's fraction of total counts
            for that sample (log scale).

            A steep drop indicates the library is dominated by a few variants; a flatter
            curve indicates more even representation. Samples with similar curves have
            similar diversity profiles.
            """,
            plot=linegraph.plot(plot_data, pconfig),
        )

    def variant_count_distribution_plot(self) -> None:
        """Create variant count distribution histogram."""
        if self.fusion_counts_data is None:
            return

        # Get fusion counts (exclude type column and sample columns for now)
        fusion_data = self.fusion_counts_data[
            self.fusion_counts_data["type"] == "fusion"
        ].copy()

        if fusion_data.empty:
            return

        # Get sample columns
        sample_cols = [
            c for c in fusion_data.columns
            if c not in ["fusion_id", "type", "total"]
        ]

        if not sample_cols:
            return

        # Aggregate counts across all samples
        all_counts = []
        for sample in sample_cols:
            counts = fusion_data[sample].fillna(0).values
            all_counts.extend(counts[counts > 0].tolist())

        if not all_counts:
            return

        # Create histogram data (log scale bins)
        log_counts = np.log10(np.array(all_counts) + 1)
        bins = np.linspace(0, np.ceil(log_counts.max()), 20)
        hist, bin_edges = np.histogram(log_counts, bins=bins)

        # Convert to plot format
        plot_data = OrderedDict()
        for i in range(len(hist)):
            bin_label = f"10^{bin_edges[i]:.1f} - 10^{bin_edges[i+1]:.1f}"
            plot_data[bin_label] = {"Count": int(hist[i])}

        categories = list(set().union(*[d.keys() for d in plot_data.values()])) if plot_data else []

        pconfig = {
            "id": "fusilli_variant_distribution",
            "title": "FUSILLI: Variant Count Distribution",
            "ylab": "Number of Variants",
            "xlab": "Count (log10 scale)",
        }

        self.add_section(
            name="Variant Count Distribution",
            anchor="fusilli_variant_distribution",
            description="Distribution of variant counts (log scale).",
            helptext="""
            This histogram shows the distribution of variant abundances across all samples.
            The x-axis is on a log10 scale to better visualize the wide range of counts.
            A more uniform distribution indicates better library representation.
            """,
            plot=bargraph.plot(plot_data, cats=categories, pconfig=pconfig),
        )

    def add_summary_table(self) -> None:
        """Add diversity metrics to general stats table."""
        if self.diversity_metrics is None or self.diversity_metrics.empty:
            return

        headers = OrderedDict()
        headers["evenness"] = {
            "title": "Evenness",
            "description": "Pielou's evenness index",
            "format": "{:.3f}",
            "ymax": 1.0,
            "ymin": 0.0,
        }
        headers["observed_variants"] = {
            "title": "Observed Variants",
            "description": "Number of unique variants detected",
            "format": "{:,.0f}",
        }
        headers["top1_fraction"] = {
            "title": "Top Variant Fraction",
            "description": "Fraction of total counts from the single most abundant variant",
            "format": "{:.3f}",
            "ymax": 1.0,
            "ymin": 0.0,
        }

        if "sample" in self.diversity_metrics.columns:
            stats_data = self.diversity_metrics.set_index("sample")
        else:
            stats_data = self.diversity_metrics

        self.general_stats_addcols(stats_data.to_dict(orient="index"), headers)
