"""
MultiQC module for FUSILLI detection metrics.

This module visualizes:
- Detection efficiency (detection_efficiency, prefilter_efficiency, matching_efficiency)
- Library coverage (variant, breakpoint, partner coverage)
- Sensitivity analysis (sensitivity_index vs expected_detection_fraction)
"""

import logging
from collections import OrderedDict
from pathlib import Path

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, bargraph

from fusilli_multiqc.utils import parse_csv_file

logger = logging.getLogger(__name__)

import os


class MultiqcModule(BaseMultiqcModule):
    """
    FUSILLI Detection Metrics Module

    Parses fusion_qc_metrics.csv, sensitivity_metrics.csv, and fusion_counts_summary.csv
    to generate detection efficiency, coverage, and sensitivity visualizations.
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="FUSILLI Detection Metrics",
            anchor="fusilli_detection",
            info="Fusion detection efficiency, library coverage, and sensitivity metrics.",
        )

        # Find and parse input files
        self.fusion_qc_data = None
        self.sensitivity_data = None
        self.fusion_counts_data = None

        # Search for files using registered patterns
        for f in self.find_log_files("fusilli_detection/fusion_qc_metrics"):
            file_path = os.path.join(f["root"], f["fn"])
            self.fusion_qc_data = parse_csv_file(file_path)
            if self.fusion_qc_data is not None:
                self.add_data_source(f)
                break

        for f in self.find_log_files("fusilli_detection/sensitivity_metrics"):
            file_path = os.path.join(f["root"], f["fn"])
            self.sensitivity_data = parse_csv_file(file_path)
            if self.sensitivity_data is not None:
                self.add_data_source(f)
                break

        for f in self.find_log_files("fusilli_detection/fusion_counts"):
            file_path = os.path.join(f["root"], f["fn"])
            self.fusion_counts_data = parse_csv_file(file_path)
            if self.fusion_counts_data is not None:
                self.add_data_source(f)
                break

        # Optional: decay_metrics for variant-called read retention (vs raw)
        self.decay_data = None
        for f in self.find_log_files("fusilli_preprocessing/decay_metrics"):
            file_path = os.path.join(f["root"], f["fn"])
            self.decay_data = parse_csv_file(file_path)
            break

        # If no data found, exit
        if self.fusion_qc_data is None and self.sensitivity_data is None:
            logger.debug("No detection data found - exiting module")
            raise UserWarning

        # Log what we found
        if self.fusion_qc_data is not None:
            logger.info(f"Loaded fusion_qc_data: {len(self.fusion_qc_data)} rows")
        if self.sensitivity_data is not None:
            logger.info(f"Loaded sensitivity_data: {len(self.sensitivity_data)} rows")

        # Calculate detection efficiency metrics
        if self.fusion_qc_data is not None:
            self.calculate_detection_efficiency()
            self.calculate_coverage_metrics()
            self.calculate_variant_called_retention()

        # Generate plots
        if self.fusion_qc_data is not None:
            self.detection_efficiency_plot()
            self.coverage_plot()
            self.variant_called_retention_plot()

        if self.sensitivity_data is not None:
            self.sensitivity_plot()

        # Add summary table
        if self.fusion_qc_data is not None:
            self.add_summary_table()

    def calculate_detection_efficiency(self) -> None:
        """Calculate detection efficiency metrics from fusion_qc_metrics."""
        if self.fusion_qc_data is None:
            return

        # Calculate detection efficiency if not present
        if "detection_efficiency" not in self.fusion_qc_data.columns:
            if "matched_reads" in self.fusion_qc_data.columns and "reads_processed" in self.fusion_qc_data.columns:
                self.fusion_qc_data["detection_efficiency"] = (
                    self.fusion_qc_data["matched_reads"] / self.fusion_qc_data["reads_processed"]
                ).fillna(0.0)

        # Calculate prefilter efficiency
        if "prefilter_efficiency" not in self.fusion_qc_data.columns:
            if "partner_end_reads" in self.fusion_qc_data.columns and "reads_processed" in self.fusion_qc_data.columns:
                self.fusion_qc_data["prefilter_efficiency"] = (
                    self.fusion_qc_data["partner_end_reads"] / self.fusion_qc_data["reads_processed"]
                ).fillna(0.0)

        # Calculate matching efficiency
        if "matching_efficiency" not in self.fusion_qc_data.columns:
            if "matched_reads" in self.fusion_qc_data.columns and "partner_end_reads" in self.fusion_qc_data.columns:
                self.fusion_qc_data["matching_efficiency"] = (
                    self.fusion_qc_data["matched_reads"] / self.fusion_qc_data["partner_end_reads"]
                ).fillna(0.0)

    def calculate_coverage_metrics(self) -> None:
        """Calculate library coverage metrics."""
        if self.fusion_qc_data is None:
            return

        # Calculate variant coverage
        if "variant_coverage" not in self.fusion_qc_data.columns:
            if "unique_fusions" in self.fusion_qc_data.columns and "expected_fusions" in self.fusion_qc_data.columns:
                self.fusion_qc_data["variant_coverage"] = (
                    self.fusion_qc_data["unique_fusions"] / self.fusion_qc_data["expected_fusions"]
                ).fillna(0.0)

        # Calculate breakpoint coverage
        # breakpoint_coverage uses a dedicated column if the pipeline provides it;
        # otherwise it falls back to variant_coverage as a reasonable proxy.
        if "breakpoint_coverage" not in self.fusion_qc_data.columns:
            if "unique_breakpoints" in self.fusion_qc_data.columns and "expected_breakpoints" in self.fusion_qc_data.columns:
                self.fusion_qc_data["breakpoint_coverage"] = (
                    self.fusion_qc_data["unique_breakpoints"] / self.fusion_qc_data["expected_breakpoints"]
                ).fillna(0.0)
            elif "variant_coverage" in self.fusion_qc_data.columns:
                self.fusion_qc_data["breakpoint_coverage"] = self.fusion_qc_data["variant_coverage"]

        # Calculate partner coverage
        if "partner_coverage" not in self.fusion_qc_data.columns:
            if "unique_partners_detected" in self.fusion_qc_data.columns:
                # Use expected_partners if available; otherwise fall back to the
                # per-sample max of unique_partners_detected (best available proxy).
                if "expected_partners" in self.fusion_qc_data.columns:
                    divisor = self.fusion_qc_data["expected_partners"]
                else:
                    max_detected = self.fusion_qc_data["unique_partners_detected"].max()
                    divisor = max_detected if max_detected > 0 else 1
                self.fusion_qc_data["partner_coverage"] = (
                    self.fusion_qc_data["unique_partners_detected"] / divisor
                ).clip(upper=1.0)

    def calculate_variant_called_retention(self) -> None:
        """
        Fraction of raw reads that had a variant (fusion) called.
        Requires decay_metrics (raw read count per sample) and fusion_qc (matched_reads).
        """
        if self.fusion_qc_data is None or self.decay_data is None:
            return
        raw = self.decay_data[self.decay_data["step"] == "raw"]
        if raw.empty or "sample" not in raw.columns or "reads" not in raw.columns:
            return
        raw_per_sample = dict(zip(raw["sample"], raw["reads"]))
        if "matched_reads" not in self.fusion_qc_data.columns or "sample" not in self.fusion_qc_data.columns:
            return
        retention = []
        for _, row in self.fusion_qc_data.iterrows():
            sample = row["sample"]
            matched = row.get("matched_reads", 0) or 0
            raw_reads = raw_per_sample.get(sample)
            if raw_reads and raw_reads > 0:
                retention.append(float(matched) / float(raw_reads))
            else:
                retention.append(0.0)
        self.fusion_qc_data["variant_called_retention"] = retention

    def variant_called_retention_plot(self) -> None:
        """Plot fraction of raw reads that had a variant called (decay metric for variant-called reads)."""
        if self.fusion_qc_data is None or "variant_called_retention" not in self.fusion_qc_data.columns:
            return
        plot_data = OrderedDict()
        for _, row in self.fusion_qc_data.iterrows():
            sample = row["sample"]
            r = row.get("variant_called_retention")
            if r is not None:
                plot_data[sample] = {"Variant-called retention": float(r)}
        if not plot_data:
            return
        pconfig = {
            "id": "fusilli_variant_called_retention",
            "title": "FUSILLI: Variant-called read retention (vs raw)",
            "ylab": "Fraction of raw reads",
            "xlab": "Sample",
            "ymax": 1.0,
            "ymin": 0.0,
        }
        self.add_section(
            name="Variant-called read retention",
            anchor="fusilli_variant_called_retention",
            description="Fraction of raw reads that had a fusion variant called.",
            helptext="""
            This is the decay/retention metric for the subset of reads that contribute to a variant call:
            **variant_called_retention** = matched_reads / raw_reads.

            It answers: of all raw reads, what fraction ended up assigned to a fusion?
            Low values indicate most reads are lost before or during detection.
            """,
            plot=bargraph.plot(
                plot_data, cats=["Variant-called retention"], pconfig=pconfig
            ),
        )

    def detection_efficiency_plot(self) -> None:
        """Create detection efficiency multi-line plot."""
        if self.fusion_qc_data is None or "sample" not in self.fusion_qc_data.columns:
            return

        plot_data = OrderedDict()
        for _, row in self.fusion_qc_data.iterrows():
            sample = row["sample"]
            plot_data[sample] = {}
            if "detection_efficiency" in row:
                plot_data[sample]["Detection Efficiency"] = row["detection_efficiency"]
            if "prefilter_efficiency" in row:
                plot_data[sample]["Prefilter Efficiency"] = row["prefilter_efficiency"]
            if "matching_efficiency" in row:
                plot_data[sample]["Matching Efficiency"] = row["matching_efficiency"]

        if not plot_data:
            return

        pconfig = {
            "id": "fusilli_detection_efficiency",
            "title": "FUSILLI: Detection Efficiency Metrics",
            "ylab": "Efficiency",
            "xlab": "Sample",
            "data_labels": [
                {"name": "Detection Efficiency", "ylab": "Efficiency"},
            ],
        }

        self.add_section(
            name="Detection Efficiency",
            anchor="fusilli_detection_efficiency",
            description="Detection efficiency metrics across samples.",
            helptext="""
            This plot shows three efficiency metrics:
            - **Detection Efficiency**: Fraction of processed reads that matched fusion sequences
            - **Prefilter Efficiency**: Fraction of reads containing partner domain 3' ends
            - **Matching Efficiency**: Fraction of prefiltered reads that matched specific breakpoints
            """,
            plot=linegraph.plot(plot_data, pconfig),
        )

    def coverage_plot(self) -> None:
        """Create library coverage line plot."""
        if self.fusion_qc_data is None or "sample" not in self.fusion_qc_data.columns:
            return

        # Each coverage metric becomes a line; samples are on the x-axis.
        # MultiQC linegraph: data = {series_name: {x_label: y_value}}
        coverage_metrics = [
            ("variant_coverage", "Variant Coverage"),
            ("breakpoint_coverage", "Breakpoint Coverage"),
            ("partner_coverage", "Partner Coverage"),
        ]

        plot_data = OrderedDict()
        for col, label in coverage_metrics:
            if col not in self.fusion_qc_data.columns:
                continue
            series = OrderedDict()
            for _, row in self.fusion_qc_data.iterrows():
                val = row[col]
                if hasattr(val, "item"):
                    val = val.item()
                series[row["sample"]] = float(val)
            if series:
                plot_data[label] = series

        if not plot_data:
            return

        pconfig = {
            "id": "fusilli_library_coverage",
            "title": "FUSILLI: Library Coverage Metrics",
            "ylab": "Coverage Fraction",
            "xlab": "Sample",
            "ymax": 1.0,
            "ymin": 0.0,
        }

        self.add_section(
            name="Library Coverage",
            anchor="fusilli_library_coverage",
            description="Library representation coverage across samples.",
            helptext="""
            This plot shows coverage metrics for different aspects of the fusion library:
            - **Variant Coverage**: Fraction of expected fusion variants detected (unique_fusions / expected_fusions)
            - **Breakpoint Coverage**: Fraction of expected breakpoint positions detected
            - **Partner Coverage**: Fraction of expected fusion partners detected

            **observed_variants vs unique_fusions**: In the pipeline, both denote the number of
            distinct fusion variants (fusion_id) with at least one read. They are the same for
            fusion-only QC. If unfused variants were included, observed_variants would be
            unique_fusions + unique_unfused (total variant types with ≥1 read).
            """,
            plot=linegraph.plot(plot_data, pconfig),
        )

    def sensitivity_plot(self) -> None:
        """Create detection efficiency bar chart showing matching and end-to-end efficiency."""
        if self.sensitivity_data is None or "sample" not in self.sensitivity_data.columns:
            return

        has_matching = "matching_efficiency" in self.sensitivity_data.columns
        has_e2e = "end_to_end_efficiency" in self.sensitivity_data.columns
        if not has_matching and not has_e2e:
            return

        plot_data = OrderedDict()
        for _, row in self.sensitivity_data.iterrows():
            sample = row["sample"]
            entry = {}
            if has_matching:
                entry["Matching Efficiency"] = float(row.get("matching_efficiency", 0.0))
            if has_e2e:
                entry["End-to-End Efficiency"] = float(row.get("end_to_end_efficiency", 0.0))
            plot_data[sample] = entry

        if not plot_data:
            return

        categories = []
        if has_matching:
            categories.append("Matching Efficiency")
        if has_e2e:
            categories.append("End-to-End Efficiency")

        pconfig = {
            "id": "fusilli_sensitivity",
            "title": "FUSILLI: Detection Efficiency",
            "ylab": "Efficiency",
            "ymax": 1.0,
            "ymin": 0.0,
            "stacking": None,
        }

        self.add_section(
            name="Detection Efficiency",
            anchor="fusilli_sensitivity",
            description="Matching efficiency and end-to-end efficiency per sample.",
            helptext="""
            Two complementary efficiency metrics, both corrected for the fraction of fragments
            geometrically capable of covering a breakpoint k-mer:

            - **Matching Efficiency**: detections / (merged reads × coverage probability).
              Measures how well the detector performs on reads that actually reached it.
            - **End-to-End Efficiency**: detections / (raw read pairs × coverage probability).
              Measures pipeline-wide yield from input material to final detection.

            The ratio between them equals the merge rate, which is also visible in the Read Decay plot.
            Values closer to 1.0 indicate better performance.
            """,
            plot=bargraph.plot(plot_data, cats=categories, pconfig=pconfig),
        )

    def add_summary_table(self) -> None:
        """Add detection yield summary table."""
        if self.fusion_qc_data is None or "sample" not in self.fusion_qc_data.columns:
            return

        # Prepare table data
        headers = OrderedDict()
        headers["reads_processed"] = {
            "title": "Reads Processed",
            "description": "Total reads processed for fusion detection",
            "format": "{:,.0f}",
        }
        headers["matched_reads"] = {
            "title": "Matched Reads",
            "description": "Reads that matched fusion sequences",
            "format": "{:,.0f}",
        }

        if "detections_per_read" in self.fusion_qc_data.columns:
            headers["detections_per_read"] = {
                "title": "Detections/Read",
                "description": "Number of detections per read",
                "format": "{:.4f}",
            }

        if "detections_per_million" in self.fusion_qc_data.columns:
            headers["detections_per_million"] = {
                "title": "Detections/Million",
                "description": "Number of detections per million reads",
                "format": "{:.2f}",
            }

        if "detection_efficiency" in self.fusion_qc_data.columns:
            headers["detection_efficiency"] = {
                "title": "Detection Efficiency",
                "description": "Fraction of reads that matched fusions",
                "format": "{:.3f}",
                "ymax": 1.0,
                "ymin": 0.0,
            }

        if "variant_called_retention" in self.fusion_qc_data.columns:
            headers["variant_called_retention"] = {
                "title": "Variant retention (vs raw)",
                "description": "Fraction of raw reads that had a fusion variant called",
                "format": "{:.4f}",
                "ymax": 1.0,
                "ymin": 0.0,
            }

        if "unique_fusions" in self.fusion_qc_data.columns:
            headers["unique_fusions"] = {
                "title": "Unique fusions",
                "description": "Distinct fusion variants (fusion_id) with ≥1 read; same as observed_variants for fusions",
                "format": "{:,.0f}",
            }

        if "variant_coverage" in self.fusion_qc_data.columns:
            headers["variant_coverage"] = {
                "title": "Variant Coverage",
                "description": "Fraction of expected fusion variants detected",
                "format": "{:.3f}",
                "ymax": 1.0,
                "ymin": 0.0,
            }

        if "partner_coverage" in self.fusion_qc_data.columns:
            headers["partner_coverage"] = {
                "title": "Partner Coverage",
                "description": "Fraction of expected fusion partners detected",
                "format": "{:.3f}",
                "ymax": 1.0,
                "ymin": 0.0,
            }

        if "sample" in self.fusion_qc_data.columns:
            stats_data = self.fusion_qc_data.set_index("sample")
        else:
            stats_data = self.fusion_qc_data

        self.general_stats_addcols(stats_data.to_dict(orient="index"), headers)
