# fusilli-multiqc

Custom MultiQC modules for the FUSILLI fusion library pipeline.

This package is made to be run with the FUSILLI pipeline: https://github.com/user/fusilli. It provides metrics and visualizations specific to the pipeline's outputs.

- **Detection Metrics**: Detection efficiency, library coverage, and sensitivity analysis
- **Diversity Metrics**: Library diversity metrics
- **Processing Metrics**: Read and base retention and loss through processing steps
- **Partner Detection**: Partner detection metrics and coverage across samples

## Installation

### Through FUSILLI

The FUSILLI pipeline uses a conda environment to install multiqc and this module. You should not need to install this module separately if you are using FUSILLI.

### Through pip

```bash
pip install fusilli-multiqc
```

## Requirements

- Python >= 3.8
- MultiQC >= 1.0
- pandas >= 1.0
- numpy >= 1.0

## Usage

After installation, the modules will be automatically discovered by MultiQC when you run:

```bash
multiqc <results_directory>
```

The modules will automatically detect and parse the following FUSILLI output files:

- `fusion_qc_metrics.csv` - Detection efficiency and coverage metrics
- `sensitivity_metrics.csv` - Sensitivity analysis data
- `fusion_counts_summary.csv` - Fusion count data for diversity analysis
- `decay_metrics.csv` - Read decay through preprocessing steps
- `partner_counts_summary.csv` - Partner detection counts

## Module Descriptions

### Detection Metrics Module (`fusilli_detection`)

Visualizes:
- **Detection Efficiency**: Fraction of processed reads that matched fusion sequences
- **Prefilter Efficiency**: Fraction of reads containing partner domain 3' ends
- **Matching Efficiency**: Fraction of prefiltered reads that matched specific breakpoints
- **Library Coverage**: Variant, breakpoint, and partner coverage metrics
- **Sensitivity Analysis**: Sensitivity index vs expected detection fraction

### Diversity Metrics Module (`fusilli_diversity`)

Visualizes:
- **Shannon Diversity**: Measures both richness and evenness
- **Simpson Diversity**: Measures dominance
- **Evenness**: Pielou's evenness index
- **Top N Fractions**: Distribution of counts across top variants
- **Variant Count Distribution**: Histogram of variant abundances

### Preprocessing Metrics Module (`fusilli_preprocessing`)

Visualizes:
- **Read Decay**: Read counts at each preprocessing step (log scale)
- **Retention Rates**: Fraction of reads retained after each step
- **Step Loss Breakdown**: Cumulative loss of reads through preprocessing

### Partner Detection Module (`fusilli_partners`)

Visualizes:
- **Partner Detection Heatmap**: Binary detection matrix of partners across samples
- **Partner Coverage**: Number of samples with detection per partner
- **Partner End vs Linker**: Comparison of detection methods

## License

MIT License - see LICENSE file for details.