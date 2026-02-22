# PrPC CRISPRa Screen Pipeline Documentation

## 1. Purpose and Scope

This document explains:
- what this pipeline does,
- how to install and run it,
- how to map outputs to manuscript figures,
- how to reproduce manuscript-level analysis products from raw/processed data.

Manuscript context: *Genome-wide arrayed CRISPR activation screen for prion protein modulators*.

Core reported study parameters from the manuscript:
- Cell line: U-251 MG (dCas9-VPR)
- Screen type: arrayed CRISPRa gain-of-function
- Scale: 19,839 genes, 22,442 perturbations
- Plate format: 384-well
- Readout: TR-FRET (Eu-POM2 donor, APC-POM1 acceptor)
- Replicates: duplicate measurements per perturbation
- Hit threshold: `|log2FC| > 1` and `p < 0.05`

---

## 2. Repository Structure

- `prpcscreen/analysis/`
  - `calculating_scores.py`: normalization, SSMD, p-values, Z-prime/SSMD control metrics
  - `processing_data.py`: feature engineering, hit calling
- `prpcscreen/misc/`
  - `converting_plate_layouts.py`: 384->96 mapping utilities
- `prpcscreen/visualization/`
  - modules for QC, plate-series, replicate scatter, histograms, volcano/flashlight, heatmaps, violin/box plots
- `prpcscreen/scripts/`
  - CLI scripts for each workflow stage
- `prpcscreen/scripts/plot_genomic_signal_skyline.py`
  - genomic localization / Skyline plotting

---

## 3. Environment Setup

### 3.1 Requirements

- Python 3.10+ (3.11 recommended)
- OS: Windows/macOS/Linux
- `pip`

### 3.2 Installation

```powershell
# from repository root
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install --upgrade pip
pip install -r requirements.txt
```

Packages used:
- pandas
- numpy
- scipy
- matplotlib
- seaborn
- openpyxl

---

## 4. Input Data Requirements

The scripts assume tabular inputs with specific columns.

### 4.1 Layout / annotation table (CSV)

Required columns (minimum for end-to-end analysis):
- `Plate_number_384`
- `Well_number_384`
- `Is_NT_ctrl`
- `Is_pos_ctrl`
- `Entrez_ID`
- `Target_flag` (needed for own-Non-targeting-aware methods)

Recommended additional columns:
- `Gene_symbol`
- `x_position`
- `Target_ID`

### 4.2 Raw assay files

`prpcscreen/scripts/merge_assay_exports.py` searches recursively for:
- TR-FRET files matching `*TR-FRET*.csv`
- CellTiterGlo files matching `*GLO*.csv`

Defaults:
- `--skip-fret 38`
- `--skip-glo 9`

Adjust these if your export headers differ.

### 4.3 Skyline plot input

Excel file/sheet must contain:
- `Gene_symbol`
- `Mean_log2FC`
- `Chromosome`
- `Start_Position`

---

## 5. What Each Stage Does

### 5.1 Plate layout conversion

Script: `prpcscreen/scripts/remap_plate_coordinates.py`

Function:
- Converts 384-well indices to 96-well mapping using the A1/A2/B1/B2 plate scheme.

---

### 5.2 Data integration

Script: `prpcscreen/scripts/merge_assay_exports.py`

Function:
- Reads raw TR-FRET and GLO exports,
- aligns them to the layout table,
- writes integrated per-well dataset.

Key outputs added:
- `Raw_rep1`
- `Raw_rep2`
- `CellTiterGlo_raw`

---

### 5.3 Normalization and statistics

Script: `prpcscreen/scripts/compute_screen_metrics.py`

Function:
- Applies per-plate normalization (default `norm_method="all Non-targeting"`),
- derives transformed metrics,
- computes SSMD and p-values,
- writes analyzed table and optional hit list.

Important derived columns:
- `Log2FC_rep1`, `Log2FC_rep2`
- `Mean_log2`
- `p_value_log2`
- `SSMD_log2`
- GLO-adjusted variants (`*_Glo_*`)

Hit extraction (default):
- `p_value_log2 < 0.05`
- `|Mean_log2| > 1`

---

### 5.4 QC and plotting

Scripts:
- `plot_plate_health.py`
- `plot_well_trajectories.py`
- `plot_replicate_agreement.py`
- `plot_signal_distributions.py`
- `plot_candidate_landscape.py`
- `plot_spatial_and_group_views.py`
- `prpcscreen/scripts/plot_genomic_signal_skyline.py`

Functions:
- Plate-level control separation metrics (SSMD controls)
- Well-series trend inspection
- Replicate agreement plots
- Distribution histograms
- Volcano and ranking plots
- Plate heatmaps and grouped violin/box plots
- Genome-positioned Skyline-style summary

---

## 6. End-to-End Usage

Below is a practical runbook from raw data to figure outputs.

### 6.1 One-command orchestration

Use the root script `orchestrate_screen_workflow.ps1` to execute all stages in sequence:

```powershell
.\orchestrate_screen_workflow.ps1 `
  -RawDir data/raw_prp `
  -LayoutCsv data/layout/layout_384.csv `
  -GenomicsExcel data/genomics/PrP_genes_and_NT_ordered_with_chromosome.xlsx
```

Optional flags:
- `-OutputDir results`
- `-SkylineSheet skylineplot2`
- `-SkipFret 38`
- `-SkipGlo 9`
- `-HeatmapPlate 1` (or range `1-4`, series `1,2,6`, or `all`)

### 6.2 Stepwise execution

```powershell
# 0) activate env
.\.venv\Scripts\Activate.ps1

# 1) integrate raw data
python prpcscreen/scripts/merge_assay_exports.py `
  data/raw_prp `
  data/layout/layout_384.csv `
  results/01_integrated.csv

# 2) analyze / score / call hits
python prpcscreen/scripts/compute_screen_metrics.py `
  results/01_integrated.csv `
  results/02_analyzed.csv `
  --hits_csv results/03_hits.csv

# 3) QC + primary figure-like outputs
python prpcscreen/scripts/plot_plate_health.py results/02_analyzed.csv results/figures/plate_qc_ssmd_controls.png
python prpcscreen/scripts/plot_well_trajectories.py results/02_analyzed.csv results/figures/plate_well_series_raw_rep1.png --column Raw_rep1
python prpcscreen/scripts/plot_replicate_agreement.py results/02_analyzed.csv results/figures/replicate_agreement_log2fc.png --stem Log2FC
python prpcscreen/scripts/plot_signal_distributions.py results/02_analyzed.csv --output_html results/figures/distribution_log2fc_rep1_interactive.html --column Log2FC_rep1
python prpcscreen/scripts/plot_candidate_landscape.py results/02_analyzed.csv results/figures/candidate_flashlight_ranked_meanlog2.png --volcano_html results/figures/candidate_volcano_interactive.html --genomics_excel data/genomics/PrP_genes_and_NT_ordered_with_chromosome.xlsx
# Note: if --genomics_excel is omitted or has no Sublibrary column, volcano filtering falls back to prpcscreen/misc/supplementary_sublibrary_map.csv (from supplementary workbook 41551_2024_1278_MOESM4_ESM.xlsx).
python prpcscreen/scripts/plot_spatial_and_group_views.py results/02_analyzed.csv results/figures/plate_heatmap_raw_rep1.png results/figures/grouped_boxplot_raw_rep1.png --plate 1

# 4) Skyline / genomic localization
python prpcscreen/scripts/plot_genomic_signal_skyline.py data/genomics/PrP_genes_and_NT_ordered_with_chromosome.xlsx results/figures/genomic_skyline_meanlog2fc.png --sheet skylineplot2
```

---

## 7. Mapping Manuscript Figures to Pipeline Outputs

This mapping is based on manuscript descriptions and available scripts.

- Figure 3C (representative plate heatmap):
  - Use `prpcscreen/scripts/plot_spatial_and_group_views.py` (heatmap output)
- Figure 4A (volcano):
  - Use `prpcscreen/scripts/plot_candidate_landscape.py` (volcano output)
- Figure 5A (plate quality):
  - Use `prpcscreen/scripts/plot_plate_health.py`
- Figure 5B (replicate concordance):
  - Use `prpcscreen/scripts/plot_replicate_agreement.py --stem Log2FC`
- Figure 5C (histogram distributions):
  - Use `prpcscreen/scripts/plot_signal_distributions.py --column Log2FC_rep1`
- Figure 8B (Skyline):
  - Use `prpcscreen/scripts/plot_genomic_signal_skyline.py`

Notes:
- Figure 6 (secondary TR-FRET + WB validation) depends on external validation datasets not generated by this pipeline alone.
- Figure 7 (lentiviral titer heatmaps/distributions) requires dedicated lentiviral titer input tables and custom plotting not included in current CLI wrappers.
- Figure 8A/C/D (GWAS overlap and chromosome-20 focused analysis) requires additional GWAS/annotation integration scripts not bundled in current wrappers.

---

## 8. Reproducibility Checklist

Before claiming full reproduction, verify:

1. Input integrity
- Plate/well counts match expected 384-per-plate structure.
- Control annotations (`Is_NT_ctrl`, `Is_pos_ctrl`) are populated.

2. QC behavior
- Positive controls are separated from Non-targeting controls in QC plots.
- Replicate scatter has high concordance (manuscript reports high plate-wise concordance).

3. Hit statistics
- Hit count is in the same order of magnitude as manuscript summary when using matching preprocessing and thresholds.
- Thresholds used: `|log2FC| > 1`, `p < 0.05`.

4. Figure regeneration
- Volcano, replicate scatter, histogram, heatmap, Skyline files are generated without missing-data errors.

---

## 9. Known Reproduction Gaps

For strict manuscript-grade reproduction, verify and, where needed, extend:

- exact p-value model implementation details,
- exact control filtering rules for each figure,
- exact visual style, color coding, and label placement,
- secondary validation and GWAS integration steps.

If exact panel-level figure parity is required, treat this as the operational baseline and add manuscript-specific tuning scripts in `prpcscreen/scripts/`.

---

## 10. Troubleshooting

- `Python was not found`:
  - Install Python and re-run setup.

- `ValueError: Expected 384 wells for plate ...`:
  - Ensure one full plate is present after filtering and sorting.

- Empty hits table:
  - Confirm `p_value_log2` and `Mean_log2` exist and are numeric.
  - Check normalization method and control annotations.

- Volcano/Skyline missing points:
  - Verify required columns and non-null values.

---

## 11. Suggested Output Layout

```text
results/
  01_integrated.csv
  02_analyzed.csv
  03_hits.csv
  figures/
    plate_qc_ssmd_controls.png
    plate_well_series_raw_rep1.png
    replicate_agreement_log2fc.png
    distribution_log2fc_rep1_interactive.html
    candidate_volcano_interactive.html
    candidate_flashlight_ranked_meanlog2.png
    plate_heatmap_raw_rep1.png
    grouped_boxplot_raw_rep1.png
    genomic_skyline_meanlog2fc.png
```

---

## 12. Citation and Provenance

Use the manuscript citation and repository URL in reports. Include:
- commit hash used for analysis,
- exact command history,
- software versions (`python --version`, `pip freeze`).

This ensures full computational provenance for manuscript-related reanalysis.

