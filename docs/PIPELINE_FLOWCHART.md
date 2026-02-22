# PrPC Screen Pipeline Flowchart

This document summarizes the end-to-end data flow implemented in this repository.

Key correction to common assumptions:
- Plate/well-to-gene assignment comes from the `Layout CSV`.
- The genomics workbook (`Genomics XLSX`) is used for:
  - Skyline genomic plotting input, and
  - optional sublibrary mapping in interactive volcano output.

## High-Level Flow

```mermaid
flowchart LR
    A[Raw assay folder<br/>TR-FRET files + optional GLO files]
    B[Layout CSV<br/>Plate_number_384, Well_number_384,<br/>Is_NT_ctrl, Is_pos_ctrl, Entrez_ID, ...]
    C[Genomics XLSX<br/>Skyline sheet with Gene_symbol, Mean_log2FC,<br/>Chromosome, Start_Position<br/>+ optional Sublibrary]

    S1[1 merge_assay_exports.py<br/>Discover files recursively<br/>Parse 16x24 matrices<br/>Apply 2-channel TR-FRET correction when available<br/>Infer replicate A/B and plate labels<br/>Align to layout order]
    O1[results/01_integrated.csv<br/>Layout rows + Raw_rep1 + Raw_rep2<br/>+ optional CellTiterGlo_raw]

    S2[2 compute_screen_metrics.py<br/>Clip raw >= 1<br/>Per-plate normalization<br/>Compute DeltaNT/FoldNT/Log2FC/PercActivation<br/>Compute SSMD, moderated SSMD, p-values<br/>Compute Mean_*]
    O2[results/02_analyzed.csv<br/>Integrated + all derived metrics/statistics]
    O3[results/03_hits.csv<br/>Subset where p_value_log2 < 0.05 and abs Mean_log2 > 1]

    S3[3 plot_plate_health.py<br/>Plate-level SSMD controls QC]
    F1[plate_qc_ssmd_controls.png]

    S4[4 plot_well_trajectories.py<br/>Well-order trend plot]
    F2[plate_well_series_raw_rep1.png]

    S5[5 plot_replicate_agreement.py<br/>Scatter + Bland-Altman + error-vs-effect + binned corr]
    F3[replicate_agreement_log2fc.png]

    S6[6 plot_signal_distributions.py<br/>Interactive category histogram<br/>Genes vs Non-targeting vs Positive]
    F4[distribution_log2fc_rep1_interactive.html<br/>+ .html.gz]

    S7[7 plot_candidate_landscape.py<br/>Interactive volcano + flashlight ranking<br/>Optional sublibrary mapping from genomics XLSX<br/>Fallback map CSV if missing]
    F5[candidate_volcano_interactive.html<br/>+ .html.gz]
    F6[candidate_flashlight_ranked_meanlog2.png]

    S8[8 plot_spatial_and_group_views.py<br/>Per-plate 384 heatmap Raw_rep1<br/>Violin+box grouped distribution]
    F7[plate_heatmap_raw_rep1.png]
    F8[grouped_boxplot_raw_rep1.png]

    S9[9 plot_genomic_signal_skyline.py<br/>Read genomics sheet<br/>Filter valid chrom/position/log2FC<br/>Plot genome-position skyline]
    F9[genomic_skyline_meanlog2fc.png]

    A --> S1
    B --> S1
    S1 --> O1
    O1 --> S2
    S2 --> O2
    S2 --> O3
    O2 --> S3 --> F1
    O2 --> S4 --> F2
    O2 --> S5 --> F3
    O2 --> S6 --> F4
    O2 --> S7 --> F5
    O2 --> S7 --> F6
    O2 --> S8 --> F7
    O2 --> S8 --> F8
    C --> S7
    C --> S9
    S9 --> F9

    classDef scriptNode color:#ff0000,stroke:#ff0000,stroke-width:2px,fill:#fff5f5;
    classDef fileNode color:#0b2e6b,stroke:#0b2e6b,stroke-width:1.5px,fill:#f5f8ff;
    class S1,S2,S3,S4,S5,S6,S7,S8,S9 scriptNode;
    class A,B,C,O1,O2,O3,F1,F2,F3,F4,F5,F6,F7,F8,F9 fileNode;
```

## Pooled Flow (Reusing Existing Figure Logic)

For pooled screens, the upstream analysis is different, but the downstream figure stack is reused by emitting a compatible analyzed table.

```mermaid
flowchart LR
    P0[Pooled guide table<br/>CSV/TSV/XLSX with replicate columns<br/>e.g. Negative_R* / Positive_R*]
    P1[compute_pooled_metrics.py<br/>Detect/validate replicate columns<br/>Size-factor normalization<br/>Compute Log2FC_rep*, Mean_log2, p_value_log2, FDR<br/>Infer controls and call hits]
    PI[results_pooled/01_integrated.csv]
    PA[results_pooled/02_analyzed.csv<br/>Figure-compatible schema]
    PH[results_pooled/03_hits.csv]

    V1[plot_replicate_agreement.py<br/>--stem Log2FC]
    V2[plot_signal_distributions.py<br/>--column Log2FC_rep1]
    V3[plot_candidate_landscape.py<br/>Volcano + flashlight]
    V4[plot_well_trajectories.py<br/>Feature-order series]
    V5[plot_genomic_signal_skyline.py<br/>optional, genomics workbook]

    P0 --> P1
    P1 --> PI
    P1 --> PA
    P1 --> PH
    PA --> V1
    PA --> V2
    PA --> V3
    PA --> V4
    PA --> V5

    classDef scriptNode color:#ff0000,stroke:#ff0000,stroke-width:2px,fill:#fff5f5;
    classDef fileNode color:#0b2e6b,stroke:#0b2e6b,stroke-width:1.5px,fill:#f5f8ff;
    class P1,V1,V2,V3,V4,V5 scriptNode;
    class P0,PI,PA,PH fileNode;
```

## Stage Contracts

## Stage 1: Integration (`merge_assay_exports.py`)

Inputs:
- Raw assay directory/file root (`*.csv`, `*.tsv`, `*.txt`)
- Layout CSV (annotation scaffold)
- Optional GLO files

Transformations:
- Recursive discovery of measurement files.
- TR-FRET plate parsing from A..P x 1..24 tables.
- Two-channel TR-FRET correction when full channel blocks exist.
- Replicate classification from filename suffix (`A` -> rep1, `B` -> rep2).
- Plate label inference and alignment to layout plate order.
- Stale analysis columns removed from layout scaffold before attaching new signals.

Output:
- `results/01_integrated.csv`
- Contains layout metadata plus `Raw_rep1`, `Raw_rep2`, and optional `CellTiterGlo_raw`.

## Stage 2: Metrics and Hits (`compute_screen_metrics.py`)

Input:
- `results/01_integrated.csv`

Transformations:
- Numeric coercion and floor clipping for `Raw_rep1`/`Raw_rep2` to minimum 1.
- Per-plate normalization (`norm_plates`) to derive:
  - `DeltaNT_rep*`
  - `FoldNT_rep*`
  - `PercActivation_rep*`
  - `Raw_log2_rep*`
  - `Log2FC_rep*`
- GLO-adjusted variants when GLO is present:
  - `Raw_Glo_rep*`, `DeltaNT_Glo_rep*`, `FoldNT_Glo_rep*`, `Log2FC_Glo_rep*`
- Statistics per metric pair:
  - `SSMD_*`, `SSMD_mod_*`, `Mean_*`, `p_value_*`, `p_value_repro_*`
- Hit calling:
  - `p_value_log2 < 0.05` and `abs(Mean_log2) > 1.0`

Outputs:
- `results/02_analyzed.csv` (full derived table)
- `results/03_hits.csv` (hit subset of analyzed rows)

## Stages 3-9: Figures

Shared primary input:
- `results/02_analyzed.csv`

Generated outputs:
- `results/figures/plate_qc_ssmd_controls.png`
- `results/figures/plate_well_series_raw_rep1.png`
- `results/figures/replicate_agreement_log2fc.png`
- `results/figures/distribution_log2fc_rep1_interactive.html` (+ `.gz`)
- `results/figures/candidate_volcano_interactive.html` (+ `.gz`)
- `results/figures/candidate_flashlight_ranked_meanlog2.png`
- `results/figures/plate_heatmap_raw_rep1.png`
- `results/figures/grouped_boxplot_raw_rep1.png`
- `results/figures/genomic_skyline_meanlog2fc.png`

Genomics workbook usage:
- Stage 7 (`plot_candidate_landscape.py`):
  - Optional sublibrary enrichment labels for volcano filtering.
  - If `Sublibrary` is not found, fallback to `prpcscreen/misc/supplementary_sublibrary_map.csv`.
- Stage 9 (`plot_genomic_signal_skyline.py`):
  - Required source sheet containing:
    - `Gene_symbol`
    - `Mean_log2FC`
    - `Chromosome`
    - `Start_Position`
  - `--sheet` defaults to `skylineplot2`, with fallback to a matching/compatible sheet; if none is compatible, the step fails with a schema error.

## Optional Utility (not part of 9-stage core flow)

- `prpcscreen/scripts/remap_plate_coordinates.py`
  - Converts `Well_number_384` to `Well_number_96` mapping for downstream compatibility.
