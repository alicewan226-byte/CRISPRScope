# CRISPRScope: A Comprehensive Platform for Single-Cell CRISPR Data Analysis and Visualization

CRISPRScope is a flexible pipeline for single-cell CRISPR (Perturb-seq) data:
raw-data processing, quality control, perturbation modeling, and visualization.
This capsule integrates the data into **one Seurat object** and then analyses
perturbation effects with three complementary methods — **MIXSCAPE**, **MUSIC**
and **SCEPTRE** — all driven from that same object.

---

## Capsule layout

```
code/
  config.R            shared paths + parameters (edit this to retune)
  01_build_seurat.R   build the integrated Seurat object (+ gRNA→target map)
  02_mixscape.R       MIXSCAPE: per-cell knockout detection
  03_music.R          MUSIC: topic-model perturbation ranking
  04_sceptre.R        SCEPTRE: calibrated DE testing (import_data from the object)
  run                 entry point: runs 01 → 02 → 03 → 04
data/                 input (mounted read-only on Code Ocean)
results/              all outputs are written here
environment/Dockerfile
metadata/metadata.yml
original_scripts/     the pre-capsule exploratory scripts, kept for reference
```

## Input

Place a Cell Ranger *filtered* feature-barcode matrix in `data/`. The demo
dataset is the public 10x **human A549** CRISPR small pool:

```
data/SC3_v3_NextGem_DI_CellPlex_CRISPR_A549_30K_A549_Small_Pool_v2_No_Treatment_count_sample_feature_bc_matrix.h5
```

The h5 must contain **both** feature types: `Gene Expression` and
`CRISPR Guide Capture` (a `Multiplexing Capture` / CMO block, if present, is
ignored). gRNA feature names are expected as `<Gene>-<n>` (e.g. `APH1A-1`) or
`sg<Gene>-<n>`; non-targeting guides are auto-detected by name
(`Non-Targeting`, `NO-SITE`, `INTERGENIC`, …). Adjust `CONTROL_GUIDE_REGEX`
and `INPUT_H5` in `code/config.R` if your data differs.

## How to run

On Code Ocean, click **Reproducible Run** (executes `code/run`).

Locally:

```bash
DATA_DIR=data RESULTS_DIR=results bash code/run
```

## Outputs (`results/`)

| File / dir | Produced by | Contents |
|---|---|---|
| `seurat_object.rds` | 01 | integrated, QC'd, single-gRNA Seurat object |
| `grna_target_map.rds` | 01 | gRNA → target gene table (used by SCEPTRE) |
| `mixscape/` | 02 | `Mixscape.rds`, %KO-per-gene, perturb-score & DE-heatmap plots |
| `music/` | 03 | `music.rds`, overall & topic-specific perturbation rankings |
| `sceptre/` | 04 | `sceptre_object.rds`, discovery results, calibration outputs |

## Key parameters (`code/config.R`)

- `SPECIES` – `"Hs"` (human; this demo) or `"Mm"`. Drives the mitochondrial
  prefix (`^MT-` / `^mt-`) and hemoglobin gene list.
- `CONTROL_LABEL` – unified non-targeting label (`"NTC"`).
- `MOI` – SCEPTRE multiplicity of infection (`"low"` / `"high"`).
- `QC` – gene/count/mito/hemoglobin thresholds.
- `N_CPU` – parallel cores (defaults to cores − 1).
- `RUN_CELL_CYCLE` – `FALSE` by default. The cell-cycle step queries Ensembl
  via **biomaRt (online)**; enable only if the run has network access.

## Environment

R 4.2 (Bioconductor 3.16 base). **Seurat is pinned to 4.x** — the mixscape
functions used here assume the v4 assay model; do not upgrade to v5 without
re-testing. See `environment/Dockerfile` for the full pinned package list
(Seurat, MUSIC, SAVER, sceptre, ComplexHeatmap, clusterProfiler, org.Hs.eg.db …).

## Method notes

- **MIXSCAPE** operates directly on the Seurat object, computes a perturbation
  signature (PRTB) and classifies each cell as knockout (KO) / non-perturbed (NP).
- **MUSIC** takes the RNA count matrix + per-cell target from the object,
  imputes (SAVER) and ranks perturbations via topic models. The imputation step
  is compute-heavy; try a subset first.
- **SCEPTRE** is fed from the object with `import_data()` (RNA counts as the
  response matrix, raw gRNA counts as the guide matrix, and the shared
  gRNA→target map), replacing the original `import_data_from_cellranger()`.
