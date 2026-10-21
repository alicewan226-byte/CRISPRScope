# =============================================================================
# config.R  --  shared paths and parameters for the CRISPRScope capsule
# Sourced by 01_build_seurat.R, 02_mixscape.R, 03_music.R, 04_sceptre.R
# =============================================================================
# All scripts read inputs from DATA_DIR and write everything to RESULTS_DIR.
# On Code Ocean these default to /data (read-only) and /results.
# To run locally, override them, e.g.:
#   DATA_DIR=data RESULTS_DIR=results bash code/run
# -----------------------------------------------------------------------------

## ---- Paths ----------------------------------------------------------------
DATA_DIR    <- Sys.getenv("DATA_DIR",    "/data")
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", "/results")

# Input: a Cell Ranger multi/count filtered matrix h5 that contains BOTH
# "Gene Expression" and "CRISPR Guide Capture" feature types.
INPUT_H5 <- file.path(DATA_DIR,
  "SC3_v3_NextGem_DI_CellPlex_CRISPR_A549_30K_A549_Small_Pool_v2_No_Treatment_count_sample_feature_bc_matrix.h5")

# Shared intermediate object built by 01 and consumed by 02/03/04.
SEURAT_RDS   <- file.path(RESULTS_DIR, "seurat_object.rds")   # single-gRNA cells, processed
GRNA_MAP_RDS <- file.path(RESULTS_DIR, "grna_target_map.rds") # grna_id -> grna_target (for SCEPTRE)

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- Biology --------------------------------------------------------------
# This demo dataset is human A549 (10x public A549 CRISPR small pool).
SPECIES       <- "Hs"    # "Mm" (mouse) or "Hs" (human)
CONTROL_LABEL <- "NTC"   # unified label for non-targeting control cells
MOI           <- "low"   # SCEPTRE multiplicity-of-infection: "low" or "high"

# Mitochondrial gene prefix for percent.mt: human "^MT-", mouse "^mt-".
MITO_PATTERN  <- if (SPECIES == "Hs") "^MT-" else "^mt-"

# Hemoglobin genes for the QC red-blood-cell fraction (see 01).
HB_GENES <- if (SPECIES == "Hs")
  c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ") else
  c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Hbb-bh1", "Hbb-y")

## ---- QC thresholds (from the original public_sample6-16.r) -----------------
QC <- list(
  minGene    = 500,
  maxGene    = 9500,
  percentMT  = 25,
  minCount   = 1000,
  percentHB  = 0.05,
  minGuideCount = 5,     # nCount_gRNA > this
  guideDominance = 0.10  # a gRNA is "called" in a cell if it is > this fraction of the cell's gRNA counts
)

## ---- gRNA -> target gene mapping ------------------------------------------
# Target gene is the gRNA feature name with everything from the first "-"
# stripped, e.g. "APH1A-1" -> "APH1A" (or "sgPdcd1-1" -> "sgPdcd1" for
# sg-prefixed libraries). Guides whose name matches the non-targeting pattern
# below are treated as controls instead.
# This SAME mapping is reused by SCEPTRE (04) so all methods agree on targets.
CONTROL_GUIDE_REGEX <- "Non[- ]?Targeting|NO-?SITE|INTERGENIC|^NTC|Safe.?Harbor|scramble|^NT-"

guide_prefix <- function(x) sub("\\-.*", "", x)

# Map a vector of gRNA feature names -> target label. Controls -> control_label.
guide_to_target <- function(guide_names, control_label = CONTROL_LABEL) {
  is_ctrl <- grepl(CONTROL_GUIDE_REGEX, guide_names, ignore.case = TRUE)
  targets <- guide_prefix(guide_names)
  targets[is_ctrl] <- control_label
  targets
}

## ---- Compute --------------------------------------------------------------
N_CPU <- max(1L, parallel::detectCores() - 1L)

## ---- Optional / offline-sensitive steps -----------------------------------
# Cell-cycle scoring in the original script queried Ensembl via biomaRt (online).
# Code Ocean's reproducible run has no network by default, so it is OFF here.
# Turn on only if you enable network access AND want Phase-based diagnostics.
RUN_CELL_CYCLE <- FALSE

## ---- Small helpers --------------------------------------------------------
res_path <- function(...) file.path(RESULTS_DIR, ...)

log_msg <- function(...) cat(sprintf("[%s] %s\n", Sys.getenv("STEP", "capsule"),
                                      paste0(..., collapse = "")))
