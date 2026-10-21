# =============================================================================
# 04_sceptre.R
# SCEPTRE: calibrated differential-expression testing of perturbation effects.
# REWRITTEN to start from the shared Seurat object built in 01 via
# sceptre::import_data() (instead of import_data_from_cellranger() on raw 10x),
# so all three methods share exactly the same integrated object.
#
# Input:   results/seurat_object.rds, results/grna_target_map.rds
# Outputs: results/sceptre/ (sceptre_object.rds, discovery results, CO outputs)
# =============================================================================
Sys.setenv(STEP = "04_sceptre")
.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(if (length(.this)) dirname(.this) else "code", "config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(sceptre)
  library(Matrix)
})

OUT <- res_path("sceptre"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(SEURAT_RDS)
grna_target_data_frame <- readRDS(GRNA_MAP_RDS)

## ---- Build SCEPTRE inputs directly from the object ------------------------
response_matrix <- GetAssayData(obj, assay = "RNA",  layer = "counts")   # genes x cells
grna_matrix     <- GetAssayData(obj, assay = "gRNA", layer = "counts")   # guides x cells

# Keep only guides that are actually present in this object's gRNA assay.
grna_target_data_frame <-
  grna_target_data_frame[grna_target_data_frame$grna_id %in% rownames(grna_matrix), ]

log_msg("response: ", nrow(response_matrix), " genes; grna: ",
        nrow(grna_matrix), " guides; cells: ", ncol(response_matrix))

sceptre_object <- import_data(
  response_matrix        = response_matrix,
  grna_matrix            = grna_matrix,
  grna_target_data_frame = grna_target_data_frame,
  moi                    = MOI
)

## ---- Construct pairs -------------------------------------------------------
# Positive controls exist only for targets whose gene symbol is in the response
# matrix; guard so a mismatch doesn't abort the whole run.
positive_control_pairs <- tryCatch(
  construct_positive_control_pairs(sceptre_object),
  error = function(e) { log_msg("no positive-control pairs: ", conditionMessage(e))
                        data.frame() }
)
have_pc <- nrow(positive_control_pairs) > 0

discovery_pairs <- construct_trans_pairs(
  sceptre_object         = sceptre_object,
  positive_control_pairs = positive_control_pairs,
  pairs_to_exclude       = if (have_pc) "pc_pairs" else "none"
)

## ---- Run the pipeline ------------------------------------------------------
parallel_flag <- N_CPU > 1
sceptre_object <- set_analysis_parameters(
  sceptre_object,
  discovery_pairs        = discovery_pairs,
  positive_control_pairs = if (have_pc) positive_control_pairs else NULL
)
sceptre_object <- run_calibration_check(sceptre_object, parallel = parallel_flag)
if (have_pc) sceptre_object <- run_power_check(sceptre_object, parallel = parallel_flag)
sceptre_object <- run_discovery_analysis(sceptre_object, parallel = parallel_flag)

## ---- Collect + annotate results -------------------------------------------
discovery_result <- get_result(sceptre_object, analysis = "run_discovery_analysis")

# Response IDs are gene symbols here; add a symbol column for convenience.
discovery_result$response_gene <- discovery_result$response_id
significant <- subset(discovery_result, significant %in% TRUE)
log_msg("Significant discovery pairs: ", nrow(significant))

## ---- Save ------------------------------------------------------------------
write.csv(discovery_result, file.path(OUT, "discovery_result.csv"), row.names = FALSE)
write.csv(significant,      file.path(OUT, "discovery_significant.csv"), row.names = FALSE)
saveRDS(sceptre_object,     file.path(OUT, "sceptre_object.rds"))
tryCatch(write_outputs_to_directory(sceptre_object, file.path(OUT, "sceptre_outputs")),
         error = function(e) log_msg("write_outputs_to_directory: ", conditionMessage(e)))
log_msg("Wrote SCEPTRE outputs to ", OUT)
