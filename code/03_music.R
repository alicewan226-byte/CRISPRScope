# =============================================================================
# 03_music.R
# MUSIC: topic-model based perturbation-effect ranking, starting from the
# shared Seurat object built in 01 (no separate so.rds).
#
# Adapted from MUSIC.r.
# Input:   results/seurat_object.rds
# Outputs: results/music/music.rds (all intermediate objects) + plots
# =============================================================================
Sys.setenv(STEP = "03_music")
.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(if (length(.this)) dirname(.this) else "code", "config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(MUSIC)
})

OUT <- res_path("music"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(SEURAT_RDS)

## ---- Build MUSIC inputs from the shared object ----------------------------
# Drop any non-perturbed cells (defensive; the shared object is single-gRNA
# already) and relabel the non-targeting control to MUSIC's "CTRL" convention.
obj <- subset(obj, subset = gene != "NP")
obj$gene <- ifelse(obj$gene == CONTROL_LABEL, "CTRL", obj$gene)

expression_profile  <- as.data.frame(GetAssayData(obj, assay = "RNA", layer = "counts"))
perturb_information  <- as.character(obj$gene)
names(perturb_information) <- colnames(obj)
log_msg("Expression profile: ", nrow(expression_profile), " genes x ",
        ncol(expression_profile), " cells")
rm(obj); gc()

## ---- MUSIC pipeline --------------------------------------------------------
crop_seq_list <- Input_preprocess(expression_profile, perturb_information)

crop_seq_qc <- Cell_qc(crop_seq_list$expression_profile,
                       crop_seq_list$perturb_information,
                       species = SPECIES, plot = FALSE)

set.seed(234)
crop_seq_imputation <- Data_imputation(crop_seq_qc$expression_profile,
                                        crop_seq_qc$perturb_information,
                                        cpu_num = N_CPU)

crop_seq_filtered <- Cell_filtering(crop_seq_imputation$expression_profile,
                                    crop_seq_imputation$perturb_information,
                                    cell_num_threshold = 20,
                                    plot = FALSE, cpu_num = N_CPU)

crop_seq_vargene <- Get_high_varGenes(crop_seq_filtered$expression_profile,
                                      crop_seq_filtered$perturb_information,
                                      plot = FALSE)

## Topic models over a small range of topic numbers.
topic_model_list <- Get_topics(crop_seq_vargene$expression_profile,
                               crop_seq_vargene$perturb_information,
                               topic_number = c(4:6))

optimalModel <- Select_topic_number(topic_model_list$models, plot = FALSE)

topic_func   <- Topic_func_anno(optimalModel, species = SPECIES, plot = FALSE)
distri_diff  <- Diff_topic_distri(optimalModel, topic_model_list$perturb_information,
                                  plot = FALSE)

rank_overall_result       <- Rank_overall(distri_diff)
rank_topic_specific_result <- Rank_specific(distri_diff)
perturb_cor               <- Correlation_perturbation(distri_diff, plot = FALSE)

## ---- Save ------------------------------------------------------------------
write.csv(rank_overall_result,        file.path(OUT, "rank_overall.csv"),        row.names = FALSE)
write.csv(rank_topic_specific_result, file.path(OUT, "rank_topic_specific.csv"), row.names = FALSE)
saveRDS(mget(ls()), file.path(OUT, "music.rds"))
log_msg("Wrote ", file.path(OUT, "music.rds"))
