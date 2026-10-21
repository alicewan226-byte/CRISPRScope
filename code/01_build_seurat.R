# =============================================================================
# 01_build_seurat.R
# Build ONE integrated Seurat object from a Cell Ranger CRISPR h5, and the
# gRNA -> target map. This object is the single shared input for all three
# perturbation methods (MIXSCAPE / MUSIC / SCEPTRE).
#
# Adapted from the original public_sample6-16.r (exploratory script).
# Outputs:  results/seurat_object.rds, results/grna_target_map.rds
# =============================================================================
Sys.setenv(STEP = "01_build")
# Locate config.R relative to this script, regardless of working directory.
.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(if (length(.this)) dirname(.this) else "code", "config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
})

set.seed(1)   # reproducible PCA / UMAP / clustering

stopifnot("Input h5 not found (see config.R INPUT_H5)" = file.exists(INPUT_H5))
log_msg("Reading ", INPUT_H5)

## ---- 1. Load 10x data, build RNA + gRNA assays ----------------------------
raw <- Read10X_h5(INPUT_H5)
obj <- CreateSeuratObject(counts = raw[["Gene Expression"]], assay = "RNA")
obj[["gRNA"]] <- CreateAssayObject(counts = raw[["CRISPR Guide Capture"]])

## ---- 2. QC ----------------------------------------------------------------
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = MITO_PATTERN)

hb_present <- intersect(HB_GENES, rownames(obj[["RNA"]]))
obj[["percent.HB"]] <- if (length(hb_present))
  PercentageFeatureSet(obj, features = hb_present) else 0

log_msg("Cells before QC: ", ncol(obj))
obj <- subset(
  obj,
  subset = nFeature_RNA > QC$minGene & nFeature_RNA < QC$maxGene &
           percent.mt  < QC$percentMT & nCount_RNA   > QC$minCount &
           nCount_gRNA > QC$minGuideCount & percent.HB < QC$percentHB
)
log_msg("Cells after QC:  ", ncol(obj))

## ---- 3. Single-gRNA calling -----------------------------------------------
# Keep only gRNA counts that make up > guideDominance of a cell's total gRNA
# signal; from that, count how many guides each cell carries.
gcounts <- as.matrix(GetAssayData(obj, assay = "gRNA", layer = "counts"))
col_tot <- colSums(gcounts)
frac    <- sweep(gcounts, 2, ifelse(col_tot == 0, 1, col_tot), "/")
kept    <- gcounts * (frac > QC$guideDominance)   # zero out minor guides
kept[is.na(kept)] <- 0

obj[["new_gRNA"]] <- CreateAssayObject(counts = kept)

n_guides <- colSums(kept > 0)
obj$gRNA_perturbed <- dplyr::case_when(
  n_guides == 1 ~ "single",
  n_guides >  1 ~ "multiple",
  TRUE          ~ "non-perturbed"
)

## ---- 4. Assign dominant guide (gRNA_ID) and target gene --------------------
guide_names <- rownames(gcounts)
dominant_idx <- apply(gcounts, 2, which.max)          # dominant guide per cell
dominant_id  <- guide_names[dominant_idx]

obj$gRNA_ID <- ifelse(obj$gRNA_perturbed == "non-perturbed", "NP", dominant_id)
obj$gene    <- ifelse(obj$gRNA_perturbed == "non-perturbed",
                      "NP", guide_to_target(dominant_id))   # controls -> CONTROL_LABEL

log_msg("gRNA_perturbed table:")
print(table(obj$gRNA_perturbed))
log_msg("gene (target) table:")
print(table(obj$gene))

## ---- 5. Standard processing (needed for MIXSCAPE PRTB / UMAP) --------------
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:50)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:50)
obj <- FindClusters(obj, resolution = 0.5)

## ---- 6. Optional cell-cycle scoring (OFF by default; needs network) --------
if (RUN_CELL_CYCLE) {
  log_msg("Cell-cycle scoring enabled (requires biomaRt / network access)")
  library(biomaRt)
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  to_mouse <- function(g) {
    m <- getBM(attributes = c("external_gene_name",
                              "mmusculus_homolog_associated_gene_name"),
               filters = "external_gene_name", values = g, mart = ensembl)
    m <- m$mmusculus_homolog_associated_gene_name
    m[!is.na(m) & m != ""]
  }
  s.m   <- CaseMatch(to_mouse(cc.genes$s.genes),   rownames(obj))
  g2m.m <- CaseMatch(to_mouse(cc.genes$g2m.genes), rownames(obj))
  obj   <- CellCycleScoring(obj, s.features = s.m, g2m.features = g2m.m)
} else {
  log_msg("Cell-cycle scoring skipped (RUN_CELL_CYCLE = FALSE)")
}

## ---- 7. Keep single-gRNA cells (the population all 3 methods analyse) ------
obj_single <- subset(obj, subset = gRNA_perturbed == "single")
log_msg("Single-gRNA cells retained: ", ncol(obj_single))

## ---- 8. gRNA -> target map for SCEPTRE ------------------------------------
# SCEPTRE pairs a guide's target with that gene's expression, so the target
# label must match a gene symbol in the RNA assay. Strip the leading "sg"
# (e.g. "sgPdcd1" -> "Pdcd1"); controls become "non-targeting".
is_ctrl_guide <- grepl(CONTROL_GUIDE_REGEX, guide_names, ignore.case = TRUE)
grna_target_data_frame <- data.frame(
  grna_id     = guide_names,
  grna_target = ifelse(is_ctrl_guide, "non-targeting",
                       sub("^sg", "", guide_prefix(guide_names))),
  stringsAsFactors = FALSE
)

## ---- 9. Persist shared inputs ---------------------------------------------
saveRDS(obj_single, SEURAT_RDS)
saveRDS(grna_target_data_frame, GRNA_MAP_RDS)
log_msg("Wrote ", SEURAT_RDS)
log_msg("Wrote ", GRNA_MAP_RDS)
