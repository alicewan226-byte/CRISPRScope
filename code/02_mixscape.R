# =============================================================================
# 02_mixscape.R
# MIXSCAPE: identify cells with a detectable perturbation, starting from the
# shared Seurat object built in 01.
#
# Adapted from public_sample6-16.r (CalcPerturbSig / RunMixscape section).
# Input:   results/seurat_object.rds
# Outputs: results/mixscape/Mixscape.rds + diagnostic plots
# =============================================================================
Sys.setenv(STEP = "02_mixscape")
.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(if (length(.this)) dirname(.this) else "code", "config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

OUT <- res_path("mixscape"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
save_plot <- function(p, file, w = 12, h = 8) ggsave(file.path(OUT, file), p, width = w, height = h)

## ---- bound normalmixEM so a hard-to-fit gene can't spin --------------------
# RunMixscape fits a 2-component mixture per target gene with
# mixtools::normalmixEM, whose random-restart loop (defaults maxrestarts = 20,
# maxit = 1000) can spin for many minutes on a degenerate gene before it either
# converges, errors, or is abandoned. Cap the restarts/iterations so any gene
# that won't fit fails FAST -- the per-gene tryCatch below then skips it (its
# cells stay NP) instead of the whole run hanging. Well-behaved genes are
# unaffected (they converge in a handful of iterations).
loadNamespace("mixtools")
local({
  ns   <- asNamespace("mixtools")
  orig <- get("normalmixEM", ns)               # keep a handle to the real fn
  safe <- function(...) {                       # closes over `orig` (this env)
    a <- list(...); a$maxit <- 300L; a$maxrestarts <- 2L
    do.call(orig, a)
  }
  assignInNamespace("normalmixEM", safe, ns)
})

obj <- readRDS(SEURAT_RDS)
DefaultAssay(obj) <- "RNA"
log_msg("Loaded ", ncol(obj), " single-gRNA cells")

## ---- 1. Re-embed the single-gRNA population -------------------------------
obj <- NormalizeData(obj) %>% FindVariableFeatures() %>% ScaleData()
obj <- RunPCA(obj)
obj <- RunUMAP(obj, dims = 1:40)

obj$crispr <- ifelse(obj$gene == CONTROL_LABEL, "NT", "Perturbed")

## ---- 2. Perturbation signature (removes confounding variation) ------------
obj <- CalcPerturbSig(
  object        = obj,
  assay         = "RNA",
  slot          = "data",
  gd.class      = "gene",
  nt.cell.class = CONTROL_LABEL,
  reduction     = "pca",
  ndims         = 40,
  num.neighbors = 5,
  new.assay.name = "PRTB"
)

DefaultAssay(obj) <- "PRTB"
VariableFeatures(obj) <- VariableFeatures(obj[["RNA"]])
obj <- ScaleData(obj, do.scale = FALSE, do.center = TRUE)
obj <- RunPCA(obj, reduction.key = "prtbpca", reduction.name = "prtbpca")
obj <- RunUMAP(obj, dims = 1:40, reduction = "prtbpca",
               reduction.key = "prtbumap", reduction.name = "prtbumap")

save_plot(DimPlot(obj, group.by = "crispr", reduction = "prtbumap",
                  split.by = "crispr", ncol = 1) + ggtitle("Perturbation status (PRTB UMAP)"),
          "prtb_umap.png")

## ---- 3. Run Mixscape (resilient: per-gene, skip genes whose EM fails) ------
# Seurat's RunMixscape fits a 2-component mixture (mixtools::normalmixEM) for
# each target gene inside a single call. On this dataset a few low-cell targets
# produce a degenerate fit ("NA/NaN/Inf in foreign function call") that aborts
# the WHOLE step. Run each target against the NT cells on its own so one bad
# gene can't sink the rest; failures are logged and their cells left as NP.
# normalmixEM is randomly initialised: seed it identically per gene. This is
# both for reproducibility and for stability -- letting the RNG state drift
# across genes (seeding once before the loop) can push a fit into a
# pathological, non-converging spin. Seeding per gene keeps every gene fast.
mix_params <- list(assay = "PRTB", slot = "scale.data", labels = "gene",
                   nt.class.name = CONTROL_LABEL, min.de.genes = 1, min.cells = 30,
                   logfc.threshold = 0.25, iter.num = 10, de.assay = "RNA",
                   verbose = FALSE, prtb.type = "KO")

targets <- setdiff(unique(as.character(obj$gene)), CONTROL_LABEL)
mix_class        <- setNames(rep(NA_character_, ncol(obj)), colnames(obj))
mix_class_global <- setNames(rep(NA_character_, ncol(obj)), colnames(obj))
mix_p_ko         <- setNames(rep(NA_real_,      ncol(obj)), colnames(obj))
tools_list <- list()

for (g in targets) {
  sub <- subset(obj, cells = colnames(obj)[obj$gene %in% c(CONTROL_LABEL, g)])
  ok <- tryCatch({
    set.seed(1)                        # identical mixture-model init for every gene
    sub <- do.call(RunMixscape, c(list(object = sub), mix_params))
    TRUE
  }, error = function(e) {
    log_msg("mixscape skip ", g, ": ", conditionMessage(e)); FALSE })
  if (!ok) next
  gcells <- colnames(sub)[sub$gene == g]      # record this gene's cells only
  mix_class[gcells]        <- as.character(sub$mixscape_class[gcells])
  mix_class_global[gcells] <- as.character(sub$mixscape_class.global[gcells])
  if (!is.null(sub$mixscape_class_p_ko)) mix_p_ko[gcells] <- sub$mixscape_class_p_ko[gcells]
  if (!is.null(sub@tools$RunMixscape[[g]])) tools_list[[g]] <- sub@tools$RunMixscape[[g]]
}

# NT cells keep the control label; untested/failed target cells become NP.
mix_class[obj$gene == CONTROL_LABEL]        <- CONTROL_LABEL
mix_class_global[obj$gene == CONTROL_LABEL] <- CONTROL_LABEL
mix_class[is.na(mix_class)]               <- "NP"
mix_class_global[is.na(mix_class_global)] <- "NP"

obj$mixscape_class        <- mix_class
obj$mixscape_class.global <- mix_class_global
obj$mixscape_class_p_ko   <- mix_p_ko
obj@tools$RunMixscape     <- tools_list
log_msg("Mixscape tested ", length(tools_list), "/", length(targets), " targets")

log_msg("mixscape_class.global table:")
print(table(obj$mixscape_class.global))

## ---- 4. % KO per target gene ----------------------------------------------
df  <- prop.table(table(obj$mixscape_class.global, obj$gene), 2)
df2 <- reshape2::melt(df)
colnames(df2) <- c("class", "gene", "value")
df2$class <- factor(df2$class, levels = c(CONTROL_LABEL, "NP", "KO"))
df3 <- df2[df2$gene != CONTROL_LABEL, ]

p_ko <- ggplot(df3, aes(x = gene, y = value * 100, fill = class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("grey49", "grey79", "coral1")) +
  labs(x = "target gene", y = "% of cells", fill = "mixscape class") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_ko, "percent_KO_per_gene.png", w = 10, h = 6)

## ---- 5. Per-gene perturbation-score plots ---------------------------------
for (g in setdiff(names(obj@tools$RunMixscape), CONTROL_LABEL)) {
  tryCatch({
    pp <- PlotPerturbScore(obj, target.gene.ident = g,
                           mixscape.class = "mixscape_class", col = "coral2") +
      labs(fill = "mixscape class")
    save_plot(pp, paste0("perturbScore_", g, ".png"))
  }, error = function(e) log_msg("perturbScore skip ", g, ": ", conditionMessage(e)))
}

## ---- 6. DE heatmap per gene (KO vs control) -------------------------------
Idents(obj) <- "gene"
for (g in setdiff(unique(as.character(obj$gene)), c(CONTROL_LABEL, "NP"))) {
  tryCatch({
    hm <- MixscapeHeatmap(obj, ident.1 = CONTROL_LABEL, ident.2 = g,
                          logfc.threshold = 0.25, pval.cutoff = 0.05, balanced = TRUE,
                          assay = "RNA", max.genes = 10, angle = 0,
                          group.by = "mixscape_class", max.cells.group = 300, size = 6.5)
    save_plot(hm, paste0("heatmap_", g, ".png"))
  }, error = function(e) log_msg("heatmap skip ", g, ": ", conditionMessage(e)))
}

saveRDS(obj, file.path(OUT, "Mixscape.rds"))
log_msg("Wrote ", file.path(OUT, "Mixscape.rds"))
