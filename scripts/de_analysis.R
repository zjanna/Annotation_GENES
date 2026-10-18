#!/usr/bin/env Rscript
# ============================================================================
# de_analysis.R: differential expression per annotation + cross-annotation
#                 comparison with edgeR (TMM) + limma (voom), as in the pig anchor study.
# ----------------------------------------------------------------------------
# Usage:  Rscript de_analysis.R <organism>          e.g.  Rscript de_analysis.R pig
#
# Inputs (produced by the Snakefiles + one small file you provide):
#   results/<organism>/counts/<ann>/all_samples.counts.tsv   (one per annotation)
#   reference/<organism>/annotations/<ann>.gtf               (for gene_id -> symbol)
#   results/<organism>/targets.tsv                           (YOU provide)
#       tab-separated with a header, columns:
#         sample     condition   [block]
#       - sample    : must match the count-matrix column names
#       - condition : the grouping variable (e.g. control / hernia)
#       - block     : OPTIONAL random effect (e.g. litter/animal) -> triggers
#                     limma duplicateCorrelation, as in the pig study
#       - sex       : OPTIONAL covariate added to the design (pig study: ~ trt + sex)
#
# Outputs (results/<organism>/de/):
#   <ann>_toptable.tsv         full limma table per annotation
#   <ann>_DEG.tsv              significant genes per annotation
#   comparison_jaccard.tsv     pairwise Jaccard of DEG sets (by gene symbol)
#   comparison_vs_<anchor>.tsv gained/lost vs the Ensembl anchor
#   comparison_logFC_cor.tsv   Spearman correlation of logFC on shared genes
# ============================================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

# ---- parameters ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript de_analysis.R <organism>")
org      <- args[1]
ANCHOR   <- "ensembl"     # the annotation to compare the others against
FDR      <- 0.05          # significance threshold (adjusted p-value)
# |log2FC| threshold: pass as 2nd arg to match each anchor paper.
# e.g. mouse-limma: 1 ; PIG (Wozniak): 0 (no fold-change filter, adj.P<0.05 only)
LFC      <- if (length(args) >= 2) as.numeric(args[2]) else 1
# Low-expression filter, 3rd arg: "filterByExpr" (default) or "published" =
# CPM > 0.5 in at least two samples, as described by Wozniak et al. 2025 (pig)
FILTER   <- if (length(args) >= 3) args[3] else "filterByExpr"

annot_dir <- file.path("reference", org, "annotations")
count_of  <- function(a) file.path("results", org, "counts", a, "all_samples.counts.tsv")

# Only analyse annotations for which a count matrix actually exists.
candidate_ann <- c("ensembl", "refseq", "gencode", "ucsc")
ANN <- candidate_ann[file.exists(vapply(candidate_ann, count_of, ""))]
message("Annotations found: ", paste(ANN, collapse = ", "))

out_dir <- file.path("results", org, "de")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- sample sheet ----
targets <- read.delim(file.path("results", org, "targets.tsv"),
                       stringsAsFactors = FALSE)
stopifnot(all(c("sample", "condition") %in% colnames(targets)))
has_block <- "block" %in% colnames(targets)
has_sex   <- "sex"   %in% colnames(targets)

# ---- helper: gene_id -> gene_name (symbol) map from a GTF ----
# Each annotation names genes differently (Ensembl IDs, RefSeq symbols, ...),
# so to compare DEG SETS across annotations we translate every gene_id to its
# gene symbol and match on that. (A stricter coordinate-level crosswalk is the
# rigorous next step; symbol matching is the pragmatic first pass.)
gene_symbol_map <- function(a) {
  # Prefer a small precomputed 2-column map (gene_id <TAB> gene_name) if present,
  # else parse the full GTF, else fall back to gene_id (caller handles NAs).
  mapf <- file.path(annot_dir, paste0(a, ".genemap.tsv"))
  gtff <- file.path(annot_dir, paste0(a, ".gtf"))
  if (file.exists(mapf)) {
    m <- read.delim(mapf, header = FALSE, stringsAsFactors = FALSE)
    return(setNames(toupper(m[[2]]), m[[1]]))
  }
  if (file.exists(gtff)) {
    gtf <- readLines(gtff)
    gtf <- gtf[!grepl("^#", gtf)]
    gtf <- gtf[grepl("\tgene\t", gtf)]            # gene-level lines only
    gid <- sub('.*gene_id "([^"]+)".*',   "\\1", gtf)
    gnm <- ifelse(grepl('gene_name "', gtf),
                  sub('.*gene_name "([^"]+)".*', "\\1", gtf), gid)
    return(setNames(toupper(gnm), gid))
  }
  message("  note: no genemap/GTF for ", a, " -> matching by gene_id only")
  setNames(character(0), character(0))
}

# ---- differential expression for one annotation ----
run_de <- function(a) {
  m <- read.delim(count_of(a), check.names = FALSE, stringsAsFactors = FALSE)
  rownames(m) <- m$Geneid
  m$Geneid <- NULL
  m <- m[, targets$sample, drop = FALSE]          # align/order columns to targets

  grp    <- factor(targets$condition)
  # cell-means design with sex as covariate + explicit UH-vs-control contrast
  sex    <- if (has_sex) factor(targets$sex) else NULL
  design <- if (has_sex) model.matrix(~ 0 + grp + sex) else model.matrix(~ 0 + grp)
  cmat   <- matrix(0, ncol(design), 1, dimnames = list(colnames(design), "T_vs_ref"))
  cmat[paste0("grp", levels(grp)[2]), 1] <-  1          # 2nd level (e.g. UH) ...
  cmat[paste0("grp", levels(grp)[1]), 1] <- -1          # ... minus reference (e.g. C)

  y <- DGEList(counts = as.matrix(m), group = grp)
  keep <- if (FILTER == "published") rowSums(cpm(y) > 0.5) >= 2 else filterByExpr(y, group = grp)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")         # edgeR TMM normalisation
  v <- voom(y, design)                            # limma-voom

  if (has_block) {                                # duplicateCorrelation (pig study)
    corfit <- duplicateCorrelation(v, design, block = targets$block)
    # second voom pass with the block correlation (standard limma workflow)
    v <- voom(y, design, block = targets$block, correlation = corfit$consensus)
    fit <- lmFit(v, design, block = targets$block, correlation = corfit$consensus)
  } else {
    fit <- lmFit(v, design)
  }
  fit <- eBayes(contrasts.fit(fit, cmat))
  tt <- topTable(fit, coef = 1, number = Inf, sort.by = "none")
  tt$gene_id <- rownames(tt)
  tt$symbol  <- gene_symbol_map(a)[tt$gene_id]
  tt$symbol[is.na(tt$symbol)] <- tt$gene_id[is.na(tt$symbol)]
  write.table(tt, file.path(out_dir, paste0(a, "_toptable.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  deg <- tt[!is.na(tt$adj.P.Val) & tt$adj.P.Val < FDR & abs(tt$logFC) >= LFC, ]
  write.table(deg, file.path(out_dir, paste0(a, "_DEG.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("  %-8s : %d genes tested, %d DEGs", a, nrow(tt), nrow(deg)))
  list(tt = tt, deg_symbols = unique(deg$symbol),
       lfc = setNames(tt$logFC, tt$symbol))
}

message("Running DE per annotation ...")
res <- setNames(lapply(ANN, run_de), ANN)

# ---- Layer 1a: pairwise Jaccard of DEG sets (by symbol) ----
jac <- function(x, y) if (length(union(x, y)) == 0) NA else
  length(intersect(x, y)) / length(union(x, y))
J <- outer(ANN, ANN, Vectorize(function(i, j)
  jac(res[[i]]$deg_symbols, res[[j]]$deg_symbols)))
dimnames(J) <- list(ANN, ANN)
write.table(round(J, 3), file.path(out_dir, "comparison_jaccard.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

# ---- Layer 1b: gained / lost DEGs vs the Ensembl anchor ----
if (ANCHOR %in% ANN) {
  base <- res[[ANCHOR]]$deg_symbols
  cmp <- do.call(rbind, lapply(setdiff(ANN, ANCHOR), function(a) {
    s <- res[[a]]$deg_symbols
    data.frame(annotation = a,
               shared  = length(intersect(s, base)),
               gained  = length(setdiff(s, base)),   # DEG here, not in Ensembl
               lost    = length(setdiff(base, s)),   # DEG in Ensembl, not here
               jaccard = round(jac(s, base), 3))
  }))
  write.table(cmp, file.path(out_dir, paste0("comparison_vs_", ANCHOR, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# ---- Layer 1c: Spearman correlation of logFC on shared genes ----
C <- outer(ANN, ANN, Vectorize(function(i, j) {
  a <- res[[i]]$lfc; b <- res[[j]]$lfc
  g <- intersect(names(a), names(b)); g <- g[!is.na(g)]
  if (length(g) < 3) NA else round(cor(a[g], b[g], method = "spearman"), 3)
}))
dimnames(C) <- list(ANN, ANN)
write.table(C, file.path(out_dir, "comparison_logFC_cor.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

message("Done. Comparison tables written to ", out_dir)
# ============================================================================
# NEXT: feed the discordant genes (DEG in some annotations but not others) into
# the "fragility score" analysis (isoform count, exon overlap, mappability,
# biotype) to model WHICH genes flip when the annotation changes.
# ============================================================================
