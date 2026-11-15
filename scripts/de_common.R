#!/usr/bin/env Rscript
# ============================================================================
# de_common.R: sensitivity analysis: all three organisms analysed
# with ONE common DE method and threshold, and a scale for the annotation effect
# (how much the DEG list moves when only the DE method is changed).
#
# Common methods, applied identically to every organism and annotation:
#   edgeR : filterByExpr, TMM, quasi-likelihood F-test, FDR < 0.05, no fold-change cut-off
#   DESeq2: same genes as edgeR (filterByExpr), Wald test, padj < 0.05, no fold-change cut-off
# Design: ~ condition only (mouse KO vs WT; pig UH vs C; dog UD vs C, 22 libraries).
# Annotation effect (Ensembl vs RefSeq): genes matched by genomic coordinates
# (tables/coord_match/pairs.tsv.gz, criterion "shorter"); share of DEGs without a counterpart
# and Jaccard on loci, as in coord_match.R.
# Method effect (scale): Jaccard between the edgeR and DESeq2 DEG lists within one
# annotation (same gene identifiers), and between the anchor-study method and each
# common method.
# Writes: tables/de_common/{deg_lists.tsv.gz, annotation_effect.tsv, method_effect.tsv}
# ============================================================================
suppressMessages({library(edgeR); library(DESeq2); library(data.table)})
source("scripts/coord_match.R")          # components(); matching.R (read_de) via it
dir.create("tables/de_common", recursive = TRUE, showWarnings = FALSE)
P <- fread("tables/coord_match/pairs.tsv.gz")[shorter == TRUE]

samples <- function(o) {
  t <- read.delim(file.path(paste0(o, "_local"), "results", o, "targets.tsv"), colClasses = "character")
  if (o == "dog") { t <- t[t$UDvsC %in% c("T", "C"), ]; data.frame(sample = t$sample, grp = ifelse(t$UDvsC == "T", "test", "ref")) }
  else if (o == "pig") data.frame(sample = t$sample, grp = ifelse(t$condition == "UH", "test", "ref"))
  else data.frame(sample = t$sample, grp = ifelse(t$condition == "KO", "test", "ref"))
}
counts <- function(o, ann) {
  x <- read.delim(file.path(paste0(o, "_local"), "results", o, "counts", ann, "all_samples.counts.tsv"), check.names = FALSE)
  m <- as.matrix(x[, -1]); rownames(m) <- sub("\\.[0-9]+$", "", x$Geneid); m
}
de_common <- function(o, ann) {
  s <- samples(o); m <- counts(o, ann)[, s$sample]
  grp <- factor(s$grp, levels = c("ref", "test")); des <- model.matrix(~ grp)
  y <- DGEList(m); keep <- filterByExpr(y, des); y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
  fit <- glmQLFit(estimateDisp(y, des), des)
  qe <- topTags(glmQLFTest(fit, coef = 2), n = Inf)$table
  dd <- DESeq(DESeqDataSetFromMatrix(m[keep, ], data.frame(grp), ~ grp), quiet = TRUE)
  rd <- results(dd, name = "grp_test_vs_ref")
  list(edgeR = rownames(qe)[qe$FDR < 0.05], DESeq2 = rownames(rd)[!is.na(rd$padj) & rd$padj < 0.05])
}
jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
ann_effect <- function(o, degE, degR) {
  p <- P[organism == o]
  hasE <- degE %in% p$e; hasR <- degR %in% p$r
  comp <- components(p$e, p$r)
  lE <- unique(comp[paste0("E:", degE[hasE])]); lR <- unique(comp[paste0("R:", degR[hasR])])
  data.table(DEG_Ens = length(degE), DEG_Ref = length(degR), gap_Ens = 1 - mean(hasE), gap_Ref = 1 - mean(hasR),
             gap = 1 - (mean(hasE) + mean(hasR)) / 2, J = jac(lE, lR))
}

lists <- NULL; aeff <- NULL; meff <- NULL
for (o in c("mouse", "pig", "dog")) {
  cat("...", o, "\n")
  r <- list(ensembl = de_common(o, "ensembl"), refseq = de_common(o, "refseq"))
  anchor <- lapply(c(ensembl = "ensembl", refseq = "refseq"), function(a) {
    d <- read_de(o, a); unique(sub("\\.[0-9]+$", "", d$gene_id[d$is_deg])) })
  for (a in names(r)) for (mth in names(r[[a]]))
    lists <- rbind(lists, data.table(organism = o, annotation = a, method = mth, gene_id = r[[a]][[mth]]))
  for (mth in c("anchor", "edgeR", "DESeq2")) {
    degE <- if (mth == "anchor") anchor$ensembl else r$ensembl[[mth]]
    degR <- if (mth == "anchor") anchor$refseq  else r$refseq[[mth]]
    aeff <- rbind(aeff, cbind(data.table(organism = o, method = mth), ann_effect(o, degE, degR)))
  }
  for (a in names(r)) meff <- rbind(meff, data.table(organism = o, annotation = a,
    J_edgeR_DESeq2 = jac(r[[a]]$edgeR, r[[a]]$DESeq2),
    J_anchor_edgeR = jac(anchor[[a]], r[[a]]$edgeR), J_anchor_DESeq2 = jac(anchor[[a]], r[[a]]$DESeq2)))
}
num <- function(d) { d <- copy(d); n <- names(d)[vapply(d, is.double, TRUE)]; d[, (n) := lapply(.SD, round, 3), .SDcols = n]; d }
print(num(aeff)); print(num(meff))
fwrite(lists, "tables/de_common/deg_lists.tsv.gz", sep = "\t")
fwrite(aeff, "tables/de_common/annotation_effect.tsv", sep = "\t")
fwrite(meff, "tables/de_common/method_effect.tsv", sep = "\t")
