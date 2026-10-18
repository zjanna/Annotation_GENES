#!/usr/bin/env Rscript
# DESeq2 DE per annotation, Efcab4b anchor (Cheng et al. 2025).
# The paper's Methods state padj<0.05 & |log2FC|>1, but its published 69 DEGs
# include genes with |log2FC|<1 (e.g. Tat -0.68, Xbp1 0.66): the list was cut at
# padj<0.05 only. Default LFC is therefore 0; pass 1 to get the Methods-text rule.
# NB targets.tsv: conditions follow the FASTQ file names (KO1L-KO3L = KO), which
# match the paper's fold-change signs; the ArrayExpress SDRF genotype labels are
# swapped. Sample IDs (WT_1..KO_3) are kept from the first pass for continuity.
# Usage: Rscript de_deseq2.R <organism> [LFC]
suppressMessages(library(DESeq2))

.args  <- commandArgs(trailingOnly = TRUE)
org    <- .args[1]
ANCHOR <- "ensembl"; FDR <- 0.05; LFC <- if (length(.args) >= 2) as.numeric(.args[2]) else 0
annot_dir <- file.path("reference", org, "annotations")
count_of  <- function(a) file.path("results", org, "counts", a, "all_samples.counts.tsv")
ANN <- c("ensembl", "refseq", "gencode", "ucsc")
ANN <- ANN[file.exists(vapply(ANN, count_of, ""))]

targets <- read.delim(file.path("results", org, "targets.tsv"), stringsAsFactors = FALSE)
out_dir <- file.path("results", org, "de_deseq2"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Reference (control) level auto-detected; contrast = treated vs control.
ctrl_names <- c("WT","C","control","Control","ctrl","CTRL","healthy","normal")
lev <- unique(targets$condition)
REF <- intersect(ctrl_names, lev); REF <- if (length(REF)) REF[1] else sort(lev)[1]
TREAT <- setdiff(lev, REF)[1]
COND_LEVELS <- c(REF, setdiff(lev, REF))
HAS_BLOCK <- "block" %in% colnames(targets)   # e.g. litter -> ~block+condition (like limma dupCor)
cat(sprintf("Design: %s | reference=%s, treated=%s%s\n",
    ifelse(HAS_BLOCK, "~block+condition", "~condition"), REF, TREAT,
    ifelse(HAS_BLOCK, "  (block in model)", "")))

symmap <- function(a) {
  mf <- file.path(annot_dir, paste0(a, ".genemap.tsv"))
  if (file.exists(mf)) { m <- read.delim(mf, header = FALSE, stringsAsFactors = FALSE)
    return(setNames(toupper(m[[2]]), m[[1]])) }
  setNames(character(0), character(0))
}

run <- function(a) {
  m <- read.delim(count_of(a), check.names = FALSE); rownames(m) <- m$Geneid; m$Geneid <- NULL
  m <- round(as.matrix(m[, targets$sample]))
  cd <- data.frame(condition = factor(targets$condition, levels = COND_LEVELS),
                   row.names = targets$sample)
  design <- ~condition
  if (HAS_BLOCK) { cd$block <- factor(targets$block); design <- ~block + condition }
  dds <- DESeqDataSetFromMatrix(m, cd, design)
  dds <- DESeq(dds, quiet = TRUE)
  res <- as.data.frame(results(dds, contrast = c("condition", TREAT, REF)))
  res$gene_id <- rownames(res)
  sm <- symmap(a); res$symbol <- sm[res$gene_id]
  res$symbol[is.na(res$symbol)] <- res$gene_id[is.na(res$symbol)]
  write.table(res, file.path(out_dir, paste0(a, "_deseq2.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  deg <- res[!is.na(res$padj) & res$padj < FDR & abs(res$log2FoldChange) > LFC, ]
  cat(sprintf("  %-8s: %5d tested, %3d DEG\n", a, sum(!is.na(res$padj)), nrow(deg)))
  list(deg = unique(deg$symbol), lfc = setNames(res$log2FoldChange, res$symbol))
}

cat(sprintf("DESeq2 per annotation (padj<0.05 & |log2FC|>%g):\n", LFC))
res <- setNames(lapply(ANN, run), ANN)

jac <- function(x, y) if (length(union(x, y)) == 0) NA else length(intersect(x, y)) / length(union(x, y))
J <- outer(ANN, ANN, Vectorize(function(i, j) round(jac(res[[i]]$deg, res[[j]]$deg), 3)))
dimnames(J) <- list(ANN, ANN)
cat("\nJaccard of DEG sets (by symbol):\n"); print(J)

if (ANCHOR %in% ANN) {
  base <- res[[ANCHOR]]$deg
  cat(sprintf("\nvs %s:\n", ANCHOR))
  for (a in setdiff(ANN, ANCHOR))
    cat(sprintf("  %-8s: shared %d, gained %d, lost %d\n", a,
        length(intersect(res[[a]]$deg, base)),
        length(setdiff(res[[a]]$deg, base)),
        length(setdiff(base, res[[a]]$deg))))
}
write.table(J, file.path(out_dir, "comparison_jaccard.tsv"), sep = "\t", quote = FALSE, col.names = NA)
