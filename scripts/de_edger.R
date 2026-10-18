#!/usr/bin/env Rscript
# ============================================================================
# de_edger.R: dog differential expression with edgeR, as described in the
# anchor study (Stachowiak et al. 2024, PNAS 121:e2312724121). The SAME analysis
# is run for every annotation (ensembl/refseq/ucsc), so the only thing that
# varies is the annotation.
#
# Published method (Methods and supplementary Datasets S1 and S4 of the paper):
#   - libraries: the 22 libraries analysed in the study (Dataset S4, panel 1),
#     all deposited in SRA PRJNA901164
#   - filter: CPM > 0.5 in more than two samples
#   - TMM normalisation, edgeR likelihood-ratio test (Dataset S1 reports LR)
#   - DEG criteria: FDR < 0.05 and |log2FC| > 1.5
#   - the published results table (Dataset S1) contains protein-coding genes
#     only; this restriction is applied for the calibration only (optional 4th
#     argument), the main analysis keeps all biotypes.
#
# Contrasts (column in targets encodes membership; value T=test, C=reference):
#   UDvsC  undescended vs normal control     (unpaired; design ~condition)
#   DvsUD  descended vs undescended          (paired within dog; design ~dog_id+condition)
#   DvsC   descended vs normal control       (unpaired)
#
# Usage:  Rscript de_edger.R <organism> <contrast> [LFC] [protein_coding_ids]
#   e.g.  cd dog_local && Rscript ../scripts/de_edger.R dog UDvsC 1.5
# Reads : results/<org>/counts/<ann>/all_samples.counts.tsv
#         results/<org>/targets.tsv
#         reference/<org>/annotations/<ann>.genemap.tsv   (gene_id -> symbol)
# Writes: results/<org>/de/<ann>_toptable.tsv  (logFC, AveExpr, LR, P.Value, adj.P.Val, symbol)
# ============================================================================
suppressMessages(library(edgeR))

.args    <- commandArgs(trailingOnly = TRUE)
org      <- .args[1]
contrast <- if (length(.args) >= 2) .args[2] else "UDvsC"
LFC      <- if (length(.args) >= 3) as.numeric(.args[3]) else 1.5
PC_IDS   <- if (length(.args) >= 4) readLines(.args[4]) else NULL
FDR      <- 0.05
anns     <- c("ensembl", "refseq", "ucsc")

targets <- read.delim(file.path("results", org, "targets.tsv"),
                      colClasses = "character", check.names = FALSE)
stopifnot(contrast %in% colnames(targets))
out_dir <- file.path("results", org, "de"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# genemap: gene_id -> symbol (fallback to gene_id)
load_symbols <- function(ann) {
  f <- file.path("reference", org, "annotations", paste0(ann, ".genemap.tsv"))
  if (!file.exists(f)) return(NULL)
  m <- read.delim(f, header = FALSE, col.names = c("gene_id", "symbol"))
  setNames(ifelse(m$symbol == "" | is.na(m$symbol), m$gene_id, m$symbol), m$gene_id)
}

run_one <- function(ann) {
  cf <- file.path("results", org, "counts", ann, "all_samples.counts.tsv")
  if (!file.exists(cf)) { cat("  [skip]", ann, "- missing", cf, "\n"); return(invisible()) }
  cnt <- read.delim(cf, check.names = FALSE)
  rownames(cnt) <- cnt$Geneid; cnt$Geneid <- NULL
  cnt <- as.matrix(cnt)
  if (!is.null(PC_IDS)) cnt <- cnt[rownames(cnt) %in% PC_IDS, , drop = FALSE]

  # samples in this contrast: non-empty membership cell; value C=reference, T=test
  keep_s <- targets$sample[targets[[contrast]] %in% c("C", "T")]
  keep_s <- keep_s[keep_s %in% colnames(cnt)]
  ti <- targets[match(keep_s, targets$sample), ]
  mat <- cnt[, keep_s, drop = FALSE]
  cond <- factor(ti[[contrast]], levels = c("C", "T"))     # C is the reference level

  y <- DGEList(counts = mat, group = cond)
  keep <- rowSums(cpm(y) > 0.5) > 2                          # CPM > 0.5 in more than two samples
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")

  paired <- (contrast == "DvsUD")
  design <- if (paired) model.matrix(~ factor(ti$dog_id) + cond) else model.matrix(~ cond)
  y   <- estimateDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = ncol(design))                    # condT
  tt  <- topTags(lrt, n = Inf, sort.by = "none")$table

  sym <- load_symbols(ann)
  gid <- rownames(tt)
  out <- data.frame(gene_id = gid,
                    symbol  = if (!is.null(sym)) ifelse(gid %in% names(sym), sym[gid], gid) else gid,
                    logFC   = tt$logFC,
                    AveExpr = tt$logCPM,                      # expression (log2-CPM) for fragility
                    LR      = tt$LR,
                    P.Value = tt$PValue,
                    adj.P.Val = tt$FDR,
                    stringsAsFactors = FALSE)
  write.table(out, file.path(out_dir, paste0(ann, "_toptable.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  ndeg <- sum(out$adj.P.Val < FDR & abs(out$logFC) > LFC, na.rm = TRUE)
  cat(sprintf("  %-8s: %d genes after filtering | DEG (FDR<%.2f & |log2FC|>%.1f) = %d\n",
              ann, nrow(out), FDR, LFC, ndeg))
}

cat(sprintf("=== dog edgeR | contrast %s | |log2FC|>%.1f, FDR<%.2f | %s ===\n",
            contrast, LFC, FDR, if (contrast == "DvsUD") "paired (block = dog)" else "unpaired"))
for (a in anns) run_one(a)
cat("Done -> results/", org, "/de/<ann>_toptable.tsv\n", sep = "")
