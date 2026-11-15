#!/usr/bin/env Rscript
# ============================================================================
# topn_agreement.R: Ensembl vs RefSeq agreement measures that do not depend on the
# size of the DEG lists: overlap of the top-N genes ranked by the test statistic
# (|Wald stat| DESeq2, |t| limma, LR edgeR; p-values underflow to 0 for the strongest
# dog genes, so they cannot rank them), and the correlation of log2 fold changes, on
# genes matched one-to-one by genomic coordinates (tables/coord_match/pairs.tsv.gz, criterion "shorter"; one-to-one
# pairs formed greedily from the largest shared exonic length) and tested in both
# databases. Anchor-study DE results (matching.R, read_de).
#   top-N agreement = |top-N Ensembl ∩ top-N RefSeq| / N   (ranks within the paired genes)
# Writes: tables/topn_agreement/topn_agreement.tsv
# ============================================================================
suppressMessages(library(data.table))
source("scripts/matching.R")
dir.create("tables/topn_agreement", recursive = TRUE, showWarnings = FALSE)
P <- fread("tables/coord_match/pairs.tsv.gz")[shorter == TRUE]
pcol <- c(mouse = "stat", pig = "t", dog = "LR", dog10k = "LR")       # ranking statistic

one_to_one <- function(p) {          # greedy: largest shared length first
  p <- p[order(-S)]; ue <- ur <- character(0); keep <- logical(nrow(p))
  for (i in seq_len(nrow(p))) if (!(p$e[i] %in% ue) && !(p$r[i] %in% ur)) {
    keep[i] <- TRUE; ue <- c(ue, p$e[i]); ur <- c(ur, p$r[i]) }
  p[keep]
}
out <- NULL
for (k in names(pcol)) {
  dE <- as.data.table(read_de(k, "ensembl"))[, gene_id := sub("\\.[0-9]+$", "", gene_id)]
  dR <- as.data.table(read_de(k, "refseq"))
  dE <- dE[!is.na(get(pcol[k]))]; dR <- dR[!is.na(get(pcol[k]))]
  p <- one_to_one(P[organism == k & e %in% dE$gene_id & r %in% dR$gene_id])
  ie <- match(p$e, dE$gene_id); ir <- match(p$r, dR$gene_id)
  m <- data.table(pe = dE[[pcol[k]]][ie], pr = dR[[pcol[k]]][ir], le = dE$lfc_value[ie], lr = dR$lfc_value[ir],
                  de = dE$is_deg[ie], dr = dR$is_deg[ir])
  m[, `:=`(rank_e = frank(-abs(pe), ties.method = "first"), rank_r = frank(-abs(pr), ties.method = "first"))]
  top <- function(N) mean(which(m$rank_e <= N) %in% which(m$rank_r <= N))
  both <- m$de & m$dr
  out <- rbind(out, data.table(organism = k, n_pairs = nrow(m),
    top100 = round(top(100), 3), top500 = round(top(500), 3), top1000 = round(top(1000), 3),
    lfc_pearson = round(cor(m$le, m$lr), 3), lfc_spearman = round(cor(m$le, m$lr, method = "spearman"), 3),
    sign_agreement_DEG_both = round(100 * mean(sign(m$le[both]) == sign(m$lr[both])), 1)))
}
print(out)
fwrite(out, "tables/topn_agreement/topn_agreement.tsv", sep = "\t")
