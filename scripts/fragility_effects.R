#!/usr/bin/env Rscript
# ============================================================================
# fragility_effects.R: effect sizes for fragile vs stable genes (Table 3).
# For each organism and feature: AUC = P(feature of a random fragile gene >
# feature of a random stable gene), ties counted 1/2 (0.5 = no separation;
# < 0.5 = fragile genes lower). 95% CI by stratified bootstrap (2000 resamples
# within each class). Two-sided p: Wilcoxon rank-sum (continuous features) or
# Fisher's exact test (binary features). Genes are not fully independent, so
# the AUC and its CI are the primary summary; p-values are descriptive.
# Reads  : <org>_local/results/<org>/fragility/features_table.tsv
# Writes : results_fragility_effects.tsv
# Usage  : Rscript fragility_effects.R   (from the project root)
# ============================================================================
set.seed(1)
orgs <- c("mouse", "pig", "dog")
feats <- list(
  protein_coding = function(d) as.numeric(d$biotype_grp == "protein_coding"),
  expression     = function(d) d$expr,
  exon_overlap   = function(d) as.numeric(d$overlap_bp > 0),
  abs_len_ratio  = function(d) abs(d$log2_len_ratio))
binary <- c("protein_coding", "exon_overlap")

auc <- function(x, y) {                      # x = fragile, y = stable
  r <- rank(c(x, y)); nx <- length(x); ny <- length(y)
  (sum(r[seq_len(nx)]) - nx * (nx + 1) / 2) / (nx * ny)
}

out <- NULL
for (o in orgs) {
  d <- read.delim(file.path(paste0(o, "_local"), "results", o, "fragility", "features_table.tsv"))
  for (f in names(feats)) {
    v <- feats[[f]](d); x <- v[d$group == "kruchy"]; y <- v[d$group == "stabilny"]
    x <- x[is.finite(x)]; y <- y[is.finite(y)]
    a <- auc(x, y)
    bs <- replicate(2000, auc(sample(x, replace = TRUE), sample(y, replace = TRUE)))
    p <- if (f %in% binary) {
      fisher.test(table(factor(c(x, y), levels = 0:1), rep(c("f", "s"), c(length(x), length(y)))))$p.value
    } else suppressWarnings(wilcox.test(x, y)$p.value)
    out <- rbind(out, data.frame(organism = o, feature = f, n_fragile = length(x), n_stable = length(y),
                                 AUC = round(a, 2), CI_low = round(quantile(bs, 0.025), 2),
                                 CI_high = round(quantile(bs, 0.975), 2), p = signif(p, 2)))
  }
}
rownames(out) <- NULL
print(out)
write.table(out, "tables/fragility_effects.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
