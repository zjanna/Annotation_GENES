#!/usr/bin/env Rscript
# ============================================================================
# fragility_coord.R: annotation-sensitive ("fragile") DEGs, with each database in
# turn as the reference and genes matched by genomic coordinates (coord_match.R,
# criterion "shorter"); the two routes to fragility are kept apart.
#
# For a DEG of database X (Ensembl or RefSeq), compared with database Y:
#   stable           : it has a coordinate match in Y that is a DEG of Y
#   fragile, absent  : it has no coordinate match in Y
#   fragile, present : it has a coordinate match in Y, but no matched gene is a DEG of Y
# Features (all from X's own annotation and X's own DE fit):
#   protein-coding (gene biotype of X), expression (X's DE fit; scale differs by DE
#   method), exon overlap with another gene of X (>= 1 bp, strand ignored),
#   |log2 exonic length X / matched gene of Y| (best match = most shared bases;
#   matched genes only).
# Effect sizes: AUC (P(feature of a random fragile gene > that of a random stable
# gene), ties 1/2) with a 95% CI from 2000 stratified bootstrap resamples; no p-values
# (with thousands of genes they carry no information).
# Writes: tables/fragility_coord/{classes.tsv.gz, summary.tsv, effects.tsv}
# Usage : Rscript scripts/fragility_coord.R   (after scripts/coord_match.R)
# ============================================================================
suppressMessages({library(data.table); library(GenomicRanges)})
source("scripts/coord_match.R")          # read_exons(), refseq_chr_map(), SRC, G; matching.R via it
set.seed(1)
dir.create("tables/fragility_coord", recursive = TRUE, showWarnings = FALSE)
P <- fread("tables/coord_match/pairs.tsv.gz")[shorter == TRUE]

read_biotype <- function(files) {
  rbindlist(lapply(file.path(G, files), function(f) {
    cmd <- if (grepl("\\.gz$", f)) paste("gzip -dc", shQuote(f)) else paste("cat", shQuote(f))
    x <- fread(cmd = paste(cmd, "| awk -F'\\t' '$3==\"gene\"'"), sep = "\t", header = FALSE, select = 9, quote = "")
    data.table(gene_id = sub('^.*gene_id "([^"]*)".*$', "\\1", x$V9),
               biotype = sub('^.*gene_biotype "([^"]*)".*$', "\\1", x$V9))
  }))[, gene_id := sub("\\.[0-9]+$", "", gene_id)]
}
# genes with >= 1 bp of exon overlap with another gene of the same annotation
overlapping_genes <- function(x) {
  gr <- GRanges(x$chr, IRanges(x$start, x$end), gene_id = x$gene_id)
  red <- unlist(reduce(split(gr, gr$gene_id))); red$gene_id <- names(red)
  h <- findOverlaps(red, red, ignore.strand = TRUE)
  h <- h[red$gene_id[queryHits(h)] != red$gene_id[subjectHits(h)]]
  unique(red$gene_id[queryHits(h)])
}
auc <- function(x, y) { r <- rank(c(x, y)); nx <- length(x); ny <- length(y)
  (sum(r[seq_len(nx)]) - nx * (nx + 1) / 2) / (nx * ny) }
auc_ci <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 3 || length(y) < 3) return(c(n = length(x), AUC = NA, lo = NA, hi = NA))
  bs <- replicate(2000, auc(sample(x, replace = TRUE), sample(y, replace = TRUE)))
  c(n = length(x), AUC = auc(x, y), lo = unname(quantile(bs, 0.025)), hi = unname(quantile(bs, 0.975)))
}

classes <- NULL
for (k in c("mouse", "pig", "dog")) {
  cat("...", k, "\n")
  s <- SRC[[k]]
  dE <- as.data.table(read_de(k, "ensembl"))[, gene_id := sub("\\.[0-9]+$", "", gene_id)]
  dR <- as.data.table(read_de(k, "refseq"))
  exprcol <- if ("baseMean" %in% names(dE)) "baseMean" else "AveExpr"
  ex <- function(d) if (exprcol == "baseMean") log2(d$baseMean + 1) else d$AveExpr
  dE[, expr := ex(dE)]; dR[, expr := ex(dR)]
  degE <- dE[is_deg == TRUE]; degR <- dR[is_deg == TRUE]
  p <- P[organism == k]
  bE <- read_biotype(s$ens); bR <- read_biotype(s$ref)
  xE <- read_exons(s$ens)[, gene_id := sub("\\.[0-9]+$", "", gene_id)]; xR <- read_exons(s$ref)
  ovE <- overlapping_genes(xE); ovR <- overlapping_genes(xR)
  one_side <- function(deg, own, other_deg, own_col, oth_col, bt, ov, len_own, len_oth, db) {
    pm <- p[, .(g = get(own_col), h = get(oth_col), S, lo = get(len_own), lt = get(len_oth))]
    m <- pm[g %in% deg$gene_id]
    st <- m[h %in% other_deg, unique(g)]
    best <- m[order(-S)][!duplicated(g)]
    out <- data.table(organism = k, db = db, gene_id = deg$gene_id, expr = deg$expr)
    out[, class := fifelse(gene_id %in% st, "stable", fifelse(gene_id %in% m$g, "fragile, present", "fragile, absent"))]
    out[, biotype := bt$biotype[match(gene_id, bt$gene_id)]]
    out[, protein_coding := as.numeric(biotype %in% "protein_coding")]
    out[, exon_overlap := as.numeric(gene_id %in% ov)]
    out[, abs_log2_len := abs(log2(best$lo[match(gene_id, best$g)] / best$lt[match(gene_id, best$g)]))]
    out
  }
  classes <- rbind(classes,
    one_side(degE, "ensembl", degR$gene_id, "e", "r", bE, ovE, "LE", "LR", "Ensembl"),
    one_side(degR, "refseq",  degE$gene_id, "r", "e", bR, ovR, "LR", "LE", "RefSeq"))
}
fwrite(classes, "tables/fragility_coord/classes.tsv.gz", sep = "\t")

lev <- c("stable", "fragile, present", "fragile, absent")
classes[, class := factor(class, levels = lev)]
summ <- classes[, .(n = .N, pct_protein_coding = round(100 * mean(protein_coding), 0),
                    median_expr = round(median(expr), 2),
                    pct_exon_overlap = round(100 * mean(exon_overlap), 0),
                    median_abs_log2_len = round(median(abs_log2_len, na.rm = TRUE), 2)),
                by = .(organism, db, class)][order(match(organism, c("mouse", "pig", "dog")), db, class)]
print(summ)
fwrite(summ, "tables/fragility_coord/summary.tsv", sep = "\t")

feats <- c("protein_coding", "expr", "exon_overlap", "abs_log2_len")
eff <- NULL
for (o in c("mouse", "pig", "dog")) for (d in c("Ensembl", "RefSeq")) for (cl in lev[-1]) for (f in feats) {
  z <- classes[organism == o & db == d]
  r <- auc_ci(z[class == cl][[f]], z[class == "stable"][[f]])
  eff <- rbind(eff, data.table(organism = o, db = d, contrast = paste(cl, "vs stable"), feature = f,
                               n_fragile = r[["n"]], n_stable = sum(z$class == "stable" & is.finite(z[[f]])),
                               AUC = round(r[["AUC"]], 2), CI_low = round(r[["lo"]], 2), CI_high = round(r[["hi"]], 2)))
}
print(eff[feature %in% c("protein_coding", "expr")])
fwrite(eff, "tables/fragility_coord/effects.tsv", sep = "\t")
