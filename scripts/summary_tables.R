#!/usr/bin/env Rscript
# ============================================================================
# summary_tables.R: numbers quoted in the Results (run from the project root).
#   tables/agreement_ensembl_refseq.tsv : per organism, Ensembl vs RefSeq
#       DEGs per database, comparable DEGs, agreement among comparable genes
#       (Jaccard, NCBI Gene ID), share of DEGs without a counterpart (mean and
#       per side), DEG decomposition used in Fig. 4
#   tables/agreement_all_pairs.tsv      : Jaccard for every database pair (Fig. 1)
# ============================================================================
source("scripts/matching.R")
dir.create("tables", showWarnings = FALSE)

core <- do.call(rbind, lapply(c("mouse", "pig", "dog", "dog10k"), function(k) {
  r <- compare_dbs(k, "ensembl", "refseq")
  data.frame(key = k, DEG_ensembl = r$n_deg_x, DEG_refseq = r$n_deg_y,
             comparable_ensembl = sum(r$dx$is_deg & r$cx), comparable_refseq = sum(r$dy$is_deg & r$cy),
             jaccard_comparable = round(r$J, 3),
             no_counterpart_mean = round(r$coverage, 3),
             no_counterpart_ensembl = round(r$gap_x, 3), no_counterpart_refseq = round(r$gap_y, 3),
             shared_geneids = r$shared, ensembl_only_geneids = r$only_x, refseq_only_geneids = r$only_y)
}))
print(core, row.names = FALSE)
write.table(core, "tables/agreement_ensembl_refseq.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

dbs <- list(mouse = c("ensembl", "gencode", "refseq", "ucsc"), pig = c("ensembl", "refseq", "ucsc"),
            dog = c("ensembl", "refseq", "ucsc"))
pairs <- NULL
for (k in names(dbs)) { d <- dbs[[k]]
  for (i in seq_along(d)) for (j in seq_along(d)) if (i < j)
    pairs <- rbind(pairs, data.frame(key = k, x = d[i], y = d[j],
                                     jaccard = round(compare_dbs(k, d[i], d[j])$J, 3)))
}
print(pairs, row.names = FALSE)
write.table(pairs, "tables/agreement_all_pairs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Table 4: the eleven pig genes tested by qPCR in the anchor study (8 confirmed)
qpcr <- data.frame(gene = c("SIM1", "PITX1", "HOXA7", "METTL21C", "PVALB", "ALX1", "EYA2", "TBX1",
                            "IRX5", "ALX4", "OSR2"),
                   qpcr_confirmed = rep(c("yes", "no"), c(8, 3)))
pe <- read_de("pig", "ensembl"); pr <- read_de("pig", "refseq")
deg_of <- function(d, g) vapply(g, function(s) any(d$is_deg[toupper(d$symbol) == s]), TRUE)
qpcr$ensembl_deg <- deg_of(pe, qpcr$gene); qpcr$refseq_deg <- deg_of(pr, qpcr$gene)
qpcr$status <- ifelse(qpcr$ensembl_deg & qpcr$refseq_deg, "stable",
               ifelse(qpcr$ensembl_deg, "fragile", ifelse(qpcr$refseq_deg, "RefSeq only", "not DE")))
print(qpcr, row.names = FALSE)
write.table(qpcr, "tables/qpcr_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
