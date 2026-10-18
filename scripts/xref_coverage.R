#!/usr/bin/env Rscript
# ============================================================================
# xref_coverage.R: share of Ensembl genes cross-referenced to an NCBI Gene ID
# (Ensembl's own EntrezGene xrefs), overall, for protein-coding genes and for all
# other biotypes. Quoted in the Results (dog on the current reference): dog
# Ensembl 92 vs 112, and pig Ensembl 104 vs 112. Run from the project root.
# Writes tables/geneid_coverage.tsv
# ============================================================================
dir.create("tables", showWarnings = FALSE)
sets <- list(
  c("dog (Ensembl 92)",   "dog_local/reference/dog/annotations/ensembl_genes.tsv",       "xref/ens_dog_entrez.tsv.gz"),
  c("dog (Ensembl 112)",  "dog10k_local/reference/dog/annotations/ensembl_genes.tsv",    "xref/ens_dog10k_entrez.tsv.gz"),
  c("pig (Ensembl 104)",  "pig_local/reference/pig/annotations/ensembl_genes.tsv",       "xref/ens_pig_entrez.tsv.gz"),
  c("pig (Ensembl 112)",  "pig_local/reference/pig/annotations/ensembl112_genes.tsv",    "xref/ens_pig112_entrez.tsv.gz"))
res <- do.call(rbind, lapply(sets, function(s) {
  g <- read.delim(s[2], header = FALSE, stringsAsFactors = FALSE)[, 1:2]; names(g) <- c("id", "biotype")
  x <- read.delim(gzfile(s[3]), stringsAsFactors = FALSE)
  g$xref <- g$id %in% x$gene_stable_id[x$db_name == "EntrezGene"]
  pc <- g$biotype == "protein_coding"
  data.frame(annotation = s[1], genes = nrow(g), pct_geneid = round(100 * mean(g$xref)),
             pct_geneid_protein_coding = round(100 * mean(g$xref[pc])),
             pct_geneid_other = round(100 * mean(g$xref[!pc])))
}))
print(res, row.names = FALSE)
write.table(res, "tables/geneid_coverage.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
