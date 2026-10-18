#!/usr/bin/env Rscript
# ============================================================================
# fragility_features.R: characterise ANNOTATION-FRAGILE vs STABLE DE genes.
# Template: works for any organism whose DESeq2 results + extracted GTF features
# are in place. Usage:  Rscript fragility_features.R <organism>
#
# "fragile" = a DEG in the anchor annotation (Ensembl) that is NOT a DEG under
# RefSeq; "stable" = DEG in both. We then compare the two groups on:
#   - biotype (protein_coding / lncRNA / pseudogene / predicted ...)
#   - number of isoforms (transcripts per gene)
#   - exonic length (union of exon bases)
#   - exon overlap with neighbouring genes  [proxy for mapping ambiguity]
#   - distance to the significance threshold (padj, |log2FC|)
#
# Needs (from the anchor Ensembl GTF, extracted on the cluster):
#   reference/<org>/annotations/ensembl_genes.tsv   gene_id biotype chr start end strand
#   reference/<org>/annotations/ensembl_ntx.tsv     gene_id n_transcripts
#   reference/<org>/annotations/ensembl_exons.tsv   gene_id chr start end
# and DESeq2 results:
#   results/<org>/de_deseq2/{ensembl,refseq}_deseq2.tsv
# ============================================================================
suppressMessages({ library(GenomicRanges); library(ggplot2) })

.args <- commandArgs(trailingOnly = TRUE)
org <- .args[1]
LFC <- if (length(.args) >= 2) as.numeric(.args[2]) else 1   # mouse/pig: 0 ; dog: 1.5
# 3rd arg: cross-reference key (mouse, pig, dog, dog10k) -> genes are matched between
# Ensembl and RefSeq by NCBI Gene ID (../xref/ens_<key>_entrez.tsv.gz, ../xref/refseq_<key>_geneid.tsv)
XKEY <- if (length(.args) >= 3) .args[3] else org
FDR <- 0.05
ann_dir <- file.path("reference", org, "annotations")
out_dir <- file.path("results", org, "fragility"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Read the DE output, whichever method was used: DESeq2 (de_deseq2/<ann>_deseq2.tsv,
# cols log2FoldChange/padj) or limma (de/<ann>_toptable.tsv, cols logFC/adj.P.Val).
read_de <- function(ann) {
  f2 <- file.path("results", org, "de_deseq2", paste0(ann, "_deseq2.tsv"))
  f1 <- file.path("results", org, "de",        paste0(ann, "_toptable.tsv"))
  if (file.exists(f2)) { d <- read.delim(f2); d$lfc <- d$log2FoldChange; d$padjv <- d$padj
                         d$expr <- log2(d$baseMean + 1) }         # DESeq2 expression
  else                 { d <- read.delim(f1); d$lfc <- d$logFC;          d$padjv <- d$adj.P.Val
                         d$expr <- d$AveExpr }                    # limma expression (log2-CPM)
  d
}
is_deg <- function(d) !is.na(d$padjv) & d$padjv < FDR & abs(d$lfc) > LFC

# ---- Ensembl DEGs labelled stable/fragile vs RefSeq, genes matched by NCBI Gene ID ----
ex_x <- read.delim(gzfile(file.path("..", "xref", paste0("ens_", XKEY, "_entrez.tsv.gz"))), stringsAsFactors = FALSE)
ex_x <- unique(ex_x[ex_x$db_name == "EntrezGene", c("gene_stable_id", "xref")])
ens2id <- split(as.character(ex_x$xref), ex_x$gene_stable_id)
rx <- read.delim(file.path("..", "xref", paste0("refseq_", XKEY, "_geneid.tsv")), header = FALSE,
                 col.names = c("gene_id", "GeneID"), colClasses = "character")
rx <- unique(rx[rx$GeneID != "", ]); ref2id <- setNames(rx$GeneID, rx$gene_id)
ref_universe <- unique(ref2id)                                   # every GeneID annotated by RefSeq
ids_of <- function(g) { x <- ens2id[match(g, names(ens2id))]; x[vapply(x, is.null, TRUE)] <- list(character(0)); unname(x) }
e <- read_de("ensembl"); r <- read_de("refseq")
rD_ids <- unique(unname(ref2id[r$gene_id[is_deg(r)]])); rD_ids <- rD_ids[!is.na(rD_ids)]
eDEG <- e[is_deg(e), ]
eids <- ids_of(eDEG$gene_id)
eDEG$present_in_refseq <- vapply(eids, function(v) any(v %in% ref_universe), TRUE)
eDEG$group <- ifelse(vapply(eids, function(v) any(v %in% rD_ids), TRUE), "stabilny", "kruchy")   # stable / fragile

# ---- gene-level features from the extracted GTF tables ----
g   <- read.delim(file.path(ann_dir, "ensembl_genes.tsv"), header = FALSE,
                  col.names = c("gene_id","biotype","chr","start","end","strand"))
ntx <- read.delim(file.path(ann_dir, "ensembl_ntx.tsv"), header = FALSE,
                  col.names = c("gene_id","n_tx"))
ex  <- read.delim(file.path(ann_dir, "ensembl_exons.tsv"), header = FALSE,
                  col.names = c("gene_id","chr","start","end"))

# reduced (union) exons per gene -> exonic length; and exon overlap between genes
exg <- GRanges(ex$chr, IRanges(ex$start, ex$end)); exg$gene <- ex$gene_id
red <- reduce(split(exg, exg$gene))                       # GRangesList per gene
exonic_len <- sum(width(red))                             # named vector
gr <- unlist(red); gr$gene <- rep(names(red), elementNROWS(red))
hits <- findOverlaps(gr, gr, ignore.strand = TRUE)
hits <- hits[gr$gene[queryHits(hits)] != gr$gene[subjectHits(hits)]]   # different genes only
q <- queryHits(hits); s <- subjectHits(hits)
iw <- pmax(pmin(end(gr)[q], end(gr)[s]) - pmax(start(gr)[q], start(gr)[s]) + 1, 0)
ov_bp <- tapply(iw, gr$gene[q], sum)

feat <- data.frame(gene_id = names(exonic_len),
                   exonic_len = as.integer(exonic_len),
                   overlap_bp = as.integer(ifelse(is.na(ov_bp[names(exonic_len)]), 0, ov_bp[names(exonic_len)])))
feat <- merge(feat, g[, c("gene_id","biotype")], by = "gene_id", all.x = TRUE)
feat <- merge(feat, ntx, by = "gene_id", all.x = TRUE)
feat$n_tx[is.na(feat$n_tx)] <- 1

# ---- RefSeq exonic length of the matched gene (for cross-annotation diff) ----
# Chisanga's insight: the fragility signal is in the DIFFERENCE of effective length
# BETWEEN annotations, not the absolute length. Needs (both files optional):
#   reference/<org>/annotations/refseq_exons.tsv   gene_id chr start end
#   reference/<org>/annotations/refseq.genemap.tsv gene_id symbol
rsym_len <- NULL
rex_f  <- file.path(ann_dir, "refseq_exons.tsv")
rmap_f <- file.path(ann_dir, "refseq.genemap.tsv")
if (file.exists(rex_f) && file.exists(rmap_f)) {
  rex  <- read.delim(rex_f, header = FALSE, col.names = c("gene_id","chr","start","end"))
  rg   <- GRanges(rex$chr, IRanges(rex$start, rex$end)); rg$gene <- rex$gene_id
  rlen <- sum(width(reduce(split(rg, rg$gene))))                    # RefSeq exonic len / gene_id
  rid <- ref2id[names(rlen)]                                        # RefSeq gene -> GeneID
  rsym_len <- tapply(as.integer(rlen)[!is.na(rid)], rid[!is.na(rid)], max)  # per GeneID
}

# ---- join features to the labelled Ensembl DEGs ----
d <- merge(eDEG[, c("gene_id","symbol","lfc","padjv","expr","group","present_in_refseq")], feat, by = "gene_id", all.x = TRUE)
d$abs_lfc      <- abs(d$lfc)
d$neglog10padj <- -log10(d$padjv)
# % unique exonic bases (Hamaguchi): fraction of the gene's exon bases that do NOT
# overlap another gene's exons -> 1 means fully private, low means shared/ambiguous
d$pct_unique_exon <- round(pmax(0, 1 - d$overlap_bp / pmax(d$exonic_len, 1)), 3)
# cross-annotation length difference + presence in RefSeq (mechanism descriptor)
if (!is.null(rsym_len)) {
  d$refseq_len <- vapply(ids_of(d$gene_id), function(v) { x <- rsym_len[v]; x <- x[!is.na(x)]
                         if (length(x)) as.numeric(max(x)) else NA_real_ }, 0)
  d$log2_len_ratio <- ifelse(!is.na(d$refseq_len), log2(d$exonic_len / d$refseq_len), NA)
} else { d$refseq_len <- NA; d$log2_len_ratio <- NA }
d$biotype_grp  <- ifelse(d$biotype %in% c("protein_coding"), "protein_coding",
                  ifelse(grepl("pseudogene", d$biotype), "pseudogene",
                  ifelse(d$biotype %in% c("lncRNA","lincRNA"), "lncRNA", "inne")))
d <- d[order(d$group, d$padjv), ]
write.table(d, file.path(out_dir, "features_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- summary table (fragile vs stable) ----
summ <- do.call(rbind, lapply(split(d, d$group), function(x) data.frame(
  grupa = x$group[1], n = nrow(x),
  med_izoform = median(x$n_tx), med_dl_egzon = median(x$exonic_len),
  med_overlap_bp = median(x$overlap_bp),
  frakcja_z_overlap = round(mean(x$overlap_bp > 0), 2),
  med_pct_unik_egzon = round(median(x$pct_unique_exon), 3),
  pct_bez_refseq = round(mean(!x$present_in_refseq, na.rm = TRUE) * 100),
  med_abs_dl_ratio = round(median(abs(x$log2_len_ratio), na.rm = TRUE), 2),
  pct_pseudogen = round(mean(x$biotype_grp == "pseudogene") * 100),
  pct_protein   = round(mean(x$biotype_grp == "protein_coding") * 100),
  med_ekspresja = round(median(x$expr), 2),
  med_neglog10padj = round(median(x$neglog10padj), 2))))
cat("\n=== Podsumowanie: kruchy vs stabilny (Ensembl DEG) ===\n"); print(summ, row.names = FALSE)
write.table(summ, file.path(out_dir, "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n=== Skład biotypów ===\n"); print(table(d$group, d$biotype_grp))

# ---- plots ----
pal <- c(kruchy = "#D85A30", stabilny = "#1D9E75")
# biotype composition
p1 <- ggplot(d, aes(group, fill = biotype_grp)) + geom_bar(position = "fill") +
  labs(x = NULL, y = "udział", fill = "biotyp", title = "Biotyp: kruchy vs stabilny") +
  theme_minimal(base_size = 13)
ggsave(file.path(out_dir, "fragility_biotype.png"), p1, width = 5.5, height = 4, dpi = 130)

# numeric features (boxplots, faceted)
library(stats)
long <- rbind(
  data.frame(group = d$group, feature = "liczba izoform",        value = d$n_tx),
  data.frame(group = d$group, feature = "dlugosc egzonowa (bp)", value = d$exonic_len),
  data.frame(group = d$group, feature = "nakladanie eksonow (bp)", value = d$overlap_bp + 1),
  data.frame(group = d$group, feature = "% unikalnych eksonow",   value = pmax(d$pct_unique_exon, 0.001)),
  data.frame(group = d$group, feature = "|log2 dl. Ens/RefSeq|",  value = pmax(abs(d$log2_len_ratio), 0.001)),
  data.frame(group = d$group, feature = "ekspresja (log2)",       value = pmax(d$expr, 0.001)),
  data.frame(group = d$group, feature = "-log10(padj)",          value = d$neglog10padj))
p2 <- ggplot(long, aes(group, value, fill = group)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.8) + geom_jitter(width = 0.15, size = 0.6, alpha = 0.5) +
  facet_wrap(~feature, scales = "free_y") + scale_fill_manual(values = pal) +
  scale_y_log10() + labs(x = NULL, y = "wartosc (log10)", title = "Cechy: kruchy vs stabilny") +
  theme_minimal(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(out_dir, "fragility_features.png"), p2, width = 8, height = 5, dpi = 130)

cat("\nZapisano do", out_dir, ": features_table.tsv, summary.tsv, fragility_biotype.png, fragility_features.png\n")
