#!/usr/bin/env Rscript
# ============================================================================
# coord_match.R: Ensembl vs RefSeq gene matching by genomic coordinates
# (strand-specific exon overlap), as an alternative to matching by NCBI Gene ID
# (matching.R). Answers: are DEGs "without a counterpart" really absent from the
# other database, or only missing a Gene ID cross-reference?
#
# Inputs (reference/downloads/, from workflow/reference_tables.sh): the Ensembl and RefSeq
# GTFs of the releases used, and the NCBI assembly reports. RefSeq chromosome names are
# converted to Ensembl names (assembled molecules: chromosome name; other sequences:
# GenBank accession).
#
# For every gene: union of its exons (exonic length L). For every Ensembl-RefSeq gene
# pair on the same chromosome and strand: shared exonic bp S. A pair is a match if
#   any     : S >= 1 bp
#   shorter : S >= 0.5 * min(L_E, L_R)                     (primary criterion)
#   recip   : S >= 0.5 * L_E and S >= 0.5 * L_R
# Coverage gap of X = share of DEGs of X with no matched gene anywhere in Y's GTF.
# Agreement: matched genes are joined into loci (connected components of the match
# graph); a locus is DE under X if any of its X genes is a DEG of X;
#   J = |loci DE under both| / |loci DE under either|.
# Usage: Rscript scripts/coord_match.R    (from the repository root; writes tables/coord_match/)
# ============================================================================
suppressMessages({library(data.table); library(GenomicRanges)})
source("scripts/matching.R")
G <- "reference/downloads"
SRC <- list(   # files downloaded by workflow/reference_tables.sh
  mouse  = list(ens = "mouse_ensembl.gtf.gz",    ref = "mouse_refseq.gtf.gz",  rep = "mouse_assembly_report.txt"),
  pig    = list(ens = "pig_ensembl.gtf.gz",      ref = "pig_refseq.gtf.gz",    rep = "pig_assembly_report.txt"),
  dog    = list(ens = "dog_ensembl_plus.gtf.gz", ref = "dog_refseq.gtf.gz",    rep = "dog_assembly_report.txt"),
  dog10k = list(ens = "dog10k_ensembl.gtf.gz",   ref = "dog10k_refseq.gtf.gz", rep = "dog10k_assembly_report.txt"))
CRIT <- c("any", "shorter", "recip")

read_exons <- function(files) {
  rbindlist(lapply(file.path(G, files), function(f) {
    cmd <- if (grepl("\\.gz$", f)) paste("gzip -dc", shQuote(f)) else paste("cat", shQuote(f))
    x <- fread(cmd = paste(cmd, "| awk -F'\\t' '$3==\"exon\"'"), sep = "\t", header = FALSE,
               select = c(1, 4, 5, 7, 9), col.names = c("chr", "start", "end", "strand", "attr"), quote = "")
    x[, gene_id := sub('^.*gene_id "([^"]*)".*$', "\\1", attr)][, attr := NULL]
  }))
}

refseq_chr_map <- function(f) {
  r <- fread(cmd = paste("grep -v '^#'", shQuote(file.path(G, f))), sep = "\t", header = FALSE, quote = "")
  setNames(ifelse(r$V2 == "assembled-molecule", sub("^chr", "", r$V3), r$V5), r$V7)
}

# union of exons per gene -> GRanges with gene_id; exonic length per gene
gene_exons <- function(x) {
  gr <- GRanges(x$chr, IRanges(x$start, x$end), strand = x$strand, gene_id = x$gene_id)
  red <- unlist(reduce(split(gr, gr$gene_id)))
  red$gene_id <- names(red); names(red) <- NULL
  red
}

match_pairs <- function(E, R) {
  h <- findOverlaps(E, R)                                   # strand-aware by default
  w <- width(pintersect(E[queryHits(h)], R[subjectHits(h)]))
  p <- data.table(e = E$gene_id[queryHits(h)], r = R$gene_id[subjectHits(h)], s = w)[, .(S = sum(s)), by = .(e, r)]
  LE <- tapply(width(E), E$gene_id, sum); LR <- tapply(width(R), R$gene_id, sum)
  p[, `:=`(LE = as.numeric(LE[e]), LR = as.numeric(LR[r]))]
  p[, `:=`(any = S >= 1, shorter = S >= 0.5 * pmin(LE, LR), recip = S >= 0.5 * LE & S >= 0.5 * LR)]
  p
}

# loci = connected components of the bipartite match graph (min-label propagation)
components <- function(e, r) {
  a <- paste0("E:", e); b <- paste0("R:", r)
  nodes <- unique(c(a, b)); ia <- match(a, nodes); ib <- match(b, nodes)
  lab <- seq_along(nodes)
  repeat {
    m <- pmin(lab[ia], lab[ib])
    new <- lab
    new <- pmin(new, tapply(m, ia, min)[as.character(seq_along(nodes))], na.rm = TRUE)
    new <- pmin(new, tapply(m, ib, min)[as.character(seq_along(nodes))], na.rm = TRUE)
    if (identical(new, lab)) break
    lab <- new
  }
  setNames(lab, nodes)
}

summarise_key <- function(key) {
  s <- SRC[[key]]
  ex <- read_exons(s$ens); rx <- read_exons(s$ref)
  ex[, gene_id := sub("\\.[0-9]+$", "", gene_id)]
  cm <- refseq_chr_map(s$rep)
  rx[, chr := ifelse(chr %in% names(cm), cm[chr], chr)]
  chr_ok <- mean(rx$chr %in% unique(ex$chr))
  E <- gene_exons(ex); R <- gene_exons(rx)
  P <- match_pairs(E, R)
  dE <- read_de(key, "ensembl"); dR <- read_de(key, "refseq")
  dE$gene_id <- sub("\\.[0-9]+$", "", dE$gene_id)
  idm <- compare_dbs(key, "ensembl", "refseq")              # Gene-ID matching, for comparison
  degE <- unique(dE$gene_id[dE$is_deg]); degR <- unique(dR$gene_id[dR$is_deg])
  idE <- unique(idm$dx$gene_id[idm$dx$is_deg & idm$cx]); idR <- unique(idm$dy$gene_id[idm$dy$is_deg & idm$cy])

  rows <- NULL; per_gene <- NULL
  for (cr in CRIT) {
    pm <- P[get(cr) == TRUE]
    hasE <- degE %in% pm$e; hasR <- degR %in% pm$r
    comp <- components(pm$e, pm$r)
    locE <- unique(comp[paste0("E:", degE[hasE])]); locR <- unique(comp[paste0("R:", degR[hasR])])
    rows <- rbind(rows, data.table(
      key, criterion = cr, chr_mapped = round(chr_ok, 4),
      DEG_Ens = length(degE), DEG_Ref = length(degR),
      gap_Ens = 1 - mean(hasE), gap_Ref = 1 - mean(hasR), gap = 1 - (mean(hasE) + mean(hasR)) / 2,
      J = length(intersect(locE, locR)) / length(union(locE, locR)),
      loci_shared = length(intersect(locE, locR)), loci_Ens_only = length(setdiff(locE, locR)),
      loci_Ref_only = length(setdiff(locR, locE)),
      # agreement of the two matching rules on DEGs: ID-matched but no coordinate match, and vice versa
      Ens_id_not_coord = sum(degE %in% idE & !hasE), Ens_coord_not_id = sum(!(degE %in% idE) & hasE),
      Ref_id_not_coord = sum(degR %in% idR & !hasR), Ref_coord_not_id = sum(!(degR %in% idR) & hasR)))
    per_gene <- rbind(per_gene,
      data.table(key, criterion = cr, db = "ensembl", gene_id = degE, id_match = degE %in% idE, coord_match = hasE),
      data.table(key, criterion = cr, db = "refseq",  gene_id = degR, id_match = degR %in% idR, coord_match = hasR))
  }
  rows[, `:=`(gap_id = idm$coverage, J_id = idm$J)]
  P[, organism := key]
  list(summary = rows, per_gene = per_gene, pairs = P)
}

# main: run only when called as a script (functions above can be sourced)
if (sys.nframe() == 0L) {
  dir.create("tables/coord_match", recursive = TRUE, showWarnings = FALSE)
  keys <- if (length(commandArgs(TRUE))) strsplit(commandArgs(TRUE)[1], ",")[[1]] else c("mouse", "pig", "dog", "dog10k")
  out <- lapply(keys, function(k) { cat("...", k, "\n"); summarise_key(k) })
  S <- rbindlist(lapply(out, `[[`, "summary"))
  num <- names(S)[vapply(S, is.double, TRUE)]
  S[, (num) := lapply(.SD, round, 3), .SDcols = num]
  print(S[, .(key, criterion, chr_mapped, DEG_Ens, DEG_Ref, gap_id, gap_Ens, gap_Ref, gap, J_id, J)])
  fwrite(S, "tables/coord_match/summary.tsv", sep = "\t")
  fwrite(rbindlist(lapply(out, `[[`, "per_gene")), "tables/coord_match/deg_matching.tsv.gz", sep = "\t")
  fwrite(rbindlist(lapply(out, `[[`, "pairs")), "tables/coord_match/pairs.tsv.gz", sep = "\t")
}
