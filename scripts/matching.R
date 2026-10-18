# ============================================================================
# matching.R: shared gene matching between annotation databases by NCBI Gene ID.
# Sourced by the figure scripts, fragility_features.R and gene_matching.R.
#
#   Ensembl / GENCODE gene -> GeneID(s): Ensembl's EntrezGene cross-references for
#       the release used (xref/ens_<key>_entrez.tsv.gz; GENCODE ids lose their version)
#   RefSeq / UCSC gene -> GeneID: db_xref "GeneID:" of the RefSeq GTF
#       (xref/refseq_<key>_geneid.tsv; the UCSC ncbiRefSeq track reuses RefSeq gene ids)
#
# For a pair of databases X, Y:
#   comparable gene of X  = a gene whose GeneID(s) are annotated in Y (any gene of
#                           Y's GTF, whether or not it passed Y's expression filter)
#   G_X                   = GeneIDs (shared with Y) of the comparable DEGs of X
#   agreement  J          = |G_X ∩ G_Y| / |G_X ∪ G_Y|
#   coverage gap of X     = share of DEGs of X that are not comparable
# Keys: "mouse", "pig", "dog" (CanFam3.1 / Ensembl 92 / RefSeq AR105) and
# "dog10k" (Dog10K_Boxer_Tasha / Ensembl 112 / RefSeq AR106).
# ============================================================================

DE_CFG <- list(
  mouse  = list(dir = "mouse_local/results/mouse/de_deseq2",  suf = "_deseq2.tsv",   lfc = "log2FoldChange", padj = "padj",      L = 0,
                ann = "mouse_local/reference/mouse/annotations"),
  pig    = list(dir = "pig_local/results/pig/de",             suf = "_toptable.tsv", lfc = "logFC",          padj = "adj.P.Val", L = 0,
                ann = "pig_local/reference/pig/annotations"),
  dog    = list(dir = "dog_local/results/dog/de",             suf = "_toptable.tsv", lfc = "logFC",          padj = "adj.P.Val", L = 1.5,
                ann = "dog_local/reference/dog/annotations"),
  dog10k = list(dir = "dog10k_local/results/dog/de",          suf = "_toptable.tsv", lfc = "logFC",          padj = "adj.P.Val", L = 1.5,
                ann = "dog10k_local/reference/dog/annotations"))

# All GeneIDs annotated by a database (every gene of its GTF, not only tested genes)
db_universe <- function(key, db) {
  gm <- read.delim(file.path(DE_CFG[[key]]$ann, paste0(db, ".genemap.tsv")), header = FALSE,
                   colClasses = "character")
  unique(unlist(gene_ids(gm[[1]], db, key)))
}

.map_cache <- new.env()
load_maps <- function(key) {
  if (!is.null(.map_cache[[key]])) return(.map_cache[[key]])
  ex <- read.delim(gzfile(file.path("xref", paste0("ens_", key, "_entrez.tsv.gz"))), stringsAsFactors = FALSE)
  ex <- unique(ex[ex$db_name == "EntrezGene", c("gene_stable_id", "xref")])
  rx <- read.delim(file.path("xref", paste0("refseq_", key, "_geneid.tsv")), header = FALSE,
                   col.names = c("gene_id", "GeneID"), colClasses = "character")
  rx <- unique(rx[rx$GeneID != "", ])
  m <- list(ens = split(as.character(ex$xref), ex$gene_stable_id),
            ref = setNames(rx$GeneID, rx$gene_id))
  assign(key, m, envir = .map_cache); m
}

# GeneIDs of a vector of gene ids from database `db` (list, one element per gene)
gene_ids <- function(ids, db, key) {
  m <- load_maps(key)
  if (db %in% c("ensembl", "gencode")) {
    out <- m$ens[match(sub("\\.[0-9]+$", "", ids), names(m$ens))]   # hashed lookup
    out[vapply(out, is.null, TRUE)] <- list(character(0))
    unname(out)
  } else {
    v <- unname(m$ref[match(ids, names(m$ref))])
    lapply(v, function(x) if (is.na(x)) character(0) else x)
  }
}

read_de <- function(key, db) {
  c <- DE_CFG[[key]]
  d <- read.delim(file.path(c$dir, paste0(db, c$suf)), stringsAsFactors = FALSE)
  d$is_deg <- !is.na(d[[c$padj]]) & d[[c$padj]] < 0.05 & abs(d[[c$lfc]]) > c$L
  d$lfc_value <- d[[c$lfc]]
  d
}

# Compare two databases of one organism; returns agreement, coverage and the sets
compare_dbs <- function(key, x, y, dx = read_de(key, x), dy = read_de(key, y)) {
  ix <- gene_ids(dx$gene_id, x, key); iy <- gene_ids(dy$gene_id, y, key)
  Ux <- db_universe(key, x); Uy <- db_universe(key, y)      # GeneIDs annotated by each database
  hit <- function(l, U) { n <- lengths(l); f <- rep(seq_along(l), n)  # any GeneID of the gene in U
    out <- logical(length(l)); out[unique(f[unlist(l) %in% U])] <- TRUE; out }
  cx <- hit(ix, Uy); cy <- hit(iy, Ux)                      # comparable genes
  Gx <- intersect(unique(unlist(ix[dx$is_deg & cx])), Uy)
  Gy <- intersect(unique(unlist(iy[dy$is_deg & cy])), Ux)
  list(J = length(intersect(Gx, Gy)) / length(union(Gx, Gy)),
       gap_x = 1 - sum(dx$is_deg & cx) / sum(dx$is_deg),
       gap_y = 1 - sum(dy$is_deg & cy) / sum(dy$is_deg),
       coverage = (2 - sum(dx$is_deg & cx) / sum(dx$is_deg) - sum(dy$is_deg & cy) / sum(dy$is_deg)) / 2,
       n_deg_x = sum(dx$is_deg), n_deg_y = sum(dy$is_deg),
       shared = length(intersect(Gx, Gy)), only_x = length(setdiff(Gx, Gy)), only_y = length(setdiff(Gy, Gx)),
       ix = ix, iy = iy, cx = cx, cy = cy, Gy = Gy, dx = dx, dy = dy)
}
