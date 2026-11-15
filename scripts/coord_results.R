# ============================================================================
# coord_results.R: read the coordinate-based matching results (coord_match.R)
# in the same shape as compare_dbs() of matching.R, for the figure scripts.
#   coord_summary(key, crit) -> list(J, gap_x, gap_y, coverage, n_deg_x, n_deg_y,
#                                    shared, only_x, only_y, nc_x, nc_y)
# shared / only_x / only_y are counted in loci (connected components of matched
# genes), nc_x / nc_y = DEGs without a counterpart in the other database.
# Primary criterion: "shorter" (shared exonic bp >= 50% of the shorter gene).
# ============================================================================
.coord <- read.delim("tables/coord_match/summary.tsv", stringsAsFactors = FALSE)
.coord_deg <- read.delim(gzfile("tables/coord_match/deg_matching.tsv.gz"), stringsAsFactors = FALSE)
coord_summary <- function(key, crit = "shorter") {
  s <- .coord[.coord$key == key & .coord$criterion == crit, ]
  g <- .coord_deg[.coord_deg$key == key & .coord_deg$criterion == crit, ]
  nc_x <- sum(g$db == "ensembl" & !g$coord_match); nc_y <- sum(g$db == "refseq" & !g$coord_match)
  list(J = s$J, gap_x = nc_x / s$DEG_Ens, gap_y = nc_y / s$DEG_Ref,
       coverage = (nc_x / s$DEG_Ens + nc_y / s$DEG_Ref) / 2,
       n_deg_x = s$DEG_Ens, n_deg_y = s$DEG_Ref,
       shared = s$loci_shared, only_x = s$loci_Ens_only, only_y = s$loci_Ref_only,
       nc_x = nc_x, nc_y = nc_y)
}
