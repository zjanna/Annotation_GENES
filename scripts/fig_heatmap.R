#!/usr/bin/env Rscript
# Fig. 1: per-organism DEG-list agreement between annotation databases.
# Pairwise Jaccard on comparable genes matched by NCBI Gene ID (matching.R), each unordered
# pair shown once (the matrix is symmetric and its diagonal is 1 by definition).
# Default: dot plot, one point per pair, coloured and shaped by comparison type
# (core Ensembl-RefSeq / within-family control / cross-family), pairs present in
# all three organisms joined across the maturity gradient.
# `Rscript fig_heatmap.R heatmap` writes the alternative lower-triangle heatmap.
suppressMessages({library(ggplot2); library(dplyr)})
dir.create("figures", showWarnings = FALSE)

source("scripts/matching.R")
anns <- list(mouse=c("ensembl","gencode","refseq","ucsc"), pig=c("ensembl","refseq","ucsc"), dog=c("ensembl","refseq","ucsc"))
nice <- c(ensembl="Ensembl", gencode="GENCODE", refseq="RefSeq", ucsc="UCSC")
tiles <- NULL
for (org in names(anns)) {
  d <- anns[[org]]
  for (i in seq_along(d)) for (j in seq_along(d)) if (i < j) {
    r <- compare_dbs(org, d[i], d[j])
    tiles <- rbind(tiles, data.frame(org=org, a=nice[d[i]], b=nice[d[j]], J=r$J))
  }
}
tiles$org <- factor(tiles$org, levels=c("mouse","pig","dog"))
lev <- c("Ensembl","GENCODE","RefSeq","UCSC")
# keep each unordered pair once, drop the diagonal
pairs <- tiles
pairs$pair <- paste0(pairs$a, "-", pairs$b)
fam <- c(Ensembl="E", GENCODE="E", RefSeq="R", UCSC="R")
pairs$type <- ifelse(pairs$a=="Ensembl" & pairs$b=="RefSeq", "core: Ensembl vs RefSeq",
              ifelse(fam[as.character(pairs$a)]==fam[as.character(pairs$b)], "control: same family",
                     "cross-family"))
pairs$type <- factor(pairs$type, levels=c("core: Ensembl vs RefSeq","control: same family","cross-family"))
write.table(pairs[,c("org","pair","type","J")], "figures/fig1_jaccard_values.tsv", sep="\t", quote=FALSE, row.names=FALSE)
print(pairs[,c("org","pair","type","J")], row.names=FALSE)

mode <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "dots"
if (mode == "heatmap") {
  pairs$a <- factor(pairs$a, levels=lev); pairs$b <- factor(pairs$b, levels=rev(lev))
  p <- ggplot(pairs, aes(a, b, fill=J)) +
    geom_tile(color="white", linewidth=1) +
    geom_text(aes(label=sprintf("%.2f", J)), size=3.6, color=ifelse(pairs$J>0.6,"white","grey15")) +
    facet_wrap(~org, scales="free", nrow=1) +
    scale_fill_gradient(low="#FDE9DE", high="#1D5C45", limits=c(0,1), name="Jaccard") +
    labs(x=NULL, y=NULL) + theme_minimal(base_size=12) +
    theme(panel.grid=element_blank(), axis.text.x=element_text(angle=30, hjust=1),
          strip.text=element_text(face="bold", size=11))
  w <- 8; h <- 3
} else {
  typelab <- c("core: Ensembl vs RefSeq" = "core: Ensembl vs RefSeq",
               "control: same family"    = "same family (control)",
               "cross-family"            = "cross-family")
  pal <- setNames(c("#D85A30", "#1D9E75", "#2C7FB8"), names(typelab))
  shp <- setNames(c(21, 24, 22), names(typelab))          # fillable shapes -> white ring
  # pairs observed in all three organisms are joined across the gradient
  common <- names(which(table(pairs$pair) == 3))
  dx <- c("Ensembl-RefSeq"=0, "RefSeq-UCSC"=0, "Ensembl-UCSC"=0.09,
          "Ensembl-GENCODE"=0, "GENCODE-RefSeq"=-0.09, "GENCODE-UCSC"=0.18)
  pairs$x <- as.numeric(pairs$org) + dx[pairs$pair]
  ln <- pairs[pairs$pair %in% common, ]
  labR <- pairs[pairs$org == "dog", ]
  labR$txt <- sprintf("%s  %.2f", labR$pair, labR$J)       # line name + end value
  labL <- pairs[pairs$org == "mouse", ]
  labL$txt <- sprintf("%s  %.2f", labL$pair, labL$J)
  xlabs <- c("mouse\nGRCm39 / Ens 110", "pig\nSscrofa11.1 / Ens 104", "dog\nCanFam3.1 / Ens 92")
  p <- ggplot(pairs, aes(x, J)) +
    geom_line(data=ln, aes(group=pair, colour=type, linewidth=type == "core: Ensembl vs RefSeq"), alpha=0.75) +
    geom_point(aes(shape=type, fill=type), colour="white", size=3.6, stroke=0.9) +
    ggrepel::geom_text_repel(data=labR, aes(label=txt), size=3.1, colour="grey20",
                             nudge_x=0.14, direction="y", hjust=0, min.segment.length=Inf, seed=1) +
    ggrepel::geom_text_repel(data=labL, aes(label=txt), size=3.1, colour="grey20",
                             nudge_x=-0.16, direction="y", hjust=1, box.padding=0.35,
                             min.segment.length=0.2, segment.colour="grey75", segment.size=0.3, seed=2) +
    # middle organism: values printed next to the points (core below, others above)
    geom_text(data=transform(pairs[pairs$org == "pig", ],
                             ny=ifelse(type == "core: Ensembl vs RefSeq", -0.028, 0.028)),
              aes(y=J + ny, label=sprintf("%.2f", J)), size=3.1, colour="grey20") +
    scale_colour_manual(values=pal, labels=typelab, name=NULL) +
    scale_fill_manual(values=pal, labels=typelab, name=NULL) +
    scale_shape_manual(values=shp, labels=typelab, name=NULL) +
    scale_linewidth_manual(values=c(`TRUE`=1.4, `FALSE`=0.8), guide="none") +
    guides(colour="none", fill=guide_legend(override.aes=list(size=3.6))) +
    scale_x_continuous(breaks=1:3, labels=xlabs, limits=c(0.05, 3.95), expand=c(0, 0)) +
    scale_y_continuous(limits=c(0.6, 1.04), breaks=seq(0.6, 1, 0.1)) +
    labs(x="annotation maturity (decreasing from left to right)",
         y="DEG-list agreement (Jaccard)") +
    theme_minimal(base_size=12) +
    theme(legend.position="top", legend.text=element_text(size=10),
          panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
          axis.text.x=element_text(size=10, colour="grey20"),
          axis.title.x=element_text(size=10.5, colour="grey30", margin=margin(t=6)))
  w <- 8; h <- 5
}
ggsave("figures/fig1_jaccard.pdf", p, device = cairo_pdf, width=w, height=h)
ggsave("figures/fig1_jaccard.png", p, width=w, height=h, dpi=150)
cat("zapisano figures/fig1_jaccard.{pdf,png}\n")
