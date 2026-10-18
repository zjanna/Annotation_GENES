#!/usr/bin/env Rscript
# Fig. 5: share of DEGs without a cross-database counterpart vs an objective
# annotation-maturity proxy: % of tested Ensembl genes that carry a gene symbol
# (naming completeness), computed independently of the DE result and of RefSeq.
# Genes are matched between Ensembl and RefSeq by NCBI Gene ID (matching.R).
# Dog is shown on the anchor annotation (CanFam3.1, Ensembl 92 / RefSeq AR105, filled)
# and on the current reference (canFam6, Ensembl 112 / RefSeq AR106, open).
suppressMessages({library(ggplot2); library(ggrepel)})
dir.create("figures", showWarnings = FALSE)
source("scripts/matching.R")

named <- function(key) {
  d <- read_de(key, "ensembl")
  s <- ifelse(is.na(d$symbol) | d$symbol == "", d$gene_id, d$symbol)
  100 * mean(!grepl("^ENS[A-Z]*G[0-9]", s))
}
keys <- c(mouse = "mouse", pig = "pig", dog = "dog", dog10k = "dog10k")
d <- data.frame(key = keys,
                label = c("mouse", "pig", "dog (CanFam3.1)", "dog (canFam6)"),
                maturity = sapply(keys, named),
                gap = sapply(keys, function(k) compare_dbs(k, "ensembl", "refseq")$coverage),
                current = c(TRUE, TRUE, TRUE, FALSE))
print(d[, c("label", "maturity", "gap")], row.names = FALSE)

main <- d[d$current, ]
p <- ggplot(d, aes(maturity, gap)) +
  geom_line(data = main[order(main$maturity), ], aes(group = 1), colour = "grey75", linewidth = 0.6, linetype = 2) +
  geom_segment(data = data.frame(x = d$maturity[3], y = d$gap[3], xe = d$maturity[4], ye = d$gap[4]),
               aes(x = x, y = y, xend = xe, yend = ye), colour = "#D85A30", linewidth = 0.5,
               arrow = arrow(length = unit(2.2, "mm"), type = "closed")) +
  geom_point(aes(shape = current), size = 4, colour = "#D85A30", fill = "white", stroke = 1.2) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
  ggrepel::geom_text_repel(aes(label = sprintf("%s  %.2f", label, gap)), size = 3.8,
                           point.padding = 0.6, box.padding = 0.6, min.segment.length = Inf, seed = 1) +
  labs(x = "annotation maturity (% of tested Ensembl genes with a gene symbol)",
       y = "share of DEGs without a counterpart") +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.12))) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())
ggsave("figures/fig5_maturity.pdf", p, device = cairo_pdf, width = 6.4, height = 4.8)
ggsave("figures/fig5_maturity.png", p, width = 6.4, height = 4.8, dpi = 150)
cat("fig5 done\n")
