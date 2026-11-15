#!/usr/bin/env Rscript
# Fig. 5: share of DEGs without a cross-database counterpart vs an annotation-maturity
# proxy: % of tested Ensembl genes that carry a gene symbol (naming completeness),
# computed independently of the DE result and of RefSeq.
# The share is shown for two matching rules: genomic coordinates (coord_match.R; filled
# circles, the primary measure, which does not use cross-references) and NCBI Gene ID
# (matching.R; triangles). Dog is shown on the anchor annotation (CanFam3.1, Ensembl 92 /
# RefSeq AR105) and on the current reference (canFam6, Ensembl 112 / RefSeq AR106, open).
suppressMessages({library(ggplot2); library(ggrepel)})
dir.create("figures", showWarnings = FALSE)
source("scripts/matching.R"); source("scripts/coord_results.R")

named <- function(key) {
  d <- read_de(key, "ensembl")
  s <- ifelse(is.na(d$symbol) | d$symbol == "", d$gene_id, d$symbol)
  100 * mean(!grepl("^ENS[A-Z]*G[0-9]", s))
}
keys <- c(mouse = "mouse", pig = "pig", dog = "dog", dog10k = "dog10k")
lab <- c("mouse", "pig", "dog (CanFam3.1)", "dog (canFam6)")
mat <- sapply(keys, named)
d <- rbind(
  data.frame(key = keys, label = lab, maturity = mat, rule = "genomic coordinates",
             gap = sapply(keys, function(k) coord_summary(k)$coverage)),
  data.frame(key = keys, label = lab, maturity = mat, rule = "NCBI Gene ID",
             gap = sapply(keys, function(k) compare_dbs(k, "ensembl", "refseq")$coverage)))
d$current <- d$key != "dog10k"
d$rule <- factor(d$rule, levels = c("genomic coordinates", "NCBI Gene ID"))
print(d[, c("label", "rule", "maturity", "gap")], row.names = FALSE)
write.table(d, "figures/fig5_maturity_values.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

main <- d[d$current, ]; main <- main[order(main$rule, main$maturity), ]
arr <- do.call(rbind, lapply(split(d, d$rule), function(z)
  data.frame(rule = z$rule[1], x = z$maturity[z$key == "dog"], y = z$gap[z$key == "dog"],
             xe = z$maturity[z$key == "dog10k"], ye = z$gap[z$key == "dog10k"])))
col <- c("genomic coordinates" = "#D85A30", "NCBI Gene ID" = "grey55")
p <- ggplot(d, aes(maturity, gap, colour = rule)) +
  geom_line(data = main, aes(group = rule, linetype = rule), linewidth = 0.6) +
  geom_segment(data = arr, aes(x = x, y = y, xend = xe, yend = ye), linewidth = 0.45,
               arrow = arrow(length = unit(2, "mm"), type = "closed"), show.legend = FALSE) +
  geom_point(data = d[d$current, ], aes(shape = rule, fill = rule), size = 3.6, stroke = 1.1) +
  geom_point(data = d[!d$current, ], aes(shape = rule), fill = "white", size = 3.6, stroke = 1.1, show.legend = FALSE) +
  ggrepel::geom_text_repel(aes(label = sprintf("%s  %.2f", label, gap)), size = 3.3, show.legend = FALSE,
                           point.padding = 0.5, box.padding = 0.5, min.segment.length = 0.3, seed = 1) +
  scale_colour_manual(values = col, name = "genes matched by") +
  scale_fill_manual(values = col, name = "genes matched by") +
  scale_shape_manual(values = c(21, 24), name = "genes matched by") +
  scale_linetype_manual(values = c("solid", "22"), name = "genes matched by") +
  labs(x = "annotation maturity (% of tested Ensembl genes with a gene symbol)",
       y = "share of DEGs without a counterpart") +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.12))) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave("figures/fig5_maturity.pdf", p, device = cairo_pdf, width = 6.8, height = 5)
ggsave("figures/fig5_maturity.png", p, width = 6.8, height = 5, dpi = 150)
cat("fig5 done\n")
