#!/usr/bin/env Rscript
# Fig. 2: annotation-divergence gradient across organisms, Ensembl vs RefSeq.
# Two quantities, both as fractions (0-1): agreement among comparable genes (Jaccard)
# and the share of DEGs without a counterpart in the other database. Each is shown
# for two matching rules: genomic coordinates (strand-specific exon overlap covering
# >= 50% of the shorter gene; coord_match.R; solid) and NCBI Gene ID (matching.R; dashed).
# Dog re-analysed on the current reference (canFam6, Ensembl 112 / RefSeq AR106): open symbols
# to the right of dog.
suppressMessages({library(ggplot2)})
dir.create("figures", showWarnings = FALSE)
source("scripts/matching.R"); source("scripts/coord_results.R")

keys <- c("mouse", "pig", "dog", "dog10k")
rows <- NULL
for (k in keys) {
  i <- compare_dbs(k, "ensembl", "refseq"); g <- coord_summary(k)
  rows <- rbind(rows, data.frame(key = k,
    metric = rep(c("agreement", "gap"), 2), rule = rep(c("coordinates", "Gene ID"), each = 2),
    val = c(g$J, g$coverage, i$J, i$coverage)))
}
write.table(rows, "figures/fig2_gradient_values.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(xtabs(round(val, 3) ~ key + interaction(metric, rule), rows))

lev_m <- c(agreement = "DEG agreement among comparable genes (Jaccard)",
           gap = "share of DEGs without a counterpart")
rows$metric <- factor(lev_m[rows$metric], levels = lev_m)
rows$rule <- factor(rows$rule, levels = c("coordinates", "Gene ID"),
                    labels = c("genes matched by genomic coordinates", "genes matched by NCBI Gene ID"))
orgs <- c("mouse", "pig", "dog")
main <- rows[rows$key %in% orgs, ]; main$org <- factor(main$key, levels = orgs)
d10 <- rows[rows$key == "dog10k", ]; d10$x <- 3 + 0.25
lab <- main[main$rule == levels(main$rule)[1], ]
lab$y <- lab$val + ifelse(grepl("agreement", lab$metric), 0.07, -0.07)
pal <- setNames(c("#1D9E75", "#D85A30"), lev_m)

p <- ggplot(main, aes(as.numeric(org), val, colour = metric)) +
  geom_line(aes(linetype = rule, group = interaction(metric, rule)), linewidth = 1) +
  geom_point(aes(shape = rule, fill = metric), size = 3.2) +
  geom_text(data = lab, aes(y = y, label = sprintf("%.2f", val)), size = 3.8, show.legend = FALSE) +
  geom_point(data = d10, aes(x = x, shape = rule), size = 3.2, fill = "white", stroke = 1.1, show.legend = FALSE) +
  geom_text(data = d10[d10$rule == levels(d10$rule)[1], ], aes(x = x, label = sprintf("%.2f", val)),
            hjust = -0.4, size = 3.3, show.legend = FALSE) +
  annotate("text", x = 3.25, y = -0.1, label = "open symbols: dog on canFam6\n(Ens 112 / RefSeq AR106)",
           size = 2.9, colour = "grey35", lineheight = 0.9) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_linetype_manual(values = c("solid", "22"), name = NULL) +
  scale_shape_manual(values = c(21, 24), name = NULL) +
  scale_y_continuous(limits = c(-0.16, 1.02), breaks = seq(0, 1, .2)) +
  scale_x_continuous(breaks = 1:3, labels = c("mouse\nGRCm39 / Ens 110", "pig\nSscrofa11.1 / Ens 104",
                                              "dog\nCanFam3.1 / Ens 92"), expand = expansion(add = c(0.4, 0.75))) +
  labs(x = "annotation maturity (decreasing from left to right)", y = "proportion (0-1)") +
  guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2), shape = guide_legend(order = 2)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", legend.direction = "vertical", legend.box = "horizontal",
        legend.text = element_text(size = 9), panel.grid.minor = element_blank())

ggsave("figures/fig2_gradient.pdf", p, device = cairo_pdf, width = 7.6, height = 5.2)
ggsave("figures/fig2_gradient.png", p, width = 7.6, height = 5.2, dpi = 150)
cat("written: figures/fig2_gradient.{pdf,png}\n")
