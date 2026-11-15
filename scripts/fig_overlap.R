#!/usr/bin/env Rscript
# Fig. 4: decomposition of the Ensembl vs RefSeq DEG comparison into its two parts:
#   * DEGs without a counterpart in the other database (Ensembl-only / RefSeq-only genes)
#   * among comparable genes: shared, lost under RefSeq, gained under RefSeq
# Panel A: genes matched by genomic coordinates (coord_match.R; comparable genes counted
# once per locus). Panel B: genes matched by NCBI Gene ID (matching.R; counted once per
# Gene ID). Shares of all DEGs called under either database. Every segment is labelled.
suppressMessages({library(ggplot2); library(ggrepel)})
dir.create("figures", showWarnings = FALSE)
source("scripts/matching.R"); source("scripts/coord_results.R")

lev <- c("Ensembl DEG, no RefSeq counterpart", "Ensembl only (lost under RefSeq)", "shared",
         "RefSeq only (gained)", "RefSeq DEG, no Ensembl counterpart")
rules <- c(coord = "A  genes matched by genomic coordinates", id = "B  genes matched by NCBI Gene ID")
bars <- NULL
for (k in c("mouse", "pig", "dog", "dog10k")) {
  g <- coord_summary(k)
  r <- compare_dbs(k, "ensembl", "refseq")
  bars <- rbind(bars,
    data.frame(key = k, rule = "coord", cat = lev, n = c(g$nc_x, g$only_x, g$shared, g$only_y, g$nc_y)),
    data.frame(key = k, rule = "id", cat = lev,
               n = c(sum(r$dx$is_deg & !r$cx), r$only_x, r$shared, r$only_y, sum(r$dy$is_deg & !r$cy))))
}
bars$org <- factor(c(mouse = "mouse", pig = "pig", dog = "dog\n(CanFam3.1)", dog10k = "dog\n(canFam6)")[bars$key],
                   levels = c("mouse", "pig", "dog\n(CanFam3.1)", "dog\n(canFam6)"))
bars$cat <- factor(bars$cat, levels = lev)
bars$rule <- factor(rules[bars$rule], levels = rules)
bars$pct <- ave(bars$n, bars$key, bars$rule, FUN = function(x) x / sum(x))
write.table(bars[, c("key", "rule", "cat", "n", "pct")], "figures/fig4_overlap_values.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
print(xtabs(round(100 * pct) ~ key + cat + rule, bars))

# label positions: centre of each stacked segment (stack order = reverse factor order)
bars <- bars[order(bars$rule, bars$org, -as.integer(bars$cat)), ]
bars$top <- ave(bars$pct, bars$rule, bars$org, FUN = cumsum)
bars$mid <- bars$top - bars$pct / 2
bars$lab <- ifelse(bars$pct < 0.005 & bars$n > 0, "<1%", sprintf("%.0f%%", bars$pct * 100))
big <- bars[bars$pct >= 0.05, ]; small <- bars[bars$pct < 0.05 & bars$n > 0, ]

pal <- c("Ensembl DEG, no RefSeq counterpart" = "#E8A58C", "Ensembl only (lost under RefSeq)" = "#D85A30",
         "shared" = "#1D9E75", "RefSeq only (gained)" = "#2C7FB8", "RefSeq DEG, no Ensembl counterpart" = "#9DC3E2")
p <- ggplot(bars, aes(org, pct, fill = cat)) +
  geom_col(width = 0.6, position = position_stack(reverse = FALSE), colour = "white", linewidth = 0.3) +
  geom_text(data = big, aes(y = mid, label = lab), size = 3, colour = "grey10") +
  geom_text_repel(data = small, aes(y = mid, label = lab, colour = cat), size = 2.8, nudge_x = 0.45,
                  direction = "y", hjust = 0, segment.size = 0.3, min.segment.length = 0,
                  box.padding = 0.12, show.legend = FALSE, seed = 1) +
  facet_wrap(~ rule, ncol = 1) +
  scale_fill_manual(values = pal, name = NULL) +
  scale_colour_manual(values = c(pal[1:2], shared = "#137A5A", pal[4:5]), guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"), expand = expansion(mult = c(0, 0.02))) +
  scale_x_discrete(expand = expansion(add = c(0.5, 0.9))) +
  guides(fill = guide_legend(ncol = 2)) +
  labs(x = NULL, y = "share of DEGs (Ensembl or RefSeq)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", legend.text = element_text(size = 9), strip.text = element_text(hjust = 0, face = "bold"),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
ggsave("figures/fig4_overlap.pdf", p, device = cairo_pdf, width = 7.4, height = 8)
ggsave("figures/fig4_overlap.png", p, width = 7.4, height = 8, dpi = 150)
cat("fig4 done\n")
