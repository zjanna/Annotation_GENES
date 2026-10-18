#!/usr/bin/env Rscript
# Fig. 4: decomposition of the Ensembl vs RefSeq DEG comparison into its two parts:
#   * DEGs without a counterpart in the other database (Ensembl-only / RefSeq-only genes)
#   * among comparable genes (matched by NCBI Gene ID): shared, lost under RefSeq, gained
# Shares of all DEGs called under either database (comparable genes counted once per GeneID).
suppressMessages({library(ggplot2)})
dir.create("figures", showWarnings = FALSE)
source("scripts/matching.R")

lev <- c("Ensembl DEG, no RefSeq counterpart", "Ensembl only (lost under RefSeq)", "shared",
         "RefSeq only (gained)", "RefSeq DEG, no Ensembl counterpart")
bars <- NULL
for (k in c("mouse", "pig", "dog", "dog10k")) {
  r <- compare_dbs(k, "ensembl", "refseq")
  n <- c(sum(r$dx$is_deg & !r$cx), r$only_x, r$shared, r$only_y, sum(r$dy$is_deg & !r$cy))
  bars <- rbind(bars, data.frame(key = k, cat = lev, n = n))
}
bars$org <- factor(c(mouse = "mouse", pig = "pig", dog = "dog\n(CanFam3.1)", dog10k = "dog\n(canFam6)")[bars$key],
                   levels = c("mouse", "pig", "dog\n(CanFam3.1)", "dog\n(canFam6)"))
bars$cat <- factor(bars$cat, levels = lev)
bars$pct <- ave(bars$n, bars$key, FUN = function(x) x / sum(x))
write.table(bars[, c("key", "cat", "n", "pct")], "figures/fig4_overlap_values.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(xtabs(round(100 * pct) ~ key + cat, bars))

pal <- c("Ensembl DEG, no RefSeq counterpart" = "#E8A58C", "Ensembl only (lost under RefSeq)" = "#D85A30",
         "shared" = "#1D9E75", "RefSeq only (gained)" = "#2C7FB8", "RefSeq DEG, no Ensembl counterpart" = "#9DC3E2")
p <- ggplot(bars, aes(org, n, fill = cat)) +
  geom_col(width = 0.7, position = "fill", colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(pct >= 0.04, sprintf("%.0f%%", pct * 100), "")),
            position = position_fill(vjust = 0.5), size = 3.2, colour = "grey10") +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
  guides(fill = guide_legend(ncol = 2)) +
  labs(x = NULL, y = "share of DEGs (Ensembl or RefSeq)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", legend.text = element_text(size = 9.5),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
ggsave("figures/fig4_overlap.pdf", p, device = cairo_pdf, width = 7.2, height = 5)
ggsave("figures/fig4_overlap.png", p, width = 7.2, height = 5, dpi = 150)
cat("fig4 done\n")
