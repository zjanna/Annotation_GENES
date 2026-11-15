#!/usr/bin/env Rscript
# Fig. 6: stable and annotation-sensitive DEGs, with each database as the reference and
# genes matched by genomic coordinates (fragility_coord.R): stable / fragile, present /
# fragile, absent. Panels: % protein-coding and median expression per class; class sizes
# above the bars. Expression scale differs by DE method (compare within an organism).
suppressMessages({library(ggplot2)})
dir.create("figures", showWarnings = FALSE)

s <- read.delim("tables/fragility_coord/summary.tsv", stringsAsFactors = FALSE)
s$grp <- factor(paste0(s$organism, "\n", s$db, " ref."),
                levels = paste0(rep(c("mouse", "pig", "dog"), each = 2), "\n", c("Ensembl", "RefSeq"), " ref."))
s$class <- factor(s$class, levels = c("stable", "fragile, present", "fragile, absent"))
long <- rbind(data.frame(s[, c("grp", "class", "n")], panel = "% protein-coding", val = s$pct_protein_coding),
              data.frame(s[, c("grp", "class", "n")], panel = "median expression (log2; scale differs by DE method)",
                         val = s$median_expr))
long$panel <- factor(long$panel, levels = unique(long$panel))
pal <- c(stable = "#1D9E75", "fragile, present" = "#E8A58C", "fragile, absent" = "#D85A30")
p <- ggplot(long, aes(grp, val, fill = class)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = n), position = position_dodge(0.8), vjust = -0.4, size = 2.5, colour = "grey30") +
  facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = pal, name = NULL,
                    labels = c("stable (DEG in both)", "fragile, present (counterpart not DE)",
                               "fragile, absent (no counterpart)")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0))
ggsave("figures/fig6_fragility.pdf", p, device = cairo_pdf, width = 8.2, height = 6.4)
ggsave("figures/fig6_fragility.png", p, width = 8.2, height = 6.4, dpi = 150)
cat("written: figures/fig6_fragility.{pdf,png}\n")
