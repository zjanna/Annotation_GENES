#!/usr/bin/env Rscript
# Fig. 2: annotation-divergence gradient across organisms.
# Two axes, both as fractions (0-1), one panel: concordance (Jaccard on genes
# present in both bases) DROPS, coverage gap (fraction of DEGs unmatchable) RISES,
# as annotation maturity decreases (mouse -> pig -> dog). Values computed from the DE tables.
suppressMessages({library(ggplot2); library(tidyr)})
dir.create("figures", showWarnings = FALSE)

# Ensembl vs RefSeq agreement (Jaccard on comparable genes) and coverage gap (share of
# DEGs without a counterpart in the other database), genes matched by NCBI Gene ID.
source("scripts/matching.R")
ens_ref <- function(key) { r <- compare_dbs(key, "ensembl", "refseq"); c(concord = r$J, coverage = r$coverage) }
m <- sapply(c("mouse","pig","dog","dog10k"), ens_ref)
print(round(m, 3))
d <- data.frame(
  org        = factor(c("mouse","pig","dog"), levels = c("mouse","pig","dog")),
  assembly   = c("GRCm39 / Ens 110","Sscrofa11.1 / Ens 104","CanFam3.1 / Ens 92"),
  concord    = m["concord", 1:3],
  coverage   = m["coverage", 1:3])
d10 <- data.frame(org = factor("dog", levels = levels(d$org)),
                  metryka = c("concord", "coverage"), val = m[, "dog10k"])
long <- pivot_longer(d, c(concord, coverage), names_to = "metryka", values_to = "val")
# labels above the concordance line, below the coverage line (avoid collision at "pies")
# label above the higher point of each organism, below the lower one
higher <- ave(long$val, long$org, FUN = function(v) v == max(v) & v != min(v))
long$lab_y <- long$val + ifelse(higher == 1 | (long$metryka == "concord" & ave(long$val, long$org, FUN = function(v) length(unique(v))) == 1), 0.075, -0.075)
lev2 <- c("DEG agreement among comparable genes (Jaccard)", "share of DEGs without a counterpart")
long$metryka <- factor(long$metryka, levels = c("concord","coverage"), labels = lev2)
d10$metryka  <- factor(d10$metryka,  levels = c("concord","coverage"), labels = lev2)
pal <- setNames(c("#1D9E75", "#D85A30"), lev2)

p <- ggplot(long, aes(org, val, color = metryka, group = metryka)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3.4) +
  geom_text(aes(y = lab_y, label = sprintf("%.2f", val)), size = 4, show.legend = FALSE) +
  # dog re-analysed on the current reference (canFam6, Ensembl 112 / RefSeq AR106): open symbols
  geom_point(data = d10, aes(x = as.numeric(org) + 0.22), shape = 21, fill = "white", size = 3.4, stroke = 1.2, show.legend = FALSE) +
  geom_text(data = d10, aes(x = as.numeric(org) + 0.22, y = val, label = sprintf("%.2f", val)),
            hjust = -0.35, size = 3.4, show.legend = FALSE) +
  annotate("text", x = 3.22, y = -0.08, label = "open: canFam6\n(Ens 112 / RefSeq AR106)",
           size = 3, colour = "grey35", lineheight = 0.9) +
  scale_color_manual(values = pal, name = NULL) +
  scale_y_continuous(limits = c(-0.14, 1.02), breaks = seq(0,1,.2)) +
  scale_x_discrete(labels = paste0(d$org, "\n", d$assembly), expand = expansion(add = c(0.4, 0.75))) +
  labs(x = "annotation maturity (decreasing from left to right)", y = "proportion (0-1)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top", legend.direction = "vertical",
        panel.grid.minor = element_blank())

ggsave("figures/fig2_gradient.pdf", p, device = cairo_pdf, width = 7.2, height = 4.8)
ggsave("figures/fig2_gradient.png", p, width = 7.2, height = 4.8, dpi = 150)
cat("zapisano figures/fig2_gradient.{pdf,png}\n")
