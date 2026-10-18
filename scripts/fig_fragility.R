#!/usr/bin/env Rscript
# Fig. 6: fragility features: kruchy (fragile) vs stabilny (stable), across organisms.
# Grouped bars from each organism's summary.tsv (medians / percentages), faceted by
# feature. Shows the universal core (biotype + expression separate everywhere) vs
# contextual features (overlap, length-diff) that only separate in mouse/pig.
suppressMessages({library(ggplot2); library(tidyr); library(dplyr)})
dir.create("figures", showWarnings = FALSE)

orgs <- c(mouse="mouse", pig="pig", dog="dog")
rows <- lapply(names(orgs), function(pl) {
  s <- read.delim(file.path(paste0(orgs[[pl]],"_local"),"results",orgs[[pl]],"fragility","summary.tsv"))
  s$org <- pl; s
})
S <- bind_rows(rows)
S$org <- factor(S$org, levels = c("mouse","pig","dog"))
S$grupa <- factor(S$grupa, levels = c("kruchy","stabilny"))

feat <- tibble(
  key = c("pct_protein","med_ekspresja","frakcja_z_overlap","med_abs_dl_ratio"),
  lab = c("% protein-coding","expression (median log2; scale differs by DE method)",
          "fraction of genes with exon overlap","median absolute log2 length ratio (Ens/RefSeq)"))
long <- S %>% select(org, grupa, all_of(feat$key)) %>%
  pivot_longer(all_of(feat$key), names_to = "key", values_to = "val") %>%
  left_join(feat, by = "key")
long$lab <- factor(long$lab, levels = feat$lab)

pal <- c(kruchy = "#D85A30", stabilny = "#1D9E75")
p <- ggplot(long, aes(org, val, fill = grupa)) +
  geom_col(position = position_dodge(0.75), width = 0.68) +
  geom_text(aes(label = ifelse(grepl("%", lab), sprintf("%.0f", val), sprintf("%.2f", val))),
            position = position_dodge(0.75), vjust = -0.35, size = 3, color = "grey20") +
  facet_wrap(~lab, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = pal, name = NULL,
                    labels = c("fragile (DEG in Ensembl only)","stable (DEG in both)")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 10.5))

ggsave("figures/fig6_fragility.pdf", p, device = cairo_pdf, width = 8.2, height = 5.6)
ggsave("figures/fig6_fragility.png", p, width = 8.2, height = 5.6, dpi = 150)
cat("zapisano figures/fig6_fragility.{pdf,png}\n")
