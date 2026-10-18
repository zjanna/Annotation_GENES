#!/usr/bin/env Rscript
# Fig. 3: Ensembl vs RefSeq log2FC on genes matched one-to-one by NCBI Gene ID.
suppressMessages({library(ggplot2)})
dir.create("figures", showWarnings=FALSE)
source("scripts/matching.R")
dat<-NULL; labs<-c()
for(org in c("mouse","pig","dog")){
  e <- read_de(org, "ensembl"); r <- read_de(org, "refseq")
  ie <- gene_ids(e$gene_id, "ensembl", org); ir <- gene_ids(r$gene_id, "refseq", org)
  # one-to-one pairs: Ensembl genes with exactly one GeneID, GeneIDs used once on each side
  se <- data.frame(id = vapply(ie, function(v) if (length(v) == 1) v else NA_character_, ""), ens = e$lfc_value)
  sr <- data.frame(id = vapply(ir, function(v) if (length(v) == 1) v else NA_character_, ""), ref = r$lfc_value)
  se <- se[!is.na(se$id) & !(se$id %in% se$id[duplicated(se$id)]), ]
  sr <- sr[!is.na(sr$id) & !(sr$id %in% sr$id[duplicated(sr$id)]), ]
  m <- merge(se, sr, by = "id"); m <- m[is.finite(m$ens) & is.finite(m$ref), ]
  labs[org]<-sprintf("%s  (r = %.2f, n = %s)", org, cor(m$ens,m$ref), format(nrow(m), big.mark = ","))
  m$org<-org; dat<-rbind(dat,m)
}
print(labs)
dat$org<-factor(dat$org, levels=c("mouse","pig","dog"), labels=labs[c("mouse","pig","dog")])
p<-ggplot(dat,aes(ens,ref))+
  geom_abline(slope=1,intercept=0,color="grey65",linetype=2)+
  geom_point(alpha=0.22,size=0.45,color="#2C7FB8")+
  facet_wrap(~org,scales="free")+
  labs(x=expression(log[2]*"FC (Ensembl)"), y=expression(log[2]*"FC (RefSeq)"))+
  theme_minimal(base_size=12)+
  theme(strip.text=element_text(face="bold"),panel.grid.minor=element_blank())
ggsave("figures/fig3_logfc.pdf", p, device = cairo_pdf, width=9,height=3.4); ggsave("figures/fig3_logfc.png",p,width=9,height=3.4,dpi=600)
cat("fig3 done\n")
