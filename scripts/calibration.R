#!/usr/bin/env Rscript
# ============================================================================
# calibration.R: reproduction of the anchor studies (Table 1; run from the root).
#   mouse : our Ensembl DESeq2 result vs the 47 DEGs with values in the tables of
#           Cheng et al. 2025 (log2FC as published; sign convention KO vs WT)
#   pig   : DEGs at FDR < 0.05 (published list: 234) and with |log2FC| > 1 (59)
#   dog   : protein-coding genes only (published Dataset S1): done in run_all.sh
#           (tables/calibration_dog.txt)
# Writes tables/calibration_mouse.tsv
# ============================================================================
dir.create("tables", showWarnings = FALSE)
pub <- data.frame(
  gene = c("Mup9","Pck1","Cyp3a11","Mup17","Aacs","Rsrp1","Hba-a2","Gpr146","B4galt1","Onecut1",
           "9030619P08Rik","Gm13775","Acnat2","Mup-ps13","Hspa1b","Tat","Atoh8","Ptgds","Cyp3a59","Gnat1",
           "Ces2c","Loxl4","Serpina4-ps1","Mvk","Gm49395","E2f8","Mup7","Dio1","Xbp1","Fgfr1",
           "Socs2","Hes1","Sik1","Vgll4","Idi1","Prss8","Odc1","Pmvk","Igfbp6","Cd163","Ly96","Cd44",
           "Tsku","Casp8","Baiap2","Slc25a47","Hspa8"),
  lfc_published = c(-2.6075,-1.4881,-1.6837,-1.4155,1.7884,-1.1571,1.1723,1.0139,0.9729,1.9776,
                    -2.3037,-1.9459,-1.8897,-1.0885,-1.2885,-0.6794,-0.9598,-2.0571,-1.4336,0.8614,
                    -0.9825,1.1242,0.7744,0.7144,-1.4180,1.4379,-0.8927,-1.0201,0.6627,1.2086,
                    2.24,1.16,1.13,0.92,0.92,0.87,0.84,0.67,-1.39,-1.18,-1.17,-1.1,-1.07,-0.72,-0.71,-0.68,-0.67))
d <- read.delim("mouse_local/results/mouse/de_deseq2/ensembl_deseq2.tsv")
k <- d[match(toupper(pub$gene), toupper(d$symbol)), ]
pub$lfc_ours <- round(k$log2FoldChange, 3); pub$padj_ours <- signif(k$padj, 3)
ok <- !is.na(pub$lfc_ours)
cat(sprintf("mouse: %d DEGs (published 69) | %d of %d matched published DEGs recovered | same sign %d | median |dLFC| %.2f | r %.2f\n",
            sum(d$padj < 0.05, na.rm = TRUE), sum(pub$padj_ours[ok] < 0.05), sum(ok),
            sum(sign(pub$lfc_ours[ok]) == sign(pub$lfc_published[ok])),
            median(abs(pub$lfc_ours[ok] - pub$lfc_published[ok])), cor(pub$lfc_ours[ok], pub$lfc_published[ok])))
write.table(pub, "tables/calibration_mouse.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

p <- read.delim("pig_local/results/pig/de/ensembl_toptable.tsv")
cat(sprintf("pig: %d DEGs at FDR<0.05 (published 234); %d with |log2FC|>1 (published 59)\n",
            sum(p$adj.P.Val < 0.05, na.rm = TRUE), sum(p$adj.P.Val < 0.05 & abs(p$logFC) > 1, na.rm = TRUE)))
