#!/bin/bash
# ============================================================================
# run_all.sh: reproduce all tables and figures of the paper from the count
# matrices (laptop part; the alignment and counting workflows are in workflow/).
#
# Inputs expected (see README):
#   <key>_local/results/<org>/counts/<database>/all_samples.counts.tsv
#     key mouse: ensembl gencode refseq ucsc | pig, dog: ensembl refseq ucsc |
#     dog10k (Dog10K_Boxer_Tasha): ensembl refseq
#   annotation tables and Gene ID cross-references: bash workflow/reference_tables.sh
#
# Usage: bash run_all.sh
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p figures tables

for f in mouse_local/results/mouse/counts/{ensembl,gencode,refseq,ucsc} \
         pig_local/results/pig/counts/{ensembl,refseq,ucsc} \
         dog_local/results/dog/counts/{ensembl,refseq,ucsc} \
         dog10k_local/results/dog/counts/{ensembl,refseq}; do
  [ -s "$f/all_samples.counts.tsv" ] || { echo "MISSING count matrix: $f/all_samples.counts.tsv"; exit 1; }
done
[ -s xref/ens_mouse_entrez.tsv.gz ] || { echo "Reference tables missing: run  bash workflow/reference_tables.sh"; exit 1; }

echo "== differential expression"
(cd mouse_local  && Rscript ../scripts/de_deseq2.R mouse)
(cd pig_local    && Rscript ../scripts/de_analysis.R pig 0 published)
(cd dog_local    && Rscript ../scripts/de_edger.R dog UDvsC 1.5)
(cd dog10k_local && Rscript ../scripts/de_edger.R dog UDvsC 1.5)

echo "== dog calibration: protein-coding genes only, as in the published results table (Dataset S1)"
C=calibration_dog/results/dog; rm -rf calibration_dog; mkdir -p $C/counts calibration_dog/reference
cp -r dog_local/results/dog/counts/ensembl $C/counts/; cp dog_local/results/dog/targets.tsv $C/
cp -r dog_local/reference/dog calibration_dog/reference/
awk -F'\t' '$2=="protein_coding"{print $1}' dog_local/reference/dog/annotations/ensembl_genes.tsv > calibration_dog/protein_coding.txt
(cd calibration_dog && Rscript ../scripts/de_edger.R dog UDvsC 1.5 protein_coding.txt > /dev/null)
Rscript -e 'd <- read.delim("calibration_dog/results/dog/de/ensembl_toptable.tsv")
  cat(sprintf("dog: %d DEGs (protein-coding, UD vs C; published 7619)\n", sum(d$adj.P.Val < 0.05 & abs(d$logFC) > 1.5, na.rm = TRUE)))' \
  | tee tables/calibration_dog.txt

echo "== fragile vs stable genes (Gene ID matching)"
(cd mouse_local  && Rscript ../scripts/fragility_features.R mouse 0   mouse)
(cd pig_local    && Rscript ../scripts/fragility_features.R pig   0   pig)
(cd dog_local    && Rscript ../scripts/fragility_features.R dog   1.5 dog)
(cd dog10k_local && Rscript ../scripts/fragility_features.R dog   1.5 dog10k)
for k in mouse_local/results/mouse pig_local/results/pig dog_local/results/dog dog10k_local/results/dog; do
  cp "$k/fragility/summary.tsv" "tables/fragility_summary_$(echo "$k" | cut -d_ -f1).tsv"
done

echo "== tables"
Rscript scripts/fragility_effects.R      # Table 3
Rscript scripts/summary_tables.R         # agreement numbers (Results, Figs 1, 2, 4), Table 4
Rscript scripts/calibration.R            # Table 1 (mouse, pig)
Rscript scripts/xref_coverage.R          # Gene ID cross-reference shares (Results)

echo "== figures"
for f in fig_heatmap fig_gradient fig_logfc fig_overlap fig_maturity fig_fragility; do
  Rscript "scripts/$f.R"
done
echo "done: figures/ and tables/"
