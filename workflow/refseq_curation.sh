#!/usr/bin/env bash
# Annotation-maturity proxy independent of cross-references: share of RefSeq
# transcripts that are curated (NM_/NR_) vs predicted (XM_/XR_), and share of genes
# with >= 1 curated transcript, from the RefSeq GTFs downloaded by
# workflow/reference_tables.sh. Output: tables/refseq_curation.tsv
# Usage: bash workflow/refseq_curation.sh   (from the repository root)
set -euo pipefail
D=reference/downloads
mkdir -p tables; out=tables/refseq_curation.tsv
printf "organism\tannotation_source\ttranscripts\tNM\tNR\tXM\tXR\tcurated_tx_pct\tgenes\tgenes_with_curated_pct\n" > $out
for o in mouse pig dog dog10k; do
  gzip -dc "$D/${o}_refseq.gtf.gz" | awk -F'\t' -v o=$o '
    /^#!annotation-source/ { sub(/^#!annotation-source /,""); src=$0; next }
    $3=="exon" {   # exon lines: older GTFs (CanFam3.1) have no transcript lines
      match($9,/transcript_id "[^"]+"/); t=substr($9,RSTART+15,RLENGTH-16)
      match($9,/gene_id "[^"]+"/);       g=substr($9,RSTART+9,RLENGTH-10)
      ag[g]=1
      if (!(t in seen)) { seen[t]=1; p=substr(t,1,3); n[p]++; tot++
                          if (p=="NM_"||p=="NR_") cg[g]=1 } }
    END { c=n["NM_"]+n["NR_"]; ng=0; for (g in ag) ng++; nc=0; for (g in cg) nc++
          printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.1f\t%d\t%.1f\n", o,src,tot,
                 n["NM_"],n["NR_"],n["XM_"],n["XR_"],100*c/tot,ng,100*nc/ng }' >> $out
done
cat $out
