#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --time=08:00:00
#SBATCH --partition=biology_standard
#SBATCH --account=pl0485-03
#SBATCH --job-name=mouse_calib
#SBATCH --output=/mnt/storage_5/scratch/pl0485-03/annotation_project/logs/mouse_calib_%j.log
#SBATCH --error=/mnt/storage_5/scratch/pl0485-03/annotation_project/logs/mouse_calib_%j.err
# ============================================================================
# Mouse calibration against Cheng et al. 2025 (Efcab4b KO liver, Novogene).
# The paper states STAR -> featureCounts -> DESeq2 on "Ensembl GRCm39" (cites
# Ensembl 2023, i.e. ~release 108-110) and TruSeq Stranded libraries; it does NOT
# give the release or the featureCounts options. Our first pass used -s 0 -O -M
# (as for pig and dog) and Ensembl 112 -> 60 DEG vs 69 published.
#
# This job re-counts the EXISTING STAR BAMs (no re-alignment) with every
# plausible combination, so the matching one can be picked on the laptop by
# comparing DESeq2 output with the published gene table:
#   Ensembl release : 108 109 110 112
#   strand          : -s 0 / 1 / 2
#   options         : "plain" (Novogene/featureCounts defaults: no -O, no -M)
#                     "OM"    (-O -M)
# Output: $PROJ/transfer/mouse_calib/ens<rel>_s<s>_<opt>.counts.tsv  (rsync to the laptop)
# Submit: sbatch --export=ALL ~/run_mouse_calib.sh
# ============================================================================
set -euo pipefail

PROJ=/mnt/storage_5/scratch/pl0485-03/annotation_project
BAMDIR=$PROJ/results/mouse/star
GTFDIR=$PROJ/reference/mouse/calib_gtf
OUT=$PROJ/transfer/mouse_calib   # home quota is full -> keep on scratch
mkdir -p "$GTFDIR" "$OUT"

FC=$(find "$PROJ/.snakemake/conda" -maxdepth 3 -name featureCounts 2>/dev/null | head -1)
echo "featureCounts: $FC"
[ -x "$FC" ] || { echo "ERROR: featureCounts not found in conda envs"; exit 1; }

SAMPLES="WT_1 WT_2 WT_3 KO_1 KO_2 KO_3"
BAMS=""
for s in $SAMPLES; do b="$BAMDIR/$s/Aligned.sortedByCoord.out.bam"; [ -f "$b" ] || { echo "MISSING BAM: $b"; exit 1; }; BAMS="$BAMS $b"; done

for rel in 108 109 110 112; do
  gtf="$GTFDIR/Mus_musculus.GRCm39.${rel}.gtf"
  [ -f "$gtf" ] || { wget -q -O "$gtf.gz" "https://ftp.ensembl.org/pub/release-${rel}/gtf/mus_musculus/Mus_musculus.GRCm39.${rel}.gtf.gz"; gunzip -f "$gtf.gz"; }
  for s in 0 1 2; do
    for opt in plain OM; do
      X=""; [ "$opt" = "OM" ] && X="-O -M"
      o="$OUT/ens${rel}_s${s}_${opt}"
      if ! "$FC" -p --countReadPairs -s $s $X -t exon -g gene_id -T 8 -a "$gtf" -o "$o.fc.txt" $BAMS > "$o.fc.log" 2>&1; then
        echo "featureCounts FAILED for ens$rel s$s $opt - last lines of its output:"; tail -15 "$o.fc.log"; continue
      fi
      awk -v S="$SAMPLES" 'BEGIN{n=split(S,a," ")} NR==1{next}
           NR==2{printf "Geneid"; for(i=1;i<=n;i++) printf "\t%s", a[i]; print ""; next}
           {printf "%s",$1; for(i=7;i<=NF;i++) printf "\t%s",$i; print ""}' "$o.fc.txt" > "$o.counts.tsv"
      cp "$o.fc.txt.summary" "$o.summary.txt"; rm -f "$o.fc.txt"
      echo "done ens$rel s$s $opt: Assigned $(awk '$1=="Assigned"{print $2}' "$o.summary.txt")"
    done
  done
  # gene_id -> gene_name for this release
  awk -F'\t' '$3=="gene"{gid=gn="";
    if(match($9,/gene_id "[^"]+"/))gid=substr($9,RSTART+9,RLENGTH-10);
    if(match($9,/gene_name "[^"]+"/))gn=substr($9,RSTART+11,RLENGTH-12);
    if(gn=="")gn=gid; print gid"\t"gn}' "$gtf" > "$OUT/ens${rel}.genemap.tsv"
done
echo "ALL DONE -> $OUT (rsync to the laptop)"
