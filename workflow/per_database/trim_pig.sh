#!/bin/bash
# ============================================================================
# trim_pig.sh: raw pig reads from SRA (PRJNA1147879) -> the adapter- and
# quality-trimmed FASTQ read by Snakefile_pig (data/pig/fastq/RNA_<i>.{1,2}.trimmed.fastq.gz).
# BBDuk from BBMap v36.32, the version used by Woźniak et al. 2025 (adapter removal
# with BBDuk's standard adapter-trimming settings, quality trimming at Phred 5);
# RNA_<i> is library S3797Nr<i> (run accessions: ../../data/sample_runs.tsv).
#
# Usage (from workflow/per_database):  BBDUK=/path/to/bbmap/bbduk.sh bash trim_pig.sh
# Needs: sra-tools (prefetch, fasterq-dump), gzip, BBMap (bbduk.sh and resources/adapters.fa)
# ============================================================================
set -euo pipefail
BBDUK=${BBDUK:-bbduk.sh}
ADAPTERS=${ADAPTERS:-$(dirname "$(command -v "$BBDUK")")/resources/adapters.fa}
RAW=data/pig/raw; OUT=data/pig/fastq; QC=data/pig/trimming_qc
mkdir -p "$RAW" "$OUT" "$QC"

awk -F'\t' '$1=="pig"{print $2"\t"$3}' ../../data/sample_runs.tsv | while IFS=$'\t' read -r s run; do
  i=${s#RNA_}
  r1=$RAW/RNA_S3797Nr$i.1.fastq.gz; r2=$RAW/RNA_S3797Nr$i.2.fastq.gz
  if [ ! -s "$r1" ]; then
    prefetch "$run" -O "$RAW"
    fasterq-dump --split-files -O "$RAW" "$RAW/$run"
    gzip -c "$RAW/${run}_1.fastq" > "$r1"; gzip -c "$RAW/${run}_2.fastq" > "$r2"
    rm -rf "$RAW/${run}_1.fastq" "$RAW/${run}_2.fastq" "$RAW/$run"
  fi
  "$BBDUK" threads=16 in="$r1" in2="$r2" \
    out="$OUT/RNA_$i.1.trimmed.fastq.gz" out2="$OUT/RNA_$i.2.trimmed.fastq.gz" \
    qtrim=w trimq=5 ref="$ADAPTERS" ktrim=r k=23 mink=11 hdist=1 tbo tpe \
    minlength=50 removeifeitherbad=t overwrite=t \
    stats="$QC/${i}_adapters.txt" >> "$QC/${i}_trimming.txt" 2>&1
done
