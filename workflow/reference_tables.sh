#!/bin/bash
# ============================================================================
# reference_tables.sh: build every annotation table the R scripts need, from
# public, version-pinned sources (nothing is shipped with the repository).
#
#   <org>_local/reference/<org>/annotations/
#       <db>.genemap.tsv         gene_id <TAB> gene name (gene_id if unnamed)   all databases
#       ensembl_genes.tsv        gene_id biotype chr start end strand          Ensembl
#       ensembl_ntx.tsv          gene_id number_of_transcripts                  Ensembl
#       ensembl_exons.tsv        gene_id chr start end                          Ensembl
#       refseq_exons.tsv         gene_id chr start end                          RefSeq
#   xref/ens_<key>_entrez.tsv.gz      Ensembl gene -> NCBI Gene ID (Ensembl's own xrefs)
#   xref/refseq_<key>_geneid.tsv      RefSeq gene_id -> NCBI Gene ID (GTF db_xref)
#   pig Ensembl 112 gene table + Gene ID xrefs (cross-reference check only)
#
# Annotation versions (Table 1 of the paper):
#   mouse  : Ensembl 110, GENCODE vM33, RefSeq RS_2024_02, UCSC mm39 ncbiRefSeq
#   pig    : Ensembl 104, RefSeq Annotation Release 106, UCSC susScr11 ncbiRefSeq
#   dog    : Ensembl 92 + gene models for AMH and HMGA2 (absent from release 92)
#            (data/dog_added_gene_models/amh_hmga2.gtf), RefSeq AR105, UCSC canFam3 ncbiRefSeq
#   dog10k : Dog10K_Boxer_Tasha, Ensembl 112, RefSeq AR106
# NOTE: UCSC does not version its ncbiRefSeq tracks; the files downloaded today
# may differ slightly from those used for the paper (2026-10).
#
# Usage: bash workflow/reference_tables.sh [project_root]     (needs curl, gzip, awk)
# ============================================================================
set -euo pipefail
ROOT=$(cd "${1:-$(dirname "$0")/..}" && pwd)
TMP=$ROOT/reference/downloads; mkdir -p "$TMP" "$ROOT/xref"
ENS=https://ftp.ensembl.org/pub
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases
UCSC=https://hgdownload.soe.ucsc.edu/goldenPath

fetch() { [ -s "$2" ] || curl -s -f -m 1800 -o "$2" "$1" || { echo "download failed: $1"; exit 1; }; }

# gene_id -> gene name (first record of each gene; gene_id when no name is given)
genemap() { gzcat "$1" | awk -F'\t' '/^#/{next}{gid=gn="";
    if (match($9,/gene_id "[^"]+"/))   gid=substr($9,RSTART+9,RLENGTH-10);
    if (match($9,/gene_name "[^"]+"/)) gn=substr($9,RSTART+11,RLENGTH-12);
    if (gn=="") gn=gid; if (gid!="" && !(gid in s)) { s[gid]=1; print gid"\t"gn } }' > "$2"; }

# Ensembl gene, transcript-count and exon tables
ensembl_tables() {   # $1 gtf.gz  $2 annotation dir
  gzcat "$1" | awk -F'\t' -v A="$2" '
    /^#/ {next}
    { gid=""; if (match($9,/gene_id "[^"]+"/)) gid=substr($9,RSTART+9,RLENGTH-10) }
    $3=="gene" { bt=""; if (match($9,/gene_biotype "[^"]+"/)) bt=substr($9,RSTART+14,RLENGTH-15);
                 print gid"\t"bt"\t"$1"\t"$4"\t"$5"\t"$7 > (A"/ensembl_genes.tsv") }
    $3=="transcript" { nt[gid]++ }
    $3=="exon" { print gid"\t"$1"\t"$4"\t"$5 > (A"/ensembl_exons.tsv") }
    END { for (g in nt) print g"\t"nt[g] > (A"/ensembl_ntx.tsv") }'
}

# RefSeq exons and gene_id -> GeneID; records with an empty gene_id are dropped
refseq_tables() {    # $1 gtf.gz  $2 annotation dir  $3 xref key
  gzcat "$1" | awk -F'\t' -v A="$2" -v X="$ROOT/xref/refseq_$3_geneid.tsv" '
    /^#/ {next}
    { gid=""; if (match($9,/gene_id "[^"]+"/)) gid=substr($9,RSTART+9,RLENGTH-10) }
    gid=="" {next}
    $3=="gene" { gi=""; if (match($9,/GeneID:[0-9]+/)) gi=substr($9,RSTART+7,RLENGTH-7); print gid"\t"gi > X }
    $3=="exon" { print gid"\t"$1"\t"$4"\t"$5 > (A"/refseq_exons.tsv") }'
}

build() {   # $1 key  $2 organism dir name  $3 Ensembl GTF url  $4 Ensembl entrez url  $5 RefSeq GTF url  $6 UCSC GTF url or -  $7 GENCODE GTF url or -
  local key=$1 org=$2; local A="$ROOT/${key}_local/reference/$org/annotations"; mkdir -p "$A"
  echo "== $key"
  fetch "$3" "$TMP/${key}_ensembl.gtf.gz"
  if [ "$key" = dog ]; then   # Ensembl 92 + AMH and HMGA2 gene models (absent from release 92)
    { gzcat "$TMP/dog_ensembl.gtf.gz"; grep -hv "^#" "$ROOT"/data/dog_added_gene_models/amh_hmga2.gtf; } | gzip > "$TMP/dog_ensembl_plus.gtf.gz"
    E="$TMP/dog_ensembl_plus.gtf.gz"
  else E="$TMP/${key}_ensembl.gtf.gz"; fi
  genemap "$E" "$A/ensembl.genemap.tsv"; ensembl_tables "$E" "$A"
  fetch "$4" "$ROOT/xref/ens_${key}_entrez.tsv.gz"
  fetch "$5" "$TMP/${key}_refseq.gtf.gz"
  genemap "$TMP/${key}_refseq.gtf.gz" "$A/refseq.genemap.tsv"; refseq_tables "$TMP/${key}_refseq.gtf.gz" "$A" "$key"
  if [ "$6" != - ]; then fetch "$6" "$TMP/${key}_ucsc.gtf.gz"; genemap "$TMP/${key}_ucsc.gtf.gz" "$A/ucsc.genemap.tsv"; fi
  if [ "$7" != - ]; then fetch "$7" "$TMP/${key}_gencode.gtf.gz"; genemap "$TMP/${key}_gencode.gtf.gz" "$A/gencode.genemap.tsv"; fi
  wc -l "$A"/*.tsv | tail -1
}

build mouse mouse \
  "$ENS/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz" \
  "$ENS/release-110/tsv/mus_musculus/Mus_musculus.GRCm39.110.entrez.tsv.gz" \
  "$NCBI/10090/GCF_000001635.27-RS_2024_02/GCF_000001635.27_GRCm39_genomic.gtf.gz" \
  "$UCSC/mm39/bigZips/genes/mm39.ncbiRefSeq.gtf.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.annotation.gtf.gz"
build pig pig \
  "$ENS/release-104/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.104.gtf.gz" \
  "$ENS/release-104/tsv/sus_scrofa/Sus_scrofa.Sscrofa11.1.104.entrez.tsv.gz" \
  "$NCBI/9823/106/GCF_000003025.6_Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_genomic.gtf.gz" \
  "$UCSC/susScr11/bigZips/genes/susScr11.ncbiRefSeq.gtf.gz" -
build dog dog \
  "$ENS/release-92/gtf/canis_familiaris/Canis_familiaris.CanFam3.1.92.gtf.gz" \
  "$ENS/release-92/tsv/canis_familiaris/Canis_familiaris.CanFam3.1.92.entrez.tsv.gz" \
  "$NCBI/9615/105/GCF_000002285.3_CanFam3.1/GCF_000002285.3_CanFam3.1_genomic.gtf.gz" \
  "$UCSC/canFam3/bigZips/genes/canFam3.ncbiRefSeq.gtf.gz" -
build dog10k dog \
  "$ENS/release-112/gtf/canis_lupus_familiarisboxer/Canis_lupus_familiarisboxer.Dog10K_Boxer_Tasha.112.gtf.gz" \
  "$ENS/release-112/tsv/canis_lupus_familiarisboxer/Canis_lupus_familiarisboxer.Dog10K_Boxer_Tasha.112.entrez.tsv.gz" \
  "$NCBI/9615/106/GCF_000002285.5_Dog10K_Boxer_Tasha/GCF_000002285.5_Dog10K_Boxer_Tasha_genomic.gtf.gz" - -
# pig, Ensembl 112 (same assembly): only the gene table and the Gene ID xrefs, used to
# check whether a newer Ensembl release is better cross-referenced (Results, dog section)
A="$ROOT/pig_local/reference/pig/annotations"
fetch "$ENS/release-112/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.112.gtf.gz" "$TMP/pig112_ensembl.gtf.gz"
gzcat "$TMP/pig112_ensembl.gtf.gz" | awk -F'\t' '$3=="gene"{gid=bt="";
    if (match($9,/gene_id "[^"]+"/)) gid=substr($9,RSTART+9,RLENGTH-10);
    if (match($9,/gene_biotype "[^"]+"/)) bt=substr($9,RSTART+14,RLENGTH-15); print gid"\t"bt}' > "$A/ensembl112_genes.tsv"
fetch "$ENS/release-112/tsv/sus_scrofa/Sus_scrofa.Sscrofa11.1.112.entrez.tsv.gz" "$ROOT/xref/ens_pig112_entrez.tsv.gz"
# NCBI assembly reports: RefSeq -> Ensembl sequence names for the coordinate-based gene
# matching (scripts/coord_match.R); the GTFs above are reused from $TMP
G=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF
fetch "$G/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_assembly_report.txt" "$TMP/mouse_assembly_report.txt"
fetch "$G/000/003/025/GCF_000003025.6_Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_assembly_report.txt" "$TMP/pig_assembly_report.txt"
fetch "$G/000/002/285/GCF_000002285.3_CanFam3.1/GCF_000002285.3_CanFam3.1_assembly_report.txt" "$TMP/dog_assembly_report.txt"
fetch "$G/000/002/285/GCF_000002285.5_Dog10K_Boxer_Tasha/GCF_000002285.5_Dog10K_Boxer_Tasha_assembly_report.txt" "$TMP/dog10k_assembly_report.txt"
echo "reference tables ready"
