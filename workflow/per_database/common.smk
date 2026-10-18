# -*- coding: utf-8 -*-
# ============================================================================
# common.smk: shared engine for the annotation-comparison pipeline
# ----------------------------------------------------------------------------
# This file is `include:`-d by each per-organism Snakefile (Snakefile_mouse,
# Snakefile_pig, Snakefile_dog). It contains ONLY the organism-agnostic logic:
# download references, harmonise chromosome names, build STAR index(es), align,
# count reads with featureCounts ONCE PER ANNOTATION, merge, extract gene
# features, and run the DE + fragility R scripts.
#
# The including Snakefile MUST load `configfile: "config/pipeline.yaml"` and
# define, BEFORE `include: "common.smk"`:
#   ORGANISM         str   short label, e.g. "pig" (used in all paths)
#   REF              str   path to the (decompressed) genome FASTA
#   GENOME_URL       str   URL of the gzipped genome FASTA to download
#   ANNOTATIONS      dict  {name: {"url":..., "chrom_fix":..., ["report_url":...]}}
#   ENSEMBL_KEY      str   which ANNOTATIONS key is the Ensembl anchor
#   SAMPLES          list  sample names
#   SJDB_OVERHANG    int   read_length - 1 (STAR splice-junction DB overhang)
#   STAR_EXTRA       str   extra STAR flags for this organism
#   (optional, with defaults:)
#   GENOMES          dict  {annotation: {"path":..., "url":...}}, each annotation's
#                          OWN genome FASTA (native chromosome names); default:
#                          every annotation uses REF (Ensembl names, GTFs harmonised)
#   SJDB_AT_MAPPING  bool  True: genome-only STAR index, splice junctions inserted
#                          from the GTF at mapping time (--sjdbGTFfile in STAR
#                          alignReads); default False: GTF built into the index
#   STAR_ENV         str   conda env with the STAR version (default envs/star.yaml)
#   FC_STRAND        int   featureCounts -s (0 unstranded, 1 stranded, 2 reverse); default 0
#   FC_MULTI         str   featureCounts multi-overlap/multi-mapping flags; default "-O -M"
#                          (mouse: "", featureCounts defaults, as the final mouse counts)
#   PREBUILT_INDEX   dict  {annotation: path}, use an existing STAR index (e.g. the
#                          previously built one) instead of building it
#   ENA_FASTQ        dict  {sample: [(url_R1, md5_R1), (url_R2, md5_R2)]}, samples
#                          downloaded as ready .fastq.gz from ENA (rule download_fastq)
#   COUNTS_ONLY      bool  True: stop at the count matrices (no DE / fragility; then
#                          TARGETS_FILE, DE and FRAGILITY_LFC are not needed); default False
#   ENA_SEP          str   mate separator in the downloaded file names: "_" ->
#                          {sample}_1.fastq.gz (default), "." -> {sample}.1.fastq.gz
#   TARGETS_FILE     str   sample sheet in the repo (copied to <results>/targets.tsv)
#   DE               dict  {"script":..., "args":..., "subdir":..., "suffix":...}
#   FRAGILITY_LFC    num   |log2FC| threshold of the anchor paper (fragility step)
#   get_reads(sample)      function -> {"r1": path, "r2": path} (paired FASTQ)
#
# STAR INDEX MODE (config["star_index_mode"]):
#   per_annotation – each annotation gets its OWN STAR index (its GTF as sjdb),
#                    its own alignment and its own featureCounts. The annotation
#                    therefore acts on BOTH alignment and counting, as it would
#                    for a user who picks that database for the whole analysis.
#   ensembl        – original design: one STAR index with the fixed Ensembl GTF
#                    (one alignment for all annotations); the
#                    annotation enters only at featureCounts.
# Results go to results/<mode>/<organism>/ so both modes can coexist; the
# Ensembl-index BAMs live in results/alignments/ and are shared by both modes.
#
# All paths are RELATIVE to the working directory (the repo root).
# ============================================================================

import os
import re

STAR_INDEX_MODE = config.get("star_index_mode", "per_annotation")
# keep_bams: false -> each STAR BAM is a temp() file, deleted by Snakemake as soon
# as every featureCounts job that reads it has finished (saves ~4-5 GB per BAM)
KEEP_BAMS = bool(config.get("keep_bams", False))
if STAR_INDEX_MODE not in ("per_annotation", "ensembl"):
    raise ValueError(f"star_index_mode must be 'per_annotation' or 'ensembl', got {STAR_INDEX_MODE!r}")

REF_DIR = f"reference/{ORGANISM}"
ALN_DIR = f"results/alignments/{ORGANISM}"
RES_ROOT = f"results/{STAR_INDEX_MODE}"
RES_DIR = f"{RES_ROOT}/{ORGANISM}"
ANNOTATION_NAMES = list(ANNOTATIONS.keys())

def index_for(ann):
    """Which STAR index (= which annotation's sjdb) the BAMs for `ann` come from."""
    return ann if STAR_INDEX_MODE == "per_annotation" else ENSEMBL_KEY

INDICES = sorted({index_for(a) for a in ANNOTATION_NAMES})

# --- optional per-organism settings (see header) ---
_GENOMES = globals().get("GENOMES") or {a: {"path": REF, "url": GENOME_URL} for a in ANNOTATION_NAMES}
GENOME_URLS = {g["path"]: g["url"] for g in _GENOMES.values()}
SJDB_AT_MAPPING = globals().get("SJDB_AT_MAPPING", False)
STAR_ENV = globals().get("STAR_ENV", "envs/star.yaml")
FC_STRAND = globals().get("FC_STRAND", 0)
FC_MULTI = globals().get("FC_MULTI", "-O -M")
PREBUILT_INDEX = globals().get("PREBUILT_INDEX") or {}

def index_dir(idx):
    """STAR index used for {idx}: an existing one if given, else built by star_index."""
    return PREBUILT_INDEX.get(idx, f"{REF_DIR}/star_index/{idx}")

def genome_of(idx):
    """Genome FASTA of STAR index {idx}."""
    return _GENOMES[idx]["path"]

# With native per-annotation genomes, a GTF only matches the BAMs of its own
# genome; the shared-Ensembl-index design would silently count ~0 reads.
if globals().get("GENOMES") and STAR_INDEX_MODE != "per_annotation":
    raise ValueError(f"{ORGANISM}: per-annotation genomes (GENOMES) require star_index_mode: per_annotation")

# Constrain wildcards so {ann} cannot greedily match "<ann>.raw" (from the raw
# GTF file name) etc. Without this, DAG building emits spurious
# InputFunctionException/KeyError lines and --strict-dag-evaluation would fail.
wildcard_constraints:
    ann = "|".join(ANNOTATION_NAMES),
    idx = "|".join(ANNOTATION_NAMES),
    sample = "|".join(SAMPLES),

# featureCounts settings (Rsubread argument names in the left column):
#   isPairedEnd=T            -> -p --countReadPairs
#   strandSpecific           -> -s FC_STRAND (pig/mouse 0 = unstranded; dog 2 = reverse)
#   allowMultiOverlap = T    -> -O     (count reads overlapping >1 feature)   } FC_MULTI
#   countMultiMappingReads=T -> -M     (count multi-mapping reads)           } (pig, dog)
#   GTF.featureType="exon"   -> -t exon
#   GTF.attrType="gene_id"   -> -g gene_id
# NOTE: no mapping-quality filter (no -Q), because -Q would drop the multi-
# mapping reads that -M is meant to keep. This permissive counting is where
# annotation choice bites hardest, which is exactly what we want to measure,
# applied IDENTICALLY across all annotations.
FEATURECOUNTS_PARAMS = f"-p --countReadPairs -s {FC_STRAND} {FC_MULTI} -t exon -g gene_id"

# Every final output of one organism (used by `rule all` in each Snakefile).
COUNTS_ONLY = globals().get("COUNTS_ONLY", False)
if COUNTS_ONLY:
    FINAL_OUTPUTS = (
        expand(f"{RES_DIR}/counts/{{ann}}/all_samples.counts.tsv", ann=ANNOTATION_NAMES)
        + expand(f"{REF_DIR}/annotations/{{ann}}.{{t}}.tsv", ann=ANNOTATION_NAMES, t=["genemap", "biotype"])
        + [f"{RES_DIR}/qc/multiqc.html"]
    )
else:
    FINAL_OUTPUTS = (
        expand(f"{RES_DIR}/counts/{{ann}}/all_samples.counts.tsv", ann=ANNOTATION_NAMES)
        + [f"{RES_DIR}/qc/multiqc.html",
           f"{RES_DIR}/fragility/features_table.tsv",
           f"{RES_DIR}/fragility/summary.tsv"]
    )

# NOTE: the default target `rule all` lives at the TOP of each per-organism
# Snakefile (before `include: "common.smk"`), because Snakemake requires the
# first, wildcard-free rule to be in the main workflow file.

# ============================================================================
# REFERENCE: genome FASTA (one, or one per annotation when GENOMES is set)
# ============================================================================
# GENOME url may also be a local file (gzipped or not).
rule download_genome:
    output:
        "{genome_fa}",
    wildcard_constraints:
        genome_fa = "|".join(re.escape(p) for p in GENOME_URLS),
    params:
        url = lambda wc: GENOME_URLS[wc.genome_fa],
    log:
        "logs/download_genome/{genome_fa}.log",
    shell:
        r"""
        ( set -euo pipefail
          mkdir -p $(dirname {output})
          case "{params.url}" in
            http*|ftp*) wget -q -O {output}.gz "{params.url}"; gunzip -f {output}.gz ;;
            *)          gzip -dcf "{params.url}" > {output} ;;
          esac
        ) > {log} 2>&1
        """

# ============================================================================
# READS: paired FASTQ downloaded as .fastq.gz straight from ENA (no SRA
# conversion, so no uncompressed intermediate files), MD5-checked.
# ============================================================================
ENA_FASTQ = globals().get("ENA_FASTQ") or {}
ENA_SEP = globals().get("ENA_SEP", "_")

if ENA_FASTQ:
    rule download_fastq:
        output:
            r1 = f"data/{ORGANISM}/reads/{{sample}}{ENA_SEP}1.fastq.gz",
            r2 = f"data/{ORGANISM}/reads/{{sample}}{ENA_SEP}2.fastq.gz",
        wildcard_constraints:
            sample = "|".join(ENA_FASTQ),
        params:
            u1 = lambda wc: ENA_FASTQ[wc.sample][0][0],
            m1 = lambda wc: ENA_FASTQ[wc.sample][0][1],
            u2 = lambda wc: ENA_FASTQ[wc.sample][1][0],
            m2 = lambda wc: ENA_FASTQ[wc.sample][1][1],
        log:
            f"logs/{ORGANISM}/download/{{sample}}.log",
        threads: 1
        shell:
            r"""
            ( set -euo pipefail
              mkdir -p $(dirname {output.r1})
              wget -q --tries=5 --continue -O {output.r1}.part "https://{params.u1}"
              wget -q --tries=5 --continue -O {output.r2}.part "https://{params.u2}"
              echo "{params.m1}  {output.r1}.part" | md5sum -c -
              echo "{params.m2}  {output.r2}.part" | md5sum -c -
              mv {output.r1}.part {output.r1}
              mv {output.r2}.part {output.r2}
            ) > {log} 2>&1
            """

# ============================================================================
# REFERENCE: raw annotation GTF (one per database)
# ============================================================================
rule download_gtf:
    output:
        f"{REF_DIR}/annotations/{{ann}}.raw.gtf",
    params:
        url = lambda wc: ANNOTATIONS[wc.ann]["url"],
    log:
        f"logs/{ORGANISM}/reference/{{ann}}_download.log",
    shell:
        r"""
        ( set -euo pipefail
          mkdir -p $(dirname {output})
          wget -q -O {output}.gz "{params.url}"
          gunzip -f {output}.gz
        ) > {log} 2>&1
        """

# NCBI assembly report: only needed to translate RefSeq accessions to Ensembl
# chromosome names (see harmonize_gtf).
rule download_assembly_report:
    output:
        f"{REF_DIR}/annotations/refseq.assembly_report.txt",
    params:
        url = lambda wc: next(a["report_url"] for a in ANNOTATIONS.values()
                              if a.get("chrom_fix") == "refseq_map"),
    log:
        f"logs/{ORGANISM}/reference/refseq_report.log",
    shell:
        r"""
        ( set -euo pipefail
          mkdir -p $(dirname {output})
          wget -q -O {output} "{params.url}"
        ) > {log} 2>&1
        """

# ============================================================================
# HARMONISE CHROMOSOME NAMES  (critical; see rationale below)
# ----------------------------------------------------------------------------
# The genome FASTA (and hence every BAM) uses Ensembl names ("1","2",...,"MT").
# RefSeq uses accessions ("NC_..."), UCSC/GENCODE use "chr1". If names
# disagree, STAR ignores those junctions and featureCounts assigns ZERO reads
# to those genes -> a FALSE annotation effect. So every GTF is rewritten into
# Ensembl naming first.
# ============================================================================
def _report_input(wc):
    return (f"{REF_DIR}/annotations/refseq.assembly_report.txt"
            if ANNOTATIONS[wc.ann]["chrom_fix"] == "refseq_map" else [])

# Optional "extra_gtf" (in the repo): extra gene records appended to the
# downloaded GTF, e.g. AMH + HMGA2 for the dog Ensembl 92 annotation.
def _extra_gtf(wc):
    return ANNOTATIONS[wc.ann].get("extra_gtf", [])

rule harmonize_gtf:
    input:
        gtf = f"{REF_DIR}/annotations/{{ann}}.raw.gtf",
        report = _report_input,
        extra = _extra_gtf,
    output:
        f"{REF_DIR}/annotations/{{ann}}.gtf",
    params:
        method = lambda wc: ANNOTATIONS[wc.ann]["chrom_fix"],
    log:
        f"logs/{ORGANISM}/reference/{{ann}}_harmonize.log",
    shell:
        r"""
        ( set -euo pipefail
          case "{params.method}" in
            none)
              cp {input.gtf} {output} ;;
            strip_chr)
              awk 'BEGIN{{OFS=FS="\t"}}
                   /^#/ {{print; next}}
                   {{ sub(/^chr/,"",$1); if($1=="M") $1="MT"; print }}' \
                   {input.gtf} > {output} ;;
            refseq_map)
              # map RefSeq-Accn (col 7) -> Sequence-Name (col 1), then drop any
              # "Chr" prefix so names match Ensembl ("Chr1"->"1"; mouse is already "1")
              awk -F'\t' '!/^#/ && $7!="na" {{n=$1; sub(/^Chr/,"",n); print $7"\t"n}}' \
                   {input.report} > {output}.map
              awk 'BEGIN{{OFS=FS="\t"}}
                   NR==FNR {{m[$1]=$2; next}}
                   /^#/ {{print; next}}
                   {{ if($1 in m) $1=m[$1]; print }}' \
                   {output}.map {input.gtf} > {output}
              rm -f {output}.map ;;
          esac
          for x in {input.extra}; do grep -v '^#' "$x" >> {output}; done
          # drop records with an empty gene_id (""); NCBI RefSeq/UCSC carry tRNA/
          # ncRNA rows that make featureCounts -g gene_id abort (no-op for Ensembl)
          awk -F'\t' '/^#/{{print;next}} $9 ~ /gene_id "[^"]+"/' {output} > {output}.f
          mv {output}.f {output}
        ) > {log} 2>&1
        """

# ============================================================================
# STAR GTF: the harmonised GTF restricted to what STAR can use as sjdb:
# exon records with a non-empty transcript_id, on sequences present in the
# genome FASTA (RefSeq/UCSC also annotate scaffolds absent from Ensembl's FASTA).
# ============================================================================
rule star_gtf:
    input:
        fasta = lambda wc: genome_of(wc.ann),
        gtf   = f"{REF_DIR}/annotations/{{ann}}.gtf",
    output:
        f"{REF_DIR}/annotations/{{ann}}.star.gtf",
    log:
        f"logs/{ORGANISM}/reference/{{ann}}_star_gtf.log",
    shell:
        r"""
        ( set -euo pipefail
          grep '^>' {input.fasta} | sed 's/^>//; s/[[:space:]].*//' > {output}.chroms
          awk -F'\t' 'NR==FNR {{ok[$1]=1; next}}
                      /^#/ {{next}}
                      $3=="exon" && ($1 in ok) && $9 ~ /transcript_id "[^"]+"/' \
              {output}.chroms {input.gtf} > {output}
          echo "exons kept for STAR: $(wc -l < {output})"
          echo "exons on sequences missing from the FASTA (skipped): $(awk -F'\t' 'NR==FNR {{ok[$1]=1; next}} !/^#/ && $3=="exon" && !($1 in ok)' {output}.chroms {input.gtf} | wc -l)"
          rm -f {output}.chroms
        ) > {log} 2>&1
        """

# ============================================================================
# STAR INDEX: one per {idx} (= annotation whose genome/GTF it is built for).
# In "ensembl" mode only the Ensembl index is requested; in "per_annotation"
# mode one index per annotation. STAR version: conda env STAR_ENV.
#   SJDB_AT_MAPPING False: GTF built into the index (--sjdbGTFfile here)
#   SJDB_AT_MAPPING True : genome-only index, the GTF is
#                          inserted at mapping time instead (see star_align)
# ============================================================================
def _star_gtf(idx):
    return f"{REF_DIR}/annotations/{idx}.star.gtf"

def _star_index_input(wc):
    d = {"fasta": genome_of(wc.idx)}
    if not SJDB_AT_MAPPING:
        d["gtf"] = _star_gtf(wc.idx)
    return d

rule star_index:
    input:
        unpack(_star_index_input),
    output:
        directory(f"{REF_DIR}/star_index/{{idx}}"),
    params:
        sjdb = lambda wc: "" if SJDB_AT_MAPPING else
               f"--sjdbGTFfile {_star_gtf(wc.idx)} --sjdbOverhang {SJDB_OVERHANG}",
    log:
        f"logs/{ORGANISM}/star/index_{{idx}}.log",
    threads: 8
    resources:
        mem_mb = 64000,
    conda:
        STAR_ENV
    shell:
        # --limitGenomeGenerateRAM: STAR's default cap (31 GB) is too small for
        # ~2.5-2.7 Gb mammalian genomes; keep it below the job's mem_mb (64 GB).
        r"""
        ( set -euo pipefail
          rm -rf {output} {output}_tmp && mkdir -p {output}
          STAR --version
          STAR --runMode genomeGenerate --runThreadN {threads} \
               --genomeDir {output} --genomeFastaFiles {input.fasta} {params.sjdb} \
               --limitGenomeGenerateRAM 54000000000 \
               --outTmpDir {output}_tmp --outFileNamePrefix {output}/
          rm -rf {output}_tmp
        ) > {log} 2>&1
        """

# ============================================================================
# STAR ALIGN: align each sample against index {idx} (paired-end)
# ============================================================================
def _star_align_input(wc):
    r = get_reads(wc.sample)
    d = {"fq1": r["r1"], "fq2": r["r2"], "idx": index_dir(wc.idx)}
    if SJDB_AT_MAPPING:
        d["gtf"] = _star_gtf(wc.idx)
    return d

rule star_align:
    input:
        unpack(_star_align_input),
    output:
        aln       = (lambda p: p if KEEP_BAMS else temp(p))(
                        f"{ALN_DIR}/{{idx}}/{{sample}}/Aligned.sortedByCoord.out.bam"),
        log_final = f"{ALN_DIR}/{{idx}}/{{sample}}/Log.final.out",
        sj        = f"{ALN_DIR}/{{idx}}/{{sample}}/SJ.out.tab",
    log:
        f"logs/{ORGANISM}/star/{{idx}}/{{sample}}.log",
    params:
        outdir  = f"{ALN_DIR}/{{idx}}/{{sample}}",
        readcmd = lambda wc, input: "--readFilesCommand gunzip -c" if input.fq1.endswith(".gz") else "",
        # junctions from the GTF inserted at mapping time (dog)
        sjdb    = lambda wc: f"--sjdbGTFfile {_star_gtf(wc.idx)} --sjdbOverhang {SJDB_OVERHANG}"
                  if SJDB_AT_MAPPING else "",
        # STAR_EXTRA carries organism-specific flags (e.g. for pig:
        # --outSAMstrandField intronMotif --genomeLoad NoSharedMemory).
        extra   = STAR_EXTRA,
    threads: 8
    resources:
        mem_mb = 64000,
        # each running STAR keeps ~10 GB of sort temp files (_STARtmp) next to its
        # BAM; run_pipeline.sh caps the number of concurrent alignments with
        # --resources star_slots=N to keep the disk quota safe
        star_slots = 1,
    conda:
        STAR_ENV
    shell:
        # --limitBAMsortRAM: by default STAR reserves as much RAM for sorting as
        # the index takes (~27 GB) -> peak ~55 GB; cap it so genome + sorting
        # (~27 + 25 GB) fits in the job's 64 GB. Affects memory only, not output.
        r"""
        ( set -euo pipefail
          rm -rf {params.outdir}/_STARtmp {params.outdir}/_STARgenome
          STAR --version
          STAR --runThreadN {threads} --genomeDir {input.idx} \
               --readFilesIn {input.fq1} {input.fq2} {params.readcmd} \
               --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 25000000000 \
               {params.sjdb} {params.extra} \
               --outFileNamePrefix {params.outdir}/
          rm -rf {params.outdir}/_STARtmp {params.outdir}/_STARgenome
        ) > {log} 2>&1
        """

rule samtools_index:
    input:
        f"{ALN_DIR}/{{idx}}/{{sample}}/Aligned.sortedByCoord.out.bam",
    output:
        f"{ALN_DIR}/{{idx}}/{{sample}}/Aligned.sortedByCoord.out.bam.bai",
    log:
        f"logs/{ORGANISM}/samtools_index/{{idx}}_{{sample}}.log",
    threads: 2
    wrapper:
        "v5.8.3/bio/samtools/index"

# ============================================================================
# FEATURECOUNTS: count reads per gene, ONCE PER ANNOTATION (the comparison!)
# Annotation {ann} is counted on the BAMs from index index_for({ann}).
# ============================================================================
rule featurecounts:
    input:
        samples = lambda wc: f"{ALN_DIR}/{index_for(wc.ann)}/{wc.sample}/Aligned.sortedByCoord.out.bam",
        annotation = f"{REF_DIR}/annotations/{{ann}}.gtf",
    output:
        multiext(
            f"{RES_DIR}/counts/{{ann}}/{{sample}}",
            ".counts.txt",
            ".counts.txt.summary",
        ),
    log:
        f"logs/{ORGANISM}/featurecounts/{STAR_INDEX_MODE}/{{ann}}_{{sample}}.log",
    threads: 4
    params:
        extra = FEATURECOUNTS_PARAMS,
    wrapper:
        "v5.8.3/bio/subread/featurecounts"

# ============================================================================
# MERGE: one count matrix per annotation (rows = genes, cols = samples)
# ============================================================================
rule merge_counts:
    input:
        lambda wc: expand(f"{RES_DIR}/counts/{{ann}}/{{sample}}.counts.txt",
                          ann=wc.ann, sample=SAMPLES),
    output:
        f"{RES_DIR}/counts/{{ann}}/all_samples.counts.tsv",
    log:
        f"logs/{ORGANISM}/merge_counts_{STAR_INDEX_MODE}_{{ann}}.log",
    run:
        import pandas as pd
        merged = None
        for f in input:
            name = os.path.basename(f).replace(".counts.txt", "")
            df = pd.read_csv(f, sep="\t", comment="#")
            df = df[["Geneid", df.columns[-1]]].rename(columns={df.columns[-1]: name})
            merged = df if merged is None else merged.merge(df, on="Geneid", how="outer")
        merged.to_csv(output[0], sep="\t", index=False)

# ============================================================================
# MULTIQC: aggregate STAR logs + featureCounts summaries for ALL annotations
# ============================================================================
rule multiqc:
    input:
        expand(f"{ALN_DIR}/{{idx}}/{{sample}}/Log.final.out", idx=INDICES, sample=SAMPLES),
        expand(f"{RES_DIR}/counts/{{ann}}/{{sample}}.counts.txt.summary",
               ann=ANNOTATION_NAMES, sample=SAMPLES),
    output:
        report = f"{RES_DIR}/qc/multiqc.html",
        data   = directory(f"{RES_DIR}/qc/multiqc_data"),
    log:
        f"logs/{ORGANISM}/multiqc_{STAR_INDEX_MODE}.log",
    wrapper:
        "v9.4.2/bio/multiqc"

# ============================================================================
# GENE TABLES from the GTFs. Independent of the STAR index mode.
#   <ann>.genemap.tsv     gene_id <TAB> gene_name   (cross-annotation matching)
#   ensembl_genes.tsv     gene_id biotype chr start end strand
#   ensembl_ntx.tsv       gene_id n_transcripts
#   <ann>_exons.tsv       gene_id chr start end      (ensembl + refseq)
# ============================================================================
rule genemap:
    input:
        f"{REF_DIR}/annotations/{{ann}}.gtf",
    output:
        f"{REF_DIR}/annotations/{{ann}}.genemap.tsv",
    shell:
        r"""
        awk -F'\t' '/^#/{{next}}{{gid="";gn="";
          if(match($9,/gene_id "[^"]+"/))gid=substr($9,RSTART+9,RLENGTH-10);
          if(match($9,/gene_name "[^"]+"/))gn=substr($9,RSTART+11,RLENGTH-12);
          if(gn=="")gn=gid; if(gid!="" && !(gid in s)){{s[gid]=1;print gid"\t"gn}}}}' \
          {input} > {output}
        """

# gene_id <TAB> biotype. From gene_biotype when the GTF has it (Ensembl, NCBI
# RefSeq); otherwise (UCSC ncbiRefSeq, no gene records) "protein_coding" if the
# gene has an NM_/XM_ transcript, else "non_coding".
rule biotype:
    input:
        f"{REF_DIR}/annotations/{{ann}}.gtf",
    output:
        f"{REF_DIR}/annotations/{{ann}}.biotype.tsv",
    shell:
        r"""
        awk -F'\t' '!/^#/{{gid="";
          if(match($9,/gene_id "[^"]+"/))gid=substr($9,RSTART+9,RLENGTH-10); if(gid=="")next;
          seen[gid]=1;
          if(match($9,/gene_biotype "[^"]+"/))bt[gid]=substr($9,RSTART+14,RLENGTH-15);
          if(match($9,/transcript_id "[^"]+"/)){{t=substr($9,RSTART+15,RLENGTH-16); if(t~/^[NX]M_/)pc[gid]=1}}}}
          END{{for(g in seen) print g"\t"((g in bt)?bt[g]:((g in pc)?"protein_coding":"non_coding"))}}' \
          {input} > {output}
        """

rule exons_table:
    input:
        f"{REF_DIR}/annotations/{{ann}}.gtf",
    output:
        f"{REF_DIR}/annotations/{{ann}}_exons.tsv",
    shell:
        r"""
        awk -F'\t' '$3=="exon"{{if(match($9,/gene_id "[^"]+"/))print substr($9,RSTART+9,RLENGTH-10)"\t"$1"\t"$4"\t"$5}}' \
          {input} > {output}
        """

rule ensembl_gene_features:
    input:
        f"{REF_DIR}/annotations/{ENSEMBL_KEY}.gtf",
    output:
        genes = f"{REF_DIR}/annotations/ensembl_genes.tsv",
        ntx   = f"{REF_DIR}/annotations/ensembl_ntx.tsv",
    shell:
        r"""
        awk -F'\t' '$3=="gene"{{gid=bt="";
          if(match($9,/gene_id "[^"]+"/))gid=substr($9,RSTART+9,RLENGTH-10);
          if(match($9,/gene_biotype "[^"]+"/))bt=substr($9,RSTART+14,RLENGTH-15);
          print gid"\t"bt"\t"$1"\t"$4"\t"$5"\t"$7}}' {input} > {output.genes}
        awk -F'\t' '$3=="transcript"{{if(match($9,/gene_id "[^"]+"/))print substr($9,RSTART+9,RLENGTH-10)}}' {input} \
          | sort | uniq -c | awk '{{print $2"\t"$1}}' > {output.ntx}
        """

# ============================================================================
# DIFFERENTIAL EXPRESSION + FRAGILITY (R). The scripts read RESULTS_DIR /
# REFERENCE_DIR from the environment, so they write into results/<mode>/.
# R is run as ${CONDA_PREFIX}/bin/Rscript (the conda env's own R, not whatever
# "Rscript" comes first in PATH: the login environment, ~/.bashrc and modules
# put a system R 4.4.1 there), with R_HOME/R_LIBS*/LD_LIBRARY_PATH from the
# submitting shell removed, and --vanilla to ignore ~/.Rprofile / .Renviron
# (a personal .libPaths() there would put another R library first). The log
# starts with the R actually used.
R_RUN = (
    'echo "R: ${{CONDA_PREFIX}}/bin/Rscript" > {log}; '
    'env -u R_HOME -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE -u LD_LIBRARY_PATH '
    '"${{CONDA_PREFIX}}/bin/Rscript" --vanilla -e "print(R.home()); print(.libPaths())" >> {log} 2>&1; '
    'env -u R_HOME -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE -u LD_LIBRARY_PATH '
    'RESULTS_DIR={params.root} REFERENCE_DIR=reference '
    '"${{CONDA_PREFIX}}/bin/Rscript" --vanilla'
)
# ============================================================================
if not COUNTS_ONLY:
    rule targets:
        input:
            TARGETS_FILE,
        output:
            f"{RES_DIR}/targets.tsv",
        shell:
            "cp {input} {output}"

    rule de:
        input:
            counts   = expand(f"{RES_DIR}/counts/{{ann}}/all_samples.counts.tsv", ann=ANNOTATION_NAMES),
            genemaps = expand(f"{REF_DIR}/annotations/{{ann}}.genemap.tsv", ann=ANNOTATION_NAMES),
            biotypes = expand(f"{REF_DIR}/annotations/{{ann}}.biotype.tsv", ann=ANNOTATION_NAMES),
            targets  = f"{RES_DIR}/targets.tsv",
        output:
            expand(f"{RES_DIR}/{DE['subdir']}/{{ann}}{DE['suffix']}", ann=ANNOTATION_NAMES),
        params:
            script = os.path.join(workflow.basedir, DE["script"]),
            args   = DE["args"],
            org    = ORGANISM,
            root   = RES_ROOT,
        log:
            f"logs/{ORGANISM}/de_{STAR_INDEX_MODE}.log",
        conda:
            "envs/r-de.yaml"
        shell:
            R_RUN + " {params.script} {params.org} {params.args} >> {log} 2>&1"

    rule fragility:
        input:
            de      = expand(f"{RES_DIR}/{DE['subdir']}/{{ann}}{DE['suffix']}", ann=[ENSEMBL_KEY, "refseq"]),
            genes   = f"{REF_DIR}/annotations/ensembl_genes.tsv",
            ntx     = f"{REF_DIR}/annotations/ensembl_ntx.tsv",
            exons   = f"{REF_DIR}/annotations/{ENSEMBL_KEY}_exons.tsv",
            rexons  = f"{REF_DIR}/annotations/refseq_exons.tsv",
            rmap    = f"{REF_DIR}/annotations/refseq.genemap.tsv",
        output:
            table   = f"{RES_DIR}/fragility/features_table.tsv",
            summary = f"{RES_DIR}/fragility/summary.tsv",
        params:
            script = os.path.join(workflow.basedir, "fragility_features.R"),
            lfc    = FRAGILITY_LFC,
            org    = ORGANISM,
            root   = RES_ROOT,
        log:
            f"logs/{ORGANISM}/fragility_{STAR_INDEX_MODE}.log",
        conda:
            "envs/r-de.yaml"
        shell:
            R_RUN + " {params.script} {params.org} {params.lfc} >> {log} 2>&1"
