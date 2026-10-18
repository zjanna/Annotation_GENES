# Per-database alignment and counting (main analysis)

Snakemake workflow (written for the PCSS Eagle cluster, SLURM) that produces the count
matrices of the paper. For every organism and database the reads are aligned with STAR
to **that database's own genome sequence** and counted with featureCounts against the
same database's annotation, so each database enters at both alignment (splice
junctions) and counting, and genome and annotation always share chromosome names.
The workflow stops at the count matrices; differential expression and everything after
it are run by `../../run_all.sh`.

| Snakefile | Genome / annotation per database | STAR | featureCounts |
|---|---|---|---|
| `Snakefile_mouse` | Ensembl 110 primary assembly; RefSeq GCF_000001635.27 (RS_2024_02); GENCODE vM33 primary assembly; UCSC mm39 | 2.7.11b, GTF in the index, `--sjdbOverhang 149 --outSAMstrandField intronMotif` | `-s 0` |
| `Snakefile_pig` | Ensembl 104 toplevel; RefSeq GCF_000003025.6 (AR 106); UCSC susScr11 | 2.7.11b, GTF in the index, `--sjdbOverhang 99 --outSAMstrandField intronMotif` | `-s 0 -O -M` |
| `Snakefile_dog` | Ensembl 92 toplevel + AMH, HMGA2 models (`../../data/dog_added_gene_models/`); RefSeq GCF_000002285.3 (AR 105); UCSC canFam3 | 2.5.3a: genome-only index, junctions from the GTF at mapping (`--sjdbOverhang 100`) | `-s 2 -O -M` |
| `Snakefile_dog10k` | Ensembl 112 Dog10K_Boxer_Tasha toplevel; RefSeq GCF_000002285.5 (AR 106) | as dog | `-s 2 -O -M` |

All featureCounts runs: `-p --countReadPairs -t exon -g gene_id`, no MAPQ filter.
GTF records with an empty `gene_id` are dropped before counting.

## Reads

Sample names and the corresponding SRA/ENA runs: `../../data/sample_runs.tsv`.

- **Mouse**: downloaded from ENA automatically (MD5-checked), `data/mouse/reads/`.
- **Dog** (both assemblies): the 22 libraries deposited in SRA (PRJNA901164) are
  downloaded from ENA automatically.
- **Pig**: `trim_pig.sh` downloads the raw reads from SRA (PRJNA1147879) and trims them
  with BBDuk (BBMap v36.32, as in the anchor study): adapter removal with BBDuk's standard
  settings (`ktrim=r k=23 mink=11 hdist=1 tbo tpe`), quality trimming at Phred 5
  (`qtrim=w trimq=5`) and a minimum read length of 50 (`minlength=50`), into
  `data/pig/fastq/RNA_<i>.{1,2}.trimmed.fastq.gz`, the files the workflow reads.
  `RNA_<i>` is SRA library S3797Nr<i>.

Every input is downloaded from a public repository.

## Running

Run everything from this directory.

```bash
./setup_dirs.sh                                   # log and data directories
BBDUK=/path/to/bbmap/bbduk.sh bash trim_pig.sh   # pig reads: SRA download + trimming
snakemake -s Snakefile_pig --conda-create-envs-only
sbatch --export=ALL,SNAKEFILE=Snakefile_mouse run_pipeline.sh
sbatch --export=ALL,SNAKEFILE=Snakefile_pig   run_pipeline.sh
sbatch --export=ALL,SNAKEFILE=Snakefile_dog   run_pipeline.sh
sbatch --export=ALL,SNAKEFILE=Snakefile_dog10k,STAR_SLOTS=4 run_pipeline.sh
```

`run_pipeline.sh` and `profiles/slurm/config.yaml` contain the SLURM account and
partition of our cluster (`pl0296-02`, `biology_standard`); change them for another
cluster. Requirements: Snakemake 9 with `snakemake-executor-plugin-slurm`, conda (STAR
versions from `envs/`; featureCounts, samtools and MultiQC from snakemake-wrappers
v5.8.3 / v9.4.2). A STAR index needs about 64 GB of RAM, and a full run of one organism
several hundred GB of disk.

Output: `results/per_annotation/<organism>/counts/<database>/all_samples.counts.tsv`.
Copy the matrices to the layout read by `run_all.sh`:

```bash
for o in mouse pig dog; do mkdir -p ../../${o}_local/results/$o && cp -r results/per_annotation/$o/counts ../../${o}_local/results/$o/; done
mkdir -p ../../dog10k_local/results/dog && cp -r results/per_annotation/dog10k/counts ../../dog10k_local/results/dog/
```

The workflow was written by Alicja Beręsewicz.
