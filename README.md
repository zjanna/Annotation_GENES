# Annotation choice and RNA-seq differential expression across organisms

Code accompanying the manuscript *"Ensembl and RefSeq disagree more on differential expression
in less mature genome annotations: a three-organism re-analysis"* (Zyprych-Walczak,
Beręsewicz, Skotarczak; *Genes*, revised version). The results in the paper were produced with
release `v1.0-revision1` of this repository.
<!-- update the title and add the DOI once published -->

For three organisms of decreasing annotation maturity (mouse, pig, dog) we reproduce a
published differential-expression (DE) analysis, repeat it with each annotation
database's own genome and annotation (Ensembl, RefSeq, UCSC; GENCODE for mouse), match
genes between databases by genomic coordinates (strand-specific exon overlap) and by NCBI
Gene ID, and compare the DE gene lists. Dog is also
analysed on the current reference assembly (Dog10K_Boxer_Tasha) with the latest
annotations.

No data are distributed with this repository; everything is downloaded from public
repositories (below).

## Data

| Organism | Raw RNA-seq | Anchor study | Contrast |
|---|---|---|---|
| Mouse (liver, Efcab4b KO vs WT, n = 3 + 3) | ArrayExpress [E-MTAB-13110](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13110) (ENA ERP148613) | Cheng et al. 2025, *Biochem Biophys Rep* 41:101944 | KO vs WT |
| Pig (muscle, umbilical hernia vs control, n = 15 + 15) | NCBI SRA [PRJNA1147879](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1147879) | Woźniak et al. 2025, *BMC Genomics* 26:1096 | UH vs control |
| Dog (testis, cryptorchidism) | NCBI SRA [PRJNA901164](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA901164) | Stachowiak et al. 2024, *PNAS* 121:e2312724121 | undescended vs control (UD vs C) |

**Mouse sample labels.** In the E-MTAB-13110 sample table the genotype labels are swapped
relative to the FASTQ file names (KO1L–KO3L, WT1L–WT3L). We follow the file names, which
match the published fold-change signs. Sample ids `WT_1..KO_3` in
`mouse_local/results/mouse/targets.tsv` correspond to ENA runs ERR11593203–ERR11593208
in that order; their group (`condition`) is the one implied by the file names.

**Annotations** (all version-pinned and downloaded by `workflow/reference_tables.sh`):

| Organism | Ensembl | RefSeq | UCSC | GENCODE |
|---|---|---|---|---|
| Mouse (GRCm39) | release 110 | RS_2024_02 | mm39 ncbiRefSeq | vM33 |
| Pig (Sscrofa11.1) | release 104 | Annotation Release 106 | susScr11 ncbiRefSeq | – |
| Dog (CanFam3.1) | release 92 + AMH, HMGA2* | Annotation Release 105 | canFam3 ncbiRefSeq | – |
| Dog (Dog10K_Boxer_Tasha) | release 112 | Annotation Release 106 | – | – |

\* Ensembl 92 has no gene record for AMH or HMGA2; gene models for these two genes are
appended to the Ensembl 92 annotation and provided in `data/dog_added_gene_models/`.
UCSC does not version its ncbiRefSeq tracks, so files downloaded later may differ
slightly from those used for the paper (October 2026).

## Repository layout

```
run_all.sh                       laptop part: DE -> fragile genes -> tables -> figures
scripts/                         R scripts (see below)
workflow/
  per_database/                  alignment and counting per database (cluster, Snakemake + SLURM)
  reference_tables.sh            annotation tables + Gene ID cross-references (laptop)
  calibration/                   choice of the mouse counting settings (cluster)
data/dog_added_gene_models/      AMH and HMGA2 gene models added to Ensembl 92 (dog)
data/sample_runs.tsv             sample names used here -> SRA/ENA runs
<key>_local/results/<org>/targets.tsv   sample sheets (key = mouse, pig, dog, dog10k)
```

| Script | Purpose | Paper |
|---|---|---|
| `scripts/de_deseq2.R` | DESeq2 per database (mouse; padj < 0.05) | Table 1 |
| `scripts/de_analysis.R` | limma-voom + duplicateCorrelation per database (pig; model of the anchor study) | Table 1 |
| `scripts/de_edger.R` | edgeR likelihood-ratio test per database (dog; method of the anchor study) | Table 1 |
| `scripts/matching.R` | gene matching between databases by NCBI Gene ID; agreement and share of DEGs without a counterpart | all |
| `scripts/coord_match.R` | gene matching Ensembl–RefSeq by genomic coordinates (three overlap criteria); share of DEGs without a counterpart, Jaccard on loci | Table 3, Figs 2, 4, 5 |
| `scripts/coord_results.R` | reads the coordinate-matching results for the figure scripts | Figs 2, 4, 5 |
| `scripts/fragility_coord.R` | stable vs annotation-sensitive DEGs, each database as the reference, coordinate matching; AUC with bootstrap CI | Tables 5, 6, Fig. 6 |
| `scripts/de_common.R` | one common DE method (DESeq2, edgeR QL) for all organisms; Jaccard between DE methods (scale) | Table 4 |
| `scripts/topn_agreement.R` | top-N agreement and log2FC correlation (independent of DEG-list size) | Results |
| `workflow/refseq_curation.sh` | share of curated (NM_/NR_) RefSeq records, a maturity proxy independent of cross-references | Results |
| `scripts/summary_tables.R` | Gene-ID agreement numbers quoted in the Results; the qPCR-tested pig genes | Results, Fig. 1, Table 7 |
| `scripts/xref_coverage.R` | share of Ensembl genes with an NCBI Gene ID (dog Ensembl 92 vs 112, pig 104 vs 112) | Results |
| `scripts/calibration.R` | comparison with the published results (mouse, pig; dog in `run_all.sh`) | Table 2 |
| `scripts/fig_heatmap.R` … `fig_fragility.R` | Figures 1–6 | Figs 1–6 |
| `scripts/fragility_features.R`, `scripts/fragility_effects.R` | fragile vs stable Ensembl DEGs with Gene-ID matching, as in the first submitted version (superseded by `fragility_coord.R`; also called by the optional fragility rule of the Snakemake workflow) | – |

## How to reproduce

**1. Count matrices (cluster).** `workflow/per_database/` downloads the genomes,
annotations and (for mouse and dog) the reads, aligns each organism separately for every
database with that database's own genome and annotation, and counts the reads; see
`workflow/per_database/README.md` for STAR and featureCounts settings and for the pig
reads. The matrices are then copied to

```
<key>_local/results/<org>/counts/<database>/all_samples.counts.tsv   # Geneid + one column per sample
```

(key/org: `mouse/mouse`, `pig/pig`, `dog/dog`, `dog10k/dog`). The dog analysis uses the 22
libraries deposited in SRA (PRJNA901164).

**2. Reference tables (laptop).**

```bash
bash workflow/reference_tables.sh
```

**3. Everything else (laptop).**

```bash
bash run_all.sh
```

Outputs: `figures/fig1_jaccard.pdf` … `fig6_fragility.pdf` and `tables/` (including
`tables/coord_match/`, `tables/fragility_coord/`, `tables/de_common/`, `tables/topn_agreement/`).

**Calibration.** Before comparing databases, each anchor study's Ensembl DE result is
reproduced from the public raw data with the anchor study's DE method (`run_all.sh`,
`scripts/calibration.R`): pig 234 DEGs at FDR < 0.05 and 59 with abs(log2FC) > 1, the
published gene lists; mouse 85 DEGs vs 69 published (Ensembl 110, unstranded,
featureCounts defaults, chosen as closest to the published results;
`workflow/calibration/run_mouse_calib.sh`); dog, protein-coding genes as in the published
analysis, 7526 DEGs vs 7619 published on the 22 public libraries
(`tables/calibration_dog.txt`). All inputs are public.

## Requirements

- R 4.3 with edgeR 4.0, limma 3.58, DESeq2 1.42, GenomicRanges 1.54, ggplot2 3.5, ggrepel, tidyr, dplyr
- STAR 2.7.11b (mouse, pig) and 2.5.3a (dog), featureCounts v2.0.8 (Subread)
- Snakemake 9 with the SLURM executor plugin and conda (cluster part); curl, gzip, awk (reference tables)

## Licence and citation

Code released under the MIT licence (see `LICENSE`).
<!-- add the article citation and DOI once published -->
