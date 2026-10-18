#!/bin/bash
# ============================================================================
# SLURM "master" job for the RNA-seq annotation-comparison pipeline.
# It is itself tiny (2 CPU, 8 GB): its only job is to run Snakemake, which then
# submits every rule as a separate SLURM job via the --executor slurm profile.
# Submit it FROM THE REPO DIRECTORY with:  sbatch run_pipeline.sh
# (first time only: ./setup_dirs.sh; SLURM needs logs/snakemake/ to exist)
# ============================================================================

# ---- SLURM directives (resources for THIS coordinating job only) ----
#SBATCH --nodes=1                       # one node
#SBATCH --ntasks=1                      # one task (the snakemake process)
#SBATCH --cpus-per-task=2               # 2 CPUs are plenty for the scheduler
#SBATCH --mem=8gb                       # 8 GB RAM for the snakemake process
#SBATCH --time=72:00:00                 # 3 days wall-time for the whole run
#SBATCH --partition=biology_standard    # target partition on Eagle
#SBATCH --account=pl0296-02             # grant/accounting code
#SBATCH --mail-type=FAIL,END            # email on failure and on completion (to the submitting account)
#SBATCH --job-name=rnaseq_annotation    # name shown in `squeue`
#SBATCH --output=logs/snakemake/snakemake_%j.log   # stdout (%j = job id), relative to the submit dir
#SBATCH --error=logs/snakemake/snakemake_%j.err    # stderr

# Stop on first error, treat unset variables as errors, fail on pipe errors.
set -euo pipefail

# ---- PATHS ----
# Use the directory the job was submitted from; fall back to this script's dir.
REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"

# Caches MUST live on the shared filesystem, because sub-jobs run on different
# nodes and all need to see the same cache (/tmp is node-local and not shared).
JOB_ROOT="${REPO_DIR}/.snakemake_cache"
export XDG_CACHE_HOME="${JOB_ROOT}/xdg-cache"        # generic cache dir
export PYTHONPYCACHEPREFIX="${JOB_ROOT}/pycache"     # Python bytecode cache
export CONDA_PKGS_DIRS="${JOB_ROOT}/conda-pkgs"      # conda package cache
export TMPDIR="${JOB_ROOT}/tmp"                      # temp files (shared, not /tmp)

# Create every directory we just pointed at, plus the snakemake log dir.
mkdir -p "${XDG_CACHE_HOME}" "${PYTHONPYCACHEPREFIX}" \
         "${CONDA_PKGS_DIRS}" "${TMPDIR}" "${REPO_DIR}/logs/snakemake"

# ---- ACTIVATE CONDA ----
# The env that has snakemake + snakemake-executor-plugin-slurm. Default "base";
# choose another at submit time:  sbatch --export=ALL,CONDA_ENV=snakemake,... run_pipeline.sh
source ~/.bashrc 2>/dev/null || true
eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV:-base}"

# ---- WHICH WORKFLOW ----
# Default = mouse. Choose another organism at submit time, e.g.:
#   sbatch --export=ALL,SNAKEFILE=Snakefile_pig run_pipeline.sh
# STAR index mode comes from config/pipeline.yaml; override per run with, e.g.:
#   sbatch --export=ALL,SNAKEFILE=Snakefile_dog,STAR_INDEX_MODE=ensembl run_pipeline.sh
SNAKEFILE="${SNAKEFILE:-Snakefile_mouse}"
# At most STAR_SLOTS STAR alignments at once (each needs ~10 GB of temp space
# while sorting, besides its final ~4-5 GB BAM). Override: --export=ALL,STAR_SLOTS=4,...
STAR_SLOTS="${STAR_SLOTS:-6}"
EXTRA_CONFIG=()
[ -n "${STAR_INDEX_MODE:-}" ] && EXTRA_CONFIG=(--config "star_index_mode=${STAR_INDEX_MODE}")

# ---- RUN INFO (printed to the master log for provenance) ----
echo "RNA-seq annotation-comparison pipeline (Eagle PCSS, scheduler mode)"
echo "Master Job ID:  ${SLURM_JOB_ID:-<interactive>}"
echo "Timestamp:      $(date)"
echo "Repo dir:       ${REPO_DIR}"
echo "Snakefile:      ${SNAKEFILE}"
echo "Index mode:     ${STAR_INDEX_MODE:-<from config/pipeline.yaml>}"
echo "STAR at once:   ${STAR_SLOTS}"
echo "Snakemake:      $(snakemake --version)"
echo "Account:        pl0296-02"
echo "Partition:      biology_standard"
echo ""

cd "${REPO_DIR}"

# ---- LAUNCH ----
# -s ${SNAKEFILE}            : the per-organism workflow (NOT the deprecated default Snakefile)
# --profile profiles/slurm   : read all submission/resource settings from profiles/slurm/config.yaml
# --rerun-triggers mtime input: only re-run a job if an input file is newer / changed
snakemake -s "${SNAKEFILE}" --profile profiles/slurm --rerun-triggers mtime input \
          --resources star_slots="${STAR_SLOTS}" "${EXTRA_CONFIG[@]}"

echo ""
echo "Pipeline finished at $(date)"
