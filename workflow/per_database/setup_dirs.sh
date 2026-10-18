#!/bin/bash
# ============================================================================
# setup_dirs.sh: jednorazowe przygotowanie katalogu repo na klastrze.
# Tworzy katalogi na logi SLURM (SLURM ich NIE tworzy sam; bez nich job
# kończy się bez śladu) oraz katalogi data/<organizm>/ na pobierane dane.
# Uruchom z katalogu repo:  ./setup_dirs.sh
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"

mkdir -p logs/slurm logs/snakemake data/mouse data/pig data/dog


echo "Stan data/:"
ls -l data/*/ 2>/dev/null || true
