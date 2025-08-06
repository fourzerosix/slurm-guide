#!/bin/bash
#
# Author: Dolphin Whisperer
# Email: jeremy.bell@nih.gov
# Created: 2025-02-28
# Description: This script uses a wrapper SBATCH script to invoke Snakemake and capture signals (you will need a valid Snakefile for this script to run)
# Captures SIGUSR1 (sent 5 min before timeout) and re-runs incomplete tasks via Snakemake’s --rerun-incomplete mechanism.
# Uses --requeue, retaining the same job ID and priority if interrupted.
# Snakemake's --cluster mode automatically submits individual rules as SLURM jobs, mapping threads, mem_mb, and runtime_min to SBATCH directives

#SBATCH --job-name=snakemake-on-a-plane
#SBATCH --output=snake_wrap_%j.out
#SBATCH --error=snake_wrap_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --requeue
#SBATCH --signal=B:USR1@300

set -Eeuo pipefail
trap 'echo "[$(date)] Signal USR1 received—checkpointing and requeueing..."; snakemake --rerun-incomplete; exit 1' USR1 ERR

ml snakemake

# Run Snakemake with SLURM integration
snakemake \
  --snakefile Snakefile \
  --cores 4 \
  --jobs 100 \
#  --cluster "sbatch --parsable --partition=all --cpus-per-task={threads} --mem={resources.mem_mb//1024}G --time={resources.runtime_min//60}:00" \
#  --latency-wait 60 \
#  --rerun-incomplete
