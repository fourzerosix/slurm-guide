#!/bin/bash
#
# Author: Dolphin Whisperer
# Created: 2025-02-28
# Description: This script wraps a Nextflow pipeline in SBATCH, enabling graceful checkpoint handling when the job is signaled.
# Runs Nextflow under srun so it resides in a SLURM job step and can receive signals correctly.
# Captures SIGUSR1, sends TERM to the Nextflow process to allow graceful shutdown and checkpoint cleanup.
# With --requeue, the script will be resubmitted with the same job ID if terminated—enabling Nextflow’s -resume capability from its last checkpoint.
#
#SBATCH --job-name=nextflowing
#SBATCH --output=nextflow_wrap_%j.out
#SBATCH --error=nextflow_wrap_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --requeue
#SBATCH --signal=B:USR1@300

set -Eeuo pipefail

cancel_nextflow() {
  echo "[$(date)] Caught SIGUSR1—terminating Nextflow gracefully"
  srun --jobid="${SLURM_JOB_ID}.0" kill -TERM "$NEXTFLOW_PID"
}
trap cancel_nextflow USR1 ERR

ml nextflow

# Launch Nextflow as a background job using srun
srun --exclusive nextflow run main.nf -profile slurm &
NEXTFLOW_PID=$!

wait "$NEXTFLOW_PID"
exit_status=$?

if [[ $exit_status -ne 0 ]]; then
  echo "[$(date)] Nextflow exited with $exit_status—requeueing..."
  exit 1
fi

echo "[$(date)] Pipeline completed successfully. If you are quiet enough, you will hear the flow of the universe."
exit 0
