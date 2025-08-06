#!/bin/bash
#
# Author: Dolphin Whisperer
# Created: 2025-02-28
# Description: This script If your job uses --requeue, SLURM may send SIGTERM before requeuing. Trap it to save checkpoints before exit.
# --signal=B:USR1@60 sends SIGUSR1 60 seconds before reaching the time limit.
# The trap handles the signal and runs save_checkpoint so work isn't lost.
# --requeue ensures SLURM re-submits the job.
# On restart, the script checks for the checkpoint and resumes rather than starting over.

#SBATCH --job-name=graceful-checkpoint
#SBATCH --output=graceful-ckpt_%j.out
#SBATCH --error=graceful-ckpt_%j.err
#SBATCH --time=2:00:00
#SBATCH --signal=B:USR1@60        # send SIGUSR1 60s before timeout
#SBATCH --requeue

set -eEuo pipefail
trap 'echo "CATCHED SIGUSR1; saving checkpoint..."; save_checkpoint; exit 1' USR1

ml snakemake

save_checkpoint() {
  echo "Snakemake checkpoint at $(date)" >> checkpoint_${SLURM_JOB_ID}.log
  cp -r .snakemake .snakemake_checkpoint_${SLURM_JOB_ID}
}

# Resume if checkpoint exists
if [[ -f checkpoint_${SLURM_JOB_ID}.dat ]]; then
  echo "Resuming from checkpoint"
  # resume logic...
else
  echo "Starting fresh"
fi

# Example workload
for i in $(seq 1 100); do
  echo "Step $i"
  sleep 30
  echo $i > checkpoint_${SLURM_JOB_ID}.dat
done

echo "Completed successfully"
rm -f checkpoint_${SLURM_JOB_ID}.dat
