#!/bin/bash
#SBATCH --job-name=infinite_loop_bash
#SBATCH --output=infinite_loop_bash_%j.out
#SBATCH --error=infinite_loop_bash_%j.err
#SBATCH --partition=all
#SBATCH --ntasks=1
#SBATCH --time=7-00:00:00        # 7 days (adjust as needed)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jeremy.bell@nih.gov

EMAIL="jeremy.bell@nih.gov"

echo "Starting infinite loop job at $(date) on $HOSTNAME"
iteration=0

while true; do
    iteration=$((iteration + 1))
    echo "[$(date)] Iteration $iteration: Launching srun bash..."

    srun bash || {
        LOG="iteration_${iteration}_fail_${SLURM_JOB_ID}.log"
        echo "[$(date)] ERROR: srun bash failed on iteration $iteration" > $LOG
        echo -e "\nSLURM ENVIRONMENT:" >> $LOG
        env | grep SLURM_ >> $LOG
        echo -e "\nSLURM JOB STATUS:" >> $LOG
        scontrol show job $SLURM_JOB_ID >> $LOG
        mail -s "Iteration $iteration failed in job $SLURM_JOB_ID" "$EMAIL" < $LOG
    }

    echo "[$(date)] Sleeping for 5 minutes before next run..."
    sleep 300  # 5 minutes
done
