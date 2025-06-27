#!/bin/bash
##SBATCH --job-name=container_job
#SBATCH --output=ctr_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

set -Eeuo pipefail
trap 'echo "Caught SIGINT in script"; exit 1' INT

srun apptainer exec /data/user-utils/containers/sl.sif sl -e &
PID=$!

echo "Launched work: PID=$PID; PGID=-$PID"

# Wait a bit, then send simulated Ctrl+C
sleep 4
kill -INT -"$PID"
wait "$PID"
