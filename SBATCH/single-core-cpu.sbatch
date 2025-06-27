#!/bin/bash
#SBATCH --job-name=!NT3R3$T!NG_T1M35        # Job name
#SBATCH --output=interesting-times_%j.out   # STDOUT (%j = JobID)
#SBATCH --error=interesting-times_%j.err    # STDERR
#SBATCH --partition=all                     # Partition (all = default)
#SBATCH --ntasks=1                          # Total tasks (MPI ranks)
#SBATCH --cpus-per-task=1                   # CPU cores per task
#SBATCH --mem=1G                            # Total memory (1 GB)
#SBATCH --time=00:01:00                     # Time limit (1 minute)

echo "May you live in interesting times. Regards, $(hostname)"
