#!/bin/bash
#SBATCH --job-name=P@$$!NG_M3SS@G3S
#SBATCH --output=passing-messages_%j.out
#SBATCH --partition=all
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G               # (Optional)
#SBATCH --time=00:01:00

ml openmpi/5.0.7-rwefebl
ml osu-micro-benchmarks/7.5-23m3nxv

srun --mpi=pmix osu_bw
