#!/bin/bash
#SBATCH --job-name=P@$$!NG_M3SS@G3S
#SBATCH --output=passing-messages_%j.out
#SBATCH --partition=all
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4      # 4 CPU cores for threading
#SBATCH --mem=4G               # (Optional - have encountered issues without a --mem= specification when using MPI) (srun: fatal: SLURM_MEM_PER_CPU, SLURM_MEM_PER_GPU, and SLURM_MEM_PER_NODE are mutually exclusive.)
#SBATCH --time=00:01:00

ml openmpi/5.0.7-rwefebl
ml osu-micro-benchmarks/7.5-23m3nxv

srun --mpi=pmix osu_bw
