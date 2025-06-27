#!/bin/bash
#SBATCH --job-name=i_like_pictures
#SBATCH --output=pretty-pictures_%j.out
#SBATCH --partition=gpu      # GPU partition (if “gpu” is the correct name)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:1         # Request 1 GPU
#SBATCH --mem=16G
#SBATCH --time=01:00:00

ml cuda       # load CUDA or other GPU driver modules as needed

srun nvidia-smi
