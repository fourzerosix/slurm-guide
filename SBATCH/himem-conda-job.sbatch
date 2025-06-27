#!/bin/bash
#SBATCH --job-name=matplotlib_himem
#SBATCH --output=matplotlib_%j.out
#SBATCH --error=matplotlib_%j.err
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=24:00:00

# Enable robust bash behavior
set -Eeuo pipefail
trap 'echo "[$(date)] Error in ${BASH_SOURCE[0]} at line ${BASH_LINENO[0]}: \"$BASH_COMMAND\" failed." >&2; exit 1' ERR

echo "[$(date)] Starting job on $(hostname). SLURM_JOB_ID=$SLURM_JOB_ID"

# Load Miniconda
module load miniconda3        # adjust to exact module name if needed :contentReference[oaicite:1]{index=1}

# Create env only if not already existing
ENV_NAME=matplotlib_3.10.0
YAML_FILE=/data/user-utils/conda-envs/matplotlib_3.10.0.yaml

# Set up env path in a project or scratch directory
ENV_DIR=$HOME/conda_envs/$ENV_NAME

if [[ ! -d "$ENV_DIR" ]]; then
    echo "[$(date)] Creating Conda env from $YAML_FILE"
    conda env create --prefix "$ENV_DIR" --file "$YAML_FILE"
else
    echo "[$(date)] Using existing Conda env at $ENV_DIR"
fi

# Activate env
eval "$(conda shell.bash hook)"
conda activate "$ENV_DIR"

echo "[$(date)] Activated Conda env: $(conda info --envs | grep '*' )"

# Run high-memory Python task
echo "[$(date)] Running high-memory Python script"
python - <<'EOF'
import numpy as np, matplotlib.pyplot as plt

# Example high-memory operation: large matrix
arr = np.random.random((200000, 1000))  # ~160 GB
print("Matrix shape:", arr.shape)
# Simple plot to use matplotlib
plt.hist(arr.ravel()[:1000000], bins=50)
plt.title("Sample histogram")
plt.savefig("histogram.png")
print("Histogram saved")
EOF

echo "[$(date)] Job completed successfully."
