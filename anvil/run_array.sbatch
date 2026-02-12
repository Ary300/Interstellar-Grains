#!/bin/bash
#SBATCH -J h2kmc
#SBATCH -A YOUR_ALLOCATION
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH -o logs/h2kmc_%A_%a.out
#SBATCH -e logs/h2kmc_%A_%a.err
#SBATCH --array=0-0

set -euo pipefail

PROJECT_DIR=${PROJECT_DIR:-$PWD}
MANIFEST=${MANIFEST:-$PROJECT_DIR/anvil/generated/manifest.txt}
CONDA_ENV=${CONDA_ENV:-h2-kmc}

cd "$PROJECT_DIR"
mkdir -p logs

# Use your site-specific conda/module setup.
if [ -f "$HOME/.bashrc" ]; then
  source "$HOME/.bashrc"
fi
if command -v conda >/dev/null 2>&1; then
  conda activate "$CONDA_ENV"
fi

TASK_ID=${SLURM_ARRAY_TASK_ID}
CONFIG_PATH=$(sed -n "$((TASK_ID + 1))p" "$MANIFEST")

if [ -z "${CONFIG_PATH:-}" ]; then
  echo "No config for array index ${TASK_ID} in ${MANIFEST}"
  exit 1
fi

echo "[$(date)] TASK=${TASK_ID} CONFIG=${CONFIG_PATH}"
python run_sweep.py "$CONFIG_PATH"
echo "[$(date)] TASK=${TASK_ID} DONE"

