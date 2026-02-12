#!/bin/bash
set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 <manifest-path> [project-dir] [conda-env]"
  echo "Example: $0 anvil/generated/manifest.txt \$PWD h2-kmc"
  exit 1
fi

MANIFEST=$1
PROJECT_DIR=${2:-$PWD}
CONDA_ENV=${3:-h2-kmc}

N=$(wc -l < "$MANIFEST")
if [ "$N" -le 0 ]; then
  echo "Manifest has no entries: $MANIFEST"
  exit 1
fi

ARRAY_END=$((N - 1))
echo "Submitting array with ${N} tasks (0-${ARRAY_END})"

sbatch \
  --array=0-"$ARRAY_END" \
  --export=ALL,PROJECT_DIR="$PROJECT_DIR",MANIFEST="$MANIFEST",CONDA_ENV="$CONDA_ENV" \
  anvil/run_array.sbatch

