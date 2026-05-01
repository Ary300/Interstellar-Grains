#!/bin/bash
set -euo pipefail

REMOTE_ROOT="/anvil/projects/x-bio260046/Interstellar-Grains-main"
REMOTE_HOST="x-adas17@anvil.rcac.purdue.edu"
SSH_CMD='ssh -4 -o IdentityAgent=none -o GSSAPIAuthentication=no -o PreferredAuthentications=publickey -o PubkeyAuthentication=yes -o KbdInteractiveAuthentication=no -o PasswordAuthentication=no -o BatchMode=yes -o ConnectTimeout=8 -i ~/.ssh/anvil_ed25519 -o IdentitiesOnly=yes'

mkdir -p results results/tables

echo "Pulling support-campaign CSVs from Anvil..."

rsync -av -e "$SSH_CMD" \
  --prune-empty-dirs \
  --include='*/' \
  --include='astro_transition_deep.csv' \
  --include='astro_transition_deep_raw.csv' \
  --include='astro_lh_mode_consistency.csv' \
  --include='astro_grain_size_check.csv' \
  --include='astro_sticking_model_check.csv' \
  --include='astro_porosity_check.csv' \
  --include='astro_transport_sensitivity.csv' \
  --include='astro_sensitivity_knobs.csv' \
  --exclude='*' \
  "${REMOTE_HOST}:${REMOTE_ROOT}/results/" \
  results/

rsync -av -e "$SSH_CMD" \
  --prune-empty-dirs \
  --include='*/' \
  --include='astro_transition_deep_distribution_summary.csv' \
  --include='astro_sensitivity_knobs_summary.csv' \
  --include='astro_transport_sensitivity_summary.csv' \
  --exclude='*' \
  "${REMOTE_HOST}:${REMOTE_ROOT}/results/tables/" \
  results/tables/

echo
echo "Pulled files now present locally:"
find results -maxdepth 2 -type f | egrep 'astro_transition_deep(_raw)?\.csv|astro_lh_mode_consistency\.csv|astro_grain_size_check\.csv|astro_sticking_model_check\.csv|astro_porosity_check\.csv|astro_transport_sensitivity\.csv|astro_sensitivity_knobs\.csv|astro_transition_deep_distribution_summary\.csv|astro_sensitivity_knobs_summary\.csv|astro_transport_sensitivity_summary\.csv' | sort || true
