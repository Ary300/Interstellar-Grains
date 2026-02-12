# Anvil CPU Workflow (KMC Astro Campaign)

This folder contains scripts to run `run_sweep.py` as a Slurm array on Anvil.

## 1) Prepare configs locally

Generate shard configs from a base campaign config (`config_astro_full.yaml` by default):

```bash
python anvil/generate_combo_configs.py \
  --base-config config_astro_full.yaml \
  --shards 64 \
  --out-dir anvil/generated \
  --manifest anvil/generated/manifest.txt \
  --metadata-csv anvil/generated/metadata.csv \
  --results-dir results/anvil
```

This creates:
- `anvil/generated/config_shard_XXXX.yaml`
- `anvil/generated/manifest.txt`
- `anvil/generated/metadata.csv`

## 2) Upload to Anvil

From your local machine:

```bash
rsync -av --exclude '__pycache__' --exclude '.pytest_cache' \
  /path/to/Interstellar-Grains-main/ \
  <username>@<anvil-login-host>:/path/to/project/
```

## 3) Set up environment on Anvil

```bash
cd /path/to/project
conda env create -f environment.yml -n h2-kmc || conda env update -f environment.yml -n h2-kmc
conda activate h2-kmc
pytest -q
```

If your Anvil environment requires `module load` before conda, load those modules first.

## 4) Submit Slurm array

Edit `anvil/run_array.sbatch`:
- `#SBATCH -A YOUR_ALLOCATION`
- partition/time/memory to match your allocation policy.

Then submit:

```bash
bash anvil/submit_array.sh anvil/generated/manifest.txt "$PWD" h2-kmc
```

## 5) Monitor and merge

```bash
squeue -u $USER
sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS
python anvil/merge_results.py --input-glob 'results/anvil/aggregated_shard_*.csv' --output results/astro_full_merged.csv
```

If `use_mrn: true`, merge MRN outputs too:

```bash
python anvil/merge_results.py --input-glob 'results/anvil/mrn_shard_*.csv' --output results/astro_full_mrn_merged.csv
```

## Recommended resources

- Use **CPU** nodes for this code path.
- Start with:
  - `--cpus-per-task=1`
  - `--mem=8G`
  - `--time=12:00:00`
- Increase `--time` first if tasks time out.

