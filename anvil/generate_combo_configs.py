#!/usr/bin/env python3
import argparse
import csv
import math
from itertools import product
from pathlib import Path
from typing import Any, Dict, List

import yaml


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _build_combos(parameter_sweeps: Dict[str, List[Any]]) -> List[Dict[str, Any]]:
    if not parameter_sweeps:
        return [{}]
    keys = list(parameter_sweeps.keys())
    values = [parameter_sweeps[k] for k in keys]
    combos = []
    for combo in product(*values):
        combos.append(dict(zip(keys, combo)))
    return combos


def _chunk(items: List[Dict[str, Any]], n_chunks: int) -> List[List[Dict[str, Any]]]:
    n_chunks = max(1, int(n_chunks))
    out: List[List[Dict[str, Any]]] = [[] for _ in range(n_chunks)]
    for i, item in enumerate(items):
        out[i % n_chunks].append(item)
    return [x for x in out if x]


def main() -> None:
    p = argparse.ArgumentParser(description="Generate shard configs for Anvil SLURM array jobs.")
    p.add_argument("--base-config", default="config_astro_full.yaml")
    p.add_argument("--out-dir", default="anvil/generated")
    p.add_argument("--shards", type=int, default=32, help="Number of shard config files")
    p.add_argument("--manifest", default="anvil/generated/manifest.txt")
    p.add_argument("--metadata-csv", default="anvil/generated/metadata.csv")
    p.add_argument("--results-dir", default="results/anvil")
    args = p.parse_args()

    base_config_path = Path(args.base_config)
    out_dir = Path(args.out_dir)
    manifest_path = Path(args.manifest)
    metadata_path = Path(args.metadata_csv)
    results_dir = Path(args.results_dir)

    _ensure_dir(out_dir)
    _ensure_dir(manifest_path.parent)
    _ensure_dir(metadata_path.parent)
    _ensure_dir(results_dir)

    with base_config_path.open("r") as f:
        base = yaml.safe_load(f)
    if base is None:
        base = {}

    sweeps = base.get("parameter_sweeps", {})
    combos = _build_combos(sweeps)
    shard_chunks = _chunk(combos, int(args.shards))

    manifest_lines: List[str] = []
    metadata_rows: List[Dict[str, Any]] = []

    for shard_id, chunk in enumerate(shard_chunks):
        shard_cfg = dict(base)
        shard_cfg["parameter_sweeps"] = {}
        shard_cfg["explicit_conditions"] = chunk
        shard_cfg["output_filename"] = str(results_dir / f"aggregated_shard_{shard_id:04d}.csv")
        shard_cfg["raw_runs_output"] = str(results_dir / f"raw_shard_{shard_id:04d}.csv")
        shard_cfg["mrn_output_filename"] = str(results_dir / f"mrn_shard_{shard_id:04d}.csv")

        shard_path = out_dir / f"config_shard_{shard_id:04d}.yaml"
        with shard_path.open("w") as f:
            yaml.safe_dump(shard_cfg, f, sort_keys=False)
        manifest_lines.append(str(shard_path))

        for local_idx, condition in enumerate(chunk):
            row = {"shard_id": shard_id, "condition_index": local_idx, **condition}
            metadata_rows.append(row)

    with manifest_path.open("w") as f:
        for line in manifest_lines:
            f.write(f"{line}\n")

    if metadata_rows:
        fieldnames = sorted({k for row in metadata_rows for k in row})
        with metadata_path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(metadata_rows)

    print(f"Base config: {base_config_path}")
    print(f"Total combinations: {len(combos)}")
    print(f"Generated shard configs: {len(shard_chunks)}")
    print(f"Manifest: {manifest_path}")
    print(f"Metadata: {metadata_path}")


if __name__ == "__main__":
    main()

