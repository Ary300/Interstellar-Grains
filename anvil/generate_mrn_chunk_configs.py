#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import yaml


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _mrn_weights(min_um: float, max_um: float, bins: int) -> tuple[np.ndarray, np.ndarray]:
    sizes = np.logspace(np.log10(min_um), np.log10(max_um), bins)
    weights = sizes ** (-3.5)
    weights = weights / np.sum(weights)
    return sizes, weights


def _chunk_counts(total_runs: int, chunk_size: int) -> List[int]:
    counts: List[int] = []
    remaining = int(total_runs)
    while remaining > 0:
        n = min(int(chunk_size), remaining)
        counts.append(n)
        remaining -= n
    return counts


def main() -> None:
    p = argparse.ArgumentParser(description="Generate memory-safe MRN shard configs split by temperature, size bin, and ensemble chunk.")
    p.add_argument("--base-config", default="config_astro_mrn_integration.yaml")
    p.add_argument("--out-dir", default="anvil/generated_mrn_r9")
    p.add_argument("--manifest", default="anvil/generated_mrn_r9/manifest.txt")
    p.add_argument("--metadata-csv", default="anvil/generated_mrn_r9/metadata.csv")
    p.add_argument("--results-dir", default="results/anvil_mrn_r9")
    p.add_argument("--cache-dir", default="grain_cache_mrn_r9")
    p.add_argument("--chunk-size", type=int, default=5)
    args = p.parse_args()

    base_path = Path(args.base_config)
    out_dir = Path(args.out_dir)
    manifest_path = Path(args.manifest)
    metadata_path = Path(args.metadata_csv)
    results_dir = Path(args.results_dir)

    _ensure_dir(out_dir)
    _ensure_dir(manifest_path.parent)
    _ensure_dir(metadata_path.parent)
    _ensure_dir(results_dir)

    with base_path.open("r") as f:
        base = yaml.safe_load(f) or {}

    temps = list(base.get("parameter_sweeps", {}).get("surface_temperature_k", []))
    nh_values = list(base.get("parameter_sweeps", {}).get("h_gas_density_cm3", []))
    uv_values = list(base.get("parameter_sweeps", {}).get("uv_flux_factor", []))
    if len(nh_values) != 1 or len(uv_values) != 1:
        raise SystemExit("Expected exactly one n(H) and one UV value in the MRN support config.")

    nh = float(nh_values[0])
    uv = float(uv_values[0])
    ensemble_runs = int(base.get("ensemble_runs", 20))
    chunk_counts = _chunk_counts(ensemble_runs, int(args.chunk_size))

    mrn_min_um = float(base.get("mrn_min_um", 0.005))
    mrn_max_um = float(base.get("mrn_max_um", 0.25))
    mrn_bins = int(base.get("mrn_bins", 12))
    sizes_um, weights = _mrn_weights(mrn_min_um, mrn_max_um, mrn_bins)

    manifest_lines: List[str] = []
    metadata_rows: List[Dict[str, Any]] = []
    shard_id = 0
    seed_base = int(base.get("rng_seed", 6200))

    for temp_index, temp in enumerate(temps):
        for size_index, (size_um, weight) in enumerate(zip(sizes_um, weights)):
            run_offset = 0
            for chunk_index, chunk_runs in enumerate(chunk_counts):
                cfg = dict(base)
                cfg["use_mrn"] = False
                cfg["aggregate_across_sizes"] = False
                cfg["save_raw_runs"] = True
                cfg["ensemble_runs"] = int(chunk_runs)
                cfg["min_ensemble_runs"] = int(chunk_runs)
                cfg["grain_radius_um"] = float(size_um)
                cfg["grain_cache_dir"] = str(args.cache_dir)
                cfg["rng_seed"] = seed_base + temp_index * 10000 + size_index * 100 + run_offset
                cfg["parameter_sweeps"] = {}
                cfg["explicit_conditions"] = [
                    {
                        "surface_temperature_k": float(temp),
                        "h_gas_density_cm3": float(nh),
                        "uv_flux_factor": float(uv),
                        "grain_radius_um": float(size_um),
                        "grain_radius_um_mrn": float(size_um),
                        "_mrn_weight": float(weight),
                        "mrn_size_index": int(size_index),
                        "mrn_chunk_index": int(chunk_index),
                        "mrn_chunk_n_runs": int(chunk_runs),
                    }
                ]
                cfg["output_filename"] = str(results_dir / f"aggregated_chunk_{shard_id:04d}.csv")
                cfg["raw_runs_output"] = str(results_dir / f"raw_chunk_{shard_id:04d}.csv")
                cfg["mrn_output_filename"] = str(results_dir / f"mrn_unused_{shard_id:04d}.csv")

                shard_path = out_dir / f"config_shard_{shard_id:04d}.yaml"
                with shard_path.open("w") as f:
                    yaml.safe_dump(cfg, f, sort_keys=False)
                manifest_lines.append(str(shard_path))

                metadata_rows.append(
                    {
                        "shard_id": shard_id,
                        "surface_temperature_k": float(temp),
                        "h_gas_density_cm3": float(nh),
                        "uv_flux_factor": float(uv),
                        "grain_radius_um_mrn": float(size_um),
                        "_mrn_weight": float(weight),
                        "mrn_size_index": int(size_index),
                        "mrn_chunk_index": int(chunk_index),
                        "mrn_chunk_n_runs": int(chunk_runs),
                        "rng_seed": int(cfg["rng_seed"]),
                        "raw_runs_output": cfg["raw_runs_output"],
                        "output_filename": cfg["output_filename"],
                    }
                )

                shard_id += 1
                run_offset += int(chunk_runs)

    with manifest_path.open("w") as f:
        for line in manifest_lines:
            f.write(f"{line}\n")

    if metadata_rows:
        fieldnames = sorted({k for row in metadata_rows for k in row})
        with metadata_path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(metadata_rows)

    print(f"Generated {shard_id} MRN shard configs")
    print(f"Manifest: {manifest_path}")
    print(f"Metadata: {metadata_path}")


if __name__ == "__main__":
    main()
