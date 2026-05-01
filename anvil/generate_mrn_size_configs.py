#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd
import yaml


def mrn_weights(min_um: float, max_um: float, bins: int):
    sizes = [
        10 ** (math.log10(min_um) + (math.log10(max_um) - math.log10(min_um)) * i / (bins - 1))
        for i in range(bins)
    ]
    weights = [a ** (-3.5) for a in sizes]
    total = sum(weights)
    weights = [w / total for w in weights]
    return sizes, weights


def main() -> None:
    p = argparse.ArgumentParser(description="Generate one-size-per-job MRN configs.")
    p.add_argument("--base-config", default="config_astro_mrn_integration.yaml")
    p.add_argument("--out-dir", default="anvil/generated_mrn_r8")
    p.add_argument("--results-dir", default="results/anvil_mrn_r8")
    p.add_argument("--cache-dir", default="grain_cache_mrn_r8")
    args = p.parse_args()

    base_path = Path(args.base_config)
    out_dir = Path(args.out_dir)
    results_dir = Path(args.results_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with base_path.open() as f:
        base = yaml.safe_load(f)

    temps = list(base["parameter_sweeps"]["surface_temperature_k"])
    nh_values = list(base["parameter_sweeps"]["h_gas_density_cm3"])
    uv_values = list(base["parameter_sweeps"]["uv_flux_factor"])
    if len(nh_values) != 1 or len(uv_values) != 1:
        raise SystemExit("This generator assumes a single n_H and UV slice.")

    sizes, weights = mrn_weights(
        float(base["mrn_min_um"]),
        float(base["mrn_max_um"]),
        int(base["mrn_bins"]),
    )

    manifest = []
    metadata_rows = []
    shard_id = 0
    for temp in temps:
        for size_um, weight in zip(sizes, weights):
            cfg = dict(base)
            cfg["use_mrn"] = False
            cfg["aggregate_across_sizes"] = False
            cfg["save_raw_runs"] = True
            cfg["ensemble_runs"] = int(base["ensemble_runs"])
            cfg["min_ensemble_runs"] = int(base["min_ensemble_runs"])
            cfg["grain_radius_um"] = float(size_um)
            # These extra fields flow through run_one_condition into raw rows.
            cfg["grain_radius_um_mrn"] = float(size_um)
            cfg["_mrn_weight"] = float(weight)
            cfg["grain_cache_dir"] = args.cache_dir
            cfg["output_filename"] = str(results_dir / f"aggregated_chunk_{shard_id:04d}.csv")
            cfg["raw_runs_output"] = str(results_dir / f"raw_chunk_{shard_id:04d}.csv")
            cfg["mrn_output_filename"] = str(results_dir / f"mrn_unused_{shard_id:04d}.csv")
            cfg["parameter_sweeps"] = {
                "surface_temperature_k": [float(temp)],
                "h_gas_density_cm3": [float(nh_values[0])],
                "uv_flux_factor": [float(uv_values[0])],
            }

            cfg_path = out_dir / f"config_shard_{shard_id:04d}.yaml"
            with cfg_path.open("w") as f:
                yaml.safe_dump(cfg, f, sort_keys=False)

            manifest.append(str(cfg_path))
            metadata_rows.append(
                {
                    "shard_id": shard_id,
                    "config_path": str(cfg_path),
                    "surface_temperature_k": float(temp),
                    "h_gas_density_cm3": float(nh_values[0]),
                    "uv_flux_factor": float(uv_values[0]),
                    "grain_radius_um_mrn": float(size_um),
                    "_mrn_weight": float(weight),
                    "output_filename": cfg["output_filename"],
                    "raw_runs_output": cfg["raw_runs_output"],
                }
            )
            shard_id += 1

    (out_dir / "manifest.txt").write_text("\n".join(manifest) + "\n")
    pd.DataFrame(metadata_rows).to_csv(out_dir / "metadata.csv", index=False)
    print(f"Wrote {len(manifest)} configs to {out_dir}")


if __name__ == "__main__":
    main()
