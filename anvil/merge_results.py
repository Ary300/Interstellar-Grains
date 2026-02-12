#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path
from typing import Dict, List


def _read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", newline="") as f:
        return list(csv.DictReader(f))


def _write_csv(path: Path, rows: List[Dict[str, str]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({k for row in rows for k in row.keys()})
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    p = argparse.ArgumentParser(description="Merge shard CSV outputs from Anvil runs.")
    p.add_argument("--input-glob", required=True, help="Example: 'results/anvil/aggregated_shard_*.csv'")
    p.add_argument("--output", required=True, help="Merged CSV output path")
    args = p.parse_args()

    paths = sorted(Path(".").glob(args.input_glob))
    if not paths:
        raise SystemExit(f"No files matched: {args.input_glob}")

    all_rows: List[Dict[str, str]] = []
    for path in paths:
        all_rows.extend(_read_csv(path))

    out_path = Path(args.output)
    _write_csv(out_path, all_rows)
    print(f"Merged {len(paths)} files -> {out_path} ({len(all_rows)} rows)")


if __name__ == "__main__":
    main()

