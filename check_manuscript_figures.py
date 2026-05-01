#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from PIL import Image, ImageOps


def main() -> None:
    p = argparse.ArgumentParser(description="Generate grayscale QC versions of manuscript PNG figures.")
    p.add_argument("--input-dir", default="results/plots/manuscript")
    p.add_argument("--output-dir", default="results/plots/manuscript/grayscale_qc")
    args = p.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    generated = []
    for png in sorted(input_dir.rglob("*.png")):
        if output_dir in png.parents:
            continue
        rel = png.relative_to(input_dir)
        target = output_dir / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        img = Image.open(png)
        gray = ImageOps.grayscale(img)
        gray.save(target)
        generated.append(str(target))

    report = output_dir / "README.md"
    report.write_text(
        "# Grayscale QC Outputs\n\n"
        + "\n".join(f"- `{path}`" for path in generated)
        + ("\n" if generated else "\nNo PNG figures found.\n"),
        encoding="utf-8",
    )
    print(f"Wrote {len(generated)} grayscale figures to {output_dir}")
    print(f"Wrote {report}")


if __name__ == "__main__":
    main()
