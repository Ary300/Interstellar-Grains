#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def _pick_python() -> str:
    candidates = [
        Path(sys.executable),
        Path("/usr/local/bin/python3"),
        Path("/Library/Frameworks/Python.framework/Versions/3.10/bin/python3"),
        Path("/opt/homebrew/bin/python3"),
        Path("/usr/bin/python3"),
    ]
    for candidate in candidates:
        if not candidate.exists():
            continue
        probe = subprocess.run(
            [
                str(candidate),
                "-c",
                "import importlib.util; mods=['matplotlib','ase','PIL','scipy']; "
                "raise SystemExit(0 if all(importlib.util.find_spec(m) for m in mods) else 1)",
            ],
            cwd=ROOT,
        )
        if probe.returncode == 0:
            return str(candidate)
    raise RuntimeError("No Python interpreter with matplotlib/ase/PIL/scipy found for manuscript figure build")


def main() -> None:
    python = _pick_python()
    subprocess.run(
        [python, "build_mnras_figures.py", "--out-dir", "results/plots/manuscript"],
        cwd=ROOT,
        check=True,
    )


if __name__ == "__main__":
    main()
