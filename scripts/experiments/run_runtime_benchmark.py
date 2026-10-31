#!/usr/bin/env python3
"""Random |M'|=1000 SNP anchor under the Saleem benchmark regime (supplementary analysis; not a main-text figure)."""

from __future__ import annotations

import runpy
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

TARGET = ROOT / "scripts" / "run_runtime_benchmark.py"

if __name__ == "__main__":
    sys.argv[0] = str(TARGET)
    runpy.run_path(str(TARGET), run_name="__main__")
