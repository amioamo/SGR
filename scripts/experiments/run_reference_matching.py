#!/usr/bin/env python3
"""Population-matched reference fraction analysis (manuscript Fig. 4)."""

from __future__ import annotations

import runpy
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

TARGET = ROOT / "scripts" / "run_reference_matching.py"

if __name__ == "__main__":
    sys.argv[0] = str(TARGET)
    if "--shared-reference-panel" not in sys.argv:
        sys.argv.append("--shared-reference-panel")
    if "--target-ancestries" not in sys.argv:
        sys.argv.extend(["--target-ancestries", "EUR"])
    if "--wrong-ancestries" not in sys.argv:
        sys.argv.extend(["--wrong-ancestries", "AFR,EAS,SAS,AMR"])
    if "--run-saleem" not in sys.argv:
        sys.argv.append("--run-saleem")
    runpy.run_path(str(TARGET), run_name="__main__")
