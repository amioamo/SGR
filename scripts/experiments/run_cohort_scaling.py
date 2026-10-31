#!/usr/bin/env python3
"""Cohort-size sensitivity and mechanism controls (manuscript Fig. 1, Fig. 5; Supp. Table S2, Supp. Fig. S2)."""

from __future__ import annotations

import runpy
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

TARGET = ROOT / "scripts" / "run_cohort_scaling.py"

if __name__ == "__main__":
    sys.argv[0] = str(TARGET)
    if "--shared-reference-panel" not in sys.argv:
        sys.argv.append("--shared-reference-panel")
    runpy.run_path(str(TARGET), run_name="__main__")
