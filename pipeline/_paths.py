#!/usr/bin/env python3
"""
==============================================================================
_paths.py - Single source of truth for filesystem locations

Every Python stage imports from here so that the shell scripts (00_config.sh)
and the Python scripts agree on where things live. Before this module existed
the Python stages hardcoded ~/mva and silently ignored MVA_BASE, so overriding
the data root only half worked.

    MVA_BASE=/scratch/mva bash pipeline/run_all.sh

Python adds a script's own directory to sys.path, so `from _paths import ...`
resolves when a stage is run as `python pipeline/03_inheritance_models.py`.
==============================================================================
"""
import os

# Data root. Mirrors BASE in 00_config.sh.
BASE = os.environ.get("MVA_BASE") or os.path.expanduser("~/mva")

RAW = os.path.join(BASE, "data", "raw")
REF = os.path.join(BASE, "data", "ref")
ANNOT = os.path.join(BASE, "data", "annot")
WORK = os.path.join(BASE, "work")
RESULTS = os.path.join(BASE, "results")
LOGS = os.path.join(BASE, "logs")

# Project root: the directory containing this file's parent.
PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def ensure_dirs():
    """Create the working directories if they do not exist yet."""
    for d in (RAW, REF, ANNOT, WORK, RESULTS, LOGS):
        os.makedirs(d, exist_ok=True)


if __name__ == "__main__":
    print("BASE     =", BASE)
    print("RAW      =", RAW)
    print("WORK     =", WORK)
    print("RESULTS  =", RESULTS)
    print("ANNOT    =", ANNOT)
    print("PROJECT  =", PROJECT)
