#!/usr/bin/env python3
"""Pinned AMFlow comparison on independent children of a dependent cut inventory."""
from pathlib import Path
import sys
from upstream_oracle import main

if __name__ == "__main__":
    sys.exit(main(cases=("dependent_child_two", "dependent_child_three"),
                  driver_path=Path(__file__).with_suffix(".wl")))
