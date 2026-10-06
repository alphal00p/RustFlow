#!/usr/bin/env python3
"""Independent nonfactorized cut oracle with original AMFlow mass placements."""
from pathlib import Path
import sys
from upstream_oracle import main

if __name__ == "__main__":
    sys.exit(main(cases=("partial_three_body", "all_three_body"),
                  driver_path=Path(__file__).with_suffix(".wl")))
