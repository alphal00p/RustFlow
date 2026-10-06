#!/usr/bin/env python3
"""Run physical cut benchmarks with the unchanged pinned AMFlow computational code."""

from pathlib import Path
import sys

from upstream_oracle import main


if __name__ == "__main__":
    sys.exit(main(cases=("mixed_cut_bubble", "weighted_three_body"),
                  driver_path=Path(__file__).with_suffix(".wl")))
