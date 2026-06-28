#!/usr/bin/env python3
"""Compute balanced h11 batches for parallel CYAxiverse pipeline.

Outputs JSON for GitHub Actions matrix strategy or a Makefile snippet.

Usage:
    python tools/create_batches.py                       # print table + JSON
    python tools/create_batches.py --target 3000        # bigger batches
    python tools/create_batches.py --h11-min 4 --h11-max 50
    python tools/create_batches.py --output matrix.json # save to file
"""

import argparse
import json
import math

def compute_batches(h11_min, h11_max, geoms_per_h11, target_axions):
    """Greedy packing of h11 values into batches ~target_axions each."""
    batches = []
    current_lo = None
    current_ax = 0

    for h in range(h11_min, h11_max + 1):
        ax = h * geoms_per_h11
        if current_lo is None:
            current_lo = h
            current_ax = ax
        elif current_ax + ax > target_axions * 1.25:
            batches.append((current_lo, h - 1, current_ax))
            current_lo = h
            current_ax = ax
        else:
            current_ax += ax

    if current_lo is not None:
        batches.append((current_lo, h11_max, current_ax))

    return batches


def format_table(batches, geoms_per_h11):
    lines = []
    lines.append(f"{'Batch':>6}  {'h11 range':>10}  {'geoms':>6}  {'axions':>8}")
    lines.append("-" * 36)
    total_ax = 0
    total_geoms = 0
    for i, (lo, hi, ax) in enumerate(batches, 1):
        n_geom = (hi - lo + 1) * geoms_per_h11
        total_ax += ax
        total_geoms += n_geom
        lines.append(f"{i:>6}  {lo:>4}-{hi:<4}  {n_geom:>6}  {ax:>8}")
    lines.append("-" * 36)
    lines.append(f"{'Total':>6}  {'':>10}  {total_geoms:>6}  {total_ax:>8}")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Compute balanced h11 batches")
    parser.add_argument("--h11-min", type=int, default=4, help="Minimum h11 (default: 4)")
    parser.add_argument("--h11-max", type=int, default=100, help="Maximum h11 (default: 100)")
    parser.add_argument("--geoms-per-h11", type=int, default=10, help="Geometries per h11 (default: 10)")
    parser.add_argument("--target", type=int, default=2000, help="Target axions per batch (default: 2000)")
    parser.add_argument("--output", type=str, default=None, help="Save JSON to file")
    parser.add_argument("--format", choices=["table", "json", "make"], default="table",
                        help="Output format")
    args = parser.parse_args()

    batches = compute_batches(args.h11_min, args.h11_max,
                              args.geoms_per_h11, args.target)

    if args.format == "json" or args.output:
        result = []
        for i, (lo, hi, ax) in enumerate(batches, 1):
            result.append({
                "batch": i,
                "h11_min": lo,
                "h11_max": hi,
                "geometries": (hi - lo + 1) * args.geoms_per_h11,
                "axions": ax,
            })
        json_out = json.dumps(result, indent=2)
        if args.output:
            with open(args.output, "w") as f:
                f.write(json_out)
            print(f"Saved {len(result)} batches to {args.output}")
        else:
            print(json_out)

    if args.format == "table":
        print(format_table(batches, args.geoms_per_h11))

    if args.format == "make":
        # Print Makefile-compatible target list
        targets = " ".join(
            f"batch_{lo}_{hi}"
            for lo, hi, _ in batches
        )
        print(f"BATCH_TARGETS = {targets}")
        for lo, hi, _ in batches:
            n_geom = (hi - lo + 1) * args.geoms_per_h11
            print(f"\n# h11={lo}-{hi}: {n_geom} geoms")
            print(f"batch_{lo}_{hi}:")
            print(f"\t@echo 'Running batch h11={lo}-{hi}'")


if __name__ == "__main__":
    main()
