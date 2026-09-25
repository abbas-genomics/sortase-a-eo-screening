#!/usr/bin/env python3
"""
parse_affinities.py - Rank docked ligands by best (lowest) binding affinity
=============================================================================
Reads all virtual_screening/logs/ligand_N_log.txt files, extracts the mode 1
(best pose) affinity from each Vina run, joins with compound names from
ligand_manifest.csv, and writes a ranked CSV.
"""

import argparse
import csv
import re
from pathlib import Path


def parse_log(log_path: Path):
    """Extract the best (mode 1) binding affinity in kcal/mol from a Vina log."""
    text = log_path.read_text(errors="ignore")
    # Vina result rows look like:  "   1        -7.234      0.000      0.000"
    matches = re.findall(
        r"^\s*(\d+)\s+(-?\d+\.\d+)\s+\d+\.\d+\s+\d+\.\d+",
        text,
        re.MULTILINE,
    )
    if not matches:
        return None
    # mode 1 (first row) is always the best pose per Vina's own ranking,
    # but take the minimum defensively in case ordering ever looks odd
    best = min(matches, key=lambda m: float(m[1]))
    return float(best[1])


def load_manifest(manifest_path: Path):
    """Map numbered filename (e.g. '68.sdf') -> compound name."""
    name_map = {}
    if not manifest_path.exists():
        print(f"WARNING: manifest not found at {manifest_path} — names will show as 'unknown'")
        return name_map
    with open(manifest_path, newline="") as f:
        for row in csv.DictReader(f):
            name_map[row.get("numbered_file", "")] = row.get("compound_name", "unknown")
    return name_map


def main():
    ap = argparse.ArgumentParser(description="Rank docked ligands by binding affinity")
    ap.add_argument("--logs-dir", type=Path, default=Path("virtual_screening/logs"))
    ap.add_argument("--manifest", type=Path, default=Path("output/combined/ligand_manifest.csv"))
    ap.add_argument("--output", type=Path, default=Path("virtual_screening/results/binding_affinity_ranked.csv"))
    args = ap.parse_args()

    name_map = load_manifest(args.manifest)

    log_files = sorted(
        args.logs_dir.glob("ligand_*_log.txt"),
        key=lambda p: int(re.search(r"ligand_(\d+)_log", p.name).group(1)),
    )

    if not log_files:
        print(f"ERROR: no log files found in {args.logs_dir}")
        return 1

    results = []
    missing = []
    for log_file in log_files:
        idx = int(re.search(r"ligand_(\d+)_log", log_file.name).group(1))
        affinity = parse_log(log_file)
        name = name_map.get(f"{idx}.sdf", name_map.get(f"{idx}.mol2", "unknown"))
        if affinity is None:
            missing.append(idx)
        results.append({
            "ligand_id": idx,
            "compound_name": name,
            "best_affinity_kcal_mol": affinity,
        })

    # Sort: valid affinities first (most negative = best), missing ones last
    results.sort(key=lambda r: (r["best_affinity_kcal_mol"] is None, r["best_affinity_kcal_mol"]))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["rank", "ligand_id", "compound_name", "best_affinity_kcal_mol"])
        writer.writeheader()
        for rank, r in enumerate(results, start=1):
            r["rank"] = rank
            r["rank"] = rank
            writer.writerow({"rank": rank, **r})
    print(f"Wrote {len(results)} ligands to {args.output}")
    print("")
    print("Top 10 by predicted binding affinity:")
    print(f"{'Rank':<6}{'Ligand':<10}{'Compound':<30}{'Affinity (kcal/mol)'}")
    for r in results[:10]:
        aff = r["best_affinity_kcal_mol"]
        aff_str = f"{aff:.2f}" if aff is not None else "N/A"
        print(f"{r['rank']:<6}{r['ligand_id']:<10}{r['compound_name']:<30}{aff_str}")

    if missing:
        print("")
        print(f"WARNING: {len(missing)} ligand(s) had no parseable result: {missing}")
        print("Check these logs individually — likely a failed docking run.")

    return 0


if __name__ == "__main__":
    exit(main())
