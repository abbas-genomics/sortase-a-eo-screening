#!/usr/bin/env python3
"""
validate_outputs.py  -  Quick output validator for extract/optimize/combine pipeline
===============================================================================
Checks for expected files and basic row/count consistency.
"""

import argparse
import csv
import re
import sys
from pathlib import Path


def count_csv_rows(path: Path) -> int:
    if not path.exists():
        return -1
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f)
        rows = list(reader)
    if not rows:
        return 0
    return max(0, len(rows) - 1)


def count_numbered_sdf(path: Path) -> int:
    if not path.exists():
        return 0
    files = [p for p in path.glob("*.sdf") if re.match(r"^\d+\.sdf$", p.name)]
    return len(files)


def ok(msg: str):
    print(f"[OK]   {msg}")


def warn(msg: str):
    print(f"[WARN] {msg}")


def fail(msg: str):
    print(f"[FAIL] {msg}")


def validate_extract(base_output: Path) -> int:
    print("\n== Validate extract outputs ==")
    score = 0

    status_csv = base_output / "download_status.csv"
    props_csv = base_output / "compound_properties.csv"
    log_file = base_output / "pipeline.log"
    numbered_dir = base_output / "numbered_sdf"
    renaming_csv = numbered_dir / "renaming_log.csv"
    smiles_csv = numbered_dir / "smiles.csv"

    for p in [status_csv, props_csv, log_file, renaming_csv, smiles_csv]:
        if p.exists():
            ok(f"Found {p}")
            score += 1
        else:
            fail(f"Missing {p}")

    numbered_count = count_numbered_sdf(numbered_dir)
    if numbered_count > 0:
        ok(f"Numbered SDF count = {numbered_count}")
        score += 1
    else:
        warn("No numbered SDF files found")

    rename_rows = count_csv_rows(renaming_csv)
    smiles_rows = count_csv_rows(smiles_csv)
    if rename_rows >= 0 and smiles_rows >= 0:
        if rename_rows == smiles_rows:
            ok(f"renaming_log rows == smiles rows ({rename_rows})")
            score += 1
        else:
            warn(f"Row mismatch: renaming_log={rename_rows}, smiles={smiles_rows}")

    return score


def validate_optimize(opt_output: Path) -> int:
    print("\n== Validate optimize outputs ==")
    score = 0

    sdf_dir = opt_output / "optimized_sdf"
    mol2_dir = opt_output / "optimized_mol2"
    results_csv = opt_output / "optimization_results.csv"
    log_file = opt_output / "optimization.log"

    for p in [sdf_dir, mol2_dir, results_csv, log_file]:
        if p.exists():
            ok(f"Found {p}")
            score += 1
        else:
            fail(f"Missing {p}")

    sdf_count = len(list(sdf_dir.glob("*.sdf"))) if sdf_dir.exists() else 0
    mol2_count = len(list(mol2_dir.glob("*.mol2"))) if mol2_dir.exists() else 0
    if sdf_count > 0:
        ok(f"Optimized SDF count = {sdf_count}")
        score += 1
    else:
        warn("No optimized SDF files found")

    if mol2_count > 0:
        ok(f"Optimized MOL2 count = {mol2_count}")
        score += 1
    else:
        warn("No optimized MOL2 files found")

    results_rows = count_csv_rows(results_csv)
    if results_rows >= 0:
        ok(f"optimization_results rows = {results_rows}")
        score += 1

    return score


def validate_combine(combined_output: Path) -> int:
    print("\n== Validate combine outputs ==")
    score = 0

    combined_sdf = combined_output / "combined_ligands.sdf"
    combined_mol2 = combined_output / "combined_ligands.mol2"
    manifest_csv = combined_output / "ligand_manifest.csv"
    log_file = combined_output / "combine.log"

    for p in [combined_sdf, manifest_csv, log_file]:
        if p.exists():
            ok(f"Found {p}")
            score += 1
        else:
            fail(f"Missing {p}")

    if combined_mol2.exists():
        ok(f"Found optional {combined_mol2}")
        score += 1
    else:
        warn("combined_ligands.mol2 not found (optional unless --mol2 was used)")

    rows = count_csv_rows(manifest_csv)
    if rows >= 0:
        ok(f"ligand_manifest rows = {rows}")
        score += 1

    return score


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate extraction/optimization/combination outputs")
    parser.add_argument("--base-output", default="output", help="Base output folder (default: output)")
    parser.add_argument("--opt-output", default="output/optimized", help="Optimization output folder")
    parser.add_argument("--combined-output", default="output/combined", help="Combined output folder")
    args = parser.parse_args()

    base_output = Path(args.base_output)
    opt_output = Path(args.opt_output)
    combined_output = Path(args.combined_output)

    total = 0
    total += validate_extract(base_output)
    total += validate_optimize(opt_output)
    total += validate_combine(combined_output)

    print("\n== Validation summary ==")
    print(f"Checks passed score: {total}")
    print("Review WARN/FAIL lines above if any.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
