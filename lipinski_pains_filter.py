#!/usr/bin/env python3
"""
lipinski_pains_filter.py - Apply Lipinski Ro5 and PAINS filtering
====================================================================
Reads druglikeness properties (from extract.py's compound_properties.csv)
and SMILES (from ligand_manifest.csv), computes Lipinski Rule of Five
compliance, runs RDKit's PAINS filter catalog, and writes a combined
pass/fail table.

Requires: pip install rdkit
"""

import argparse
import csv
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams


def load_properties(path: Path):
    """Map compound name -> {MW, XLogP, HBD, HBA, RotatableBonds}."""
    props = {}
    if not path.exists():
        print(f"WARNING: properties file not found at {path}")
        return props
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            name = row.get("InputName") or row.get("input_name") or row.get("Name")
            if not name:
                # Fall back to trying every column name variant extract.py might use
                for key in row:
                    if "name" in key.lower():
                        name = row[key]
                        break
            if name:
                props[name.strip().lower()] = row
    return props


def load_manifest(path: Path):
    """Return list of dicts: numbered_file, compound_name, canonical_smiles."""
    rows = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows


def lipinski_check(mw, xlogp, hbd, hba):
    """Return (violations, pass) — Lipinski allows up to 1 violation."""
    violations = 0
    try:
        if mw is not None and float(mw) > 500:
            violations += 1
        if xlogp is not None and float(xlogp) > 5:
            violations += 1
        if hbd is not None and float(hbd) > 5:
            violations += 1
        if hba is not None and float(hba) > 10:
            violations += 1
    except (ValueError, TypeError):
        return None, "unknown"
    return violations, ("pass" if violations <= 1 else "fail")


def build_pains_catalog():
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    return FilterCatalog.FilterCatalog(params)


def check_pains(smiles, catalog):
    if not smiles:
        return "no_smiles", ""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "invalid_smiles", ""
    matches = catalog.GetMatches(mol)
    if matches:
        alerts = "; ".join(m.GetDescription() for m in matches)
        return "fail", alerts
    return "pass", ""


def main():
    ap = argparse.ArgumentParser(description="Lipinski + PAINS filter for a ligand shortlist")
    ap.add_argument("--manifest", type=Path, default=Path("output/combined/ligand_manifest.csv"))
    ap.add_argument("--properties", type=Path, default=Path("output/compound_properties.csv"),
                     help="compound_properties.csv from extract.py")
    ap.add_argument("--names", nargs="+", default=None,
                     help="Compound names to check (default: all 100 in manifest)")
    ap.add_argument("--output", type=Path, default=Path("shortlist_lipinski_pains.csv"))
    args = ap.parse_args()

    manifest_rows = load_manifest(args.manifest)
    properties = load_properties(args.properties)
    catalog = build_pains_catalog()

    if args.names:
        wanted = {n.strip().lower() for n in args.names}
        manifest_rows = [r for r in manifest_rows if r["compound_name"].strip().lower() in wanted]

    results = []
    for row in manifest_rows:
        name = row["compound_name"]
        prop = properties.get(name.strip().lower(), {})

        mw = prop.get("MolecularWeight") or prop.get("MW")
        xlogp = prop.get("XLogP")
        hbd = prop.get("HBondDonorCount") or prop.get("HBD")
        hba = prop.get("HBondAcceptorCount") or prop.get("HBA")
        rot = prop.get("RotatableBondCount") or prop.get("RotatableBonds")

        violations, lipinski_result = lipinski_check(mw, xlogp, hbd, hba)
        pains_result, pains_alerts = check_pains(row.get("canonical_smiles", ""), catalog)

        results.append({
            "ligand_id": row.get("index", ""),
            "compound_name": name,
            "MW": mw or "",
            "XLogP": xlogp or "",
            "HBD": hbd or "",
            "HBA": hba or "",
            "RotatableBonds": rot or "",
            "lipinski_violations": violations if violations is not None else "unknown",
            "lipinski_result": lipinski_result,
            "pains_result": pains_result,
            "pains_alerts": pains_alerts,
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "ligand_id", "compound_name", "MW", "XLogP", "HBD", "HBA", "RotatableBonds",
            "lipinski_violations", "lipinski_result", "pains_result", "pains_alerts",
        ])
        writer.writeheader()
        writer.writerows(results)

    print(f"Wrote {len(results)} compounds to {args.output}\n")
    print(f"{'Compound':<25}{'MW':<10}{'XLogP':<10}{'Lipinski':<12}{'PAINS':<10}")
    for r in results:
        print(f"{r['compound_name']:<25}{str(r['MW']):<10}{str(r['XLogP']):<10}{r['lipinski_result']:<12}{r['pains_result']:<10}")

    pains_fails = [r for r in results if r["pains_result"] == "fail"]
    if pains_fails:
        print(f"\nWARNING: {len(pains_fails)} compound(s) flagged by PAINS filter:")
        for r in pains_fails:
            print(f"  - {r['compound_name']}: {r['pains_alerts']}")


if __name__ == "__main__":
    main()
