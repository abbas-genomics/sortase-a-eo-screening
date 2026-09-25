"""
optimize.py  –  Geometry optimization of SDF structures for docking / DFT prep
=================================================================================
Copyright (c) Abbas Aliyu. All rights reserved.
Developed by Abbas Aliyu Abbas

Reads numbered SDF files from the numbered_sdf folder produced by extract.py,
optimizes each molecule and writes:
    output_dir/optimized_sdf/   <numbered>.sdf
    output_dir/optimized_mol2/  <numbered>.mol2
    output_dir/optimization.log
    output_dir/optimization_results.csv

Four optimization methods (user's choice via --method):
  1. mmff94   – MMFF94 force field (best for drug-like / docking)
  2. mmff94s  – MMFF94s static variant (better for aromatics / conjugated)
  3. uff      – Universal Force Field (fallback; works on any element)
  4. etkdg    – ETKDGv3 conformer embedding only (ideal starting geom for DFT)

Usage:
    python optimize.py [--input-dir output/numbered_sdf]
                       [--output-dir output/optimized]
                       [--method mmff94]
                       [--steps 2000]
                       [--rmsg 0.001]
"""

import argparse
import csv
import logging
import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# ── RDKit imports ──────────────────────────────────────────────────────────────
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    from rdkit.Chem.rdForceFieldHelpers import (
        MMFFGetMoleculeForceField,
        MMFFGetMoleculeProperties,
        UFFGetMoleculeForceField,
    )
    from rdkit.Chem.rdmolfiles import MolToMolBlock, SDWriter
except ImportError as exc:
    raise SystemExit("RDKit is required: pip install rdkit-pypi") from exc

OBABEL = "obabel"  # assumes obabel is on PATH
COMMON_ATOMS = {"C", "H", "O", "N", "S", "P", "F", "Cl", "Br", "I"}

# ──────────────────────────────────────────────────────────────────────────────

def sniff_method_from_smiles(smiles: str) -> str:
    """Auto-detect best optimization method from SMILES.
    
    Returns:
        'mmff94s' if aromatic rings detected
        'uff' if unusual atoms detected
        'mmff94' for standard organic molecules
    """
    if not smiles or not isinstance(smiles, str):
        return "mmff94"  # default

    # Try to parse with RDKit to detect aromaticity and atoms.
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "mmff94"  # default if parsing fails

        # Check for aromatic atoms.
        for atom in mol.GetAtoms():
            if atom.GetIsAromatic():
                return "mmff94s"  # Better for aromatics
            # Check for unusual atoms.
            symbol = atom.GetSymbol()
            if symbol not in COMMON_ATOMS:
                return "uff"  # Unusual atom – use UFF

        return "mmff94"  # standard organic
    except Exception:
        # If RDKit parsing fails, fall back to regex heuristic.
        if re.search(r"[c,n,o,s,p]", smiles):
            return "mmff94s"
        return "mmff94"


def setup_logger(log_file: Path) -> logging.Logger:
    logger = logging.getLogger("optimize_pipeline")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setFormatter(fmt)
    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger


def write_csv(path: Path, rows: List[Dict[str, Any]], fieldnames: List[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def sdf_to_mol2(sdf_path: Path, mol2_path: Path) -> Tuple[bool, str]:
    """Use OpenBabel CLI to convert SDF → mol2."""
    try:
        result = subprocess.run(
            [OBABEL, "-isdf", str(sdf_path), "-omol2", "-O", str(mol2_path)],
            capture_output=True,
            timeout=60,
        )
        if result.returncode == 0 and mol2_path.exists():
            return True, ""
        return False, result.stderr.decode(errors="replace").strip()
    except FileNotFoundError:
        return False, "obabel not found on PATH"
    except subprocess.TimeoutExpired:
        return False, "obabel conversion timed out"


def embed_3d(mol) -> Optional[Any]:
    """Add Hs and embed 3D if not already 3D."""
    mol = Chem.AddHs(mol)
    conf_ids = mol.GetConformers()
    if conf_ids:
        # Check if conformer already has 3D coordinates.
        conf = mol.GetConformer(0)
        pos = conf.GetPositions()
        if any(p[2] != 0.0 for p in pos):
            return mol  # already 3D

    params = AllChem.ETKDGv3()
    params.randomSeed = 0xDEAD
    code = AllChem.EmbedMolecule(mol, params)
    if code == -1:
        # Fallback to distance geometry.
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if not mol.GetConformers():
        return None
    return mol


def optimize_mmff(mol, variant: str = "MMFF94", max_iters: int = 2000,
                   rmsg: float = 0.001) -> Tuple[Optional[Any], float, float, str]:
    """Optimize with MMFF94 or MMFF94s. Returns (mol, e_before, e_after, info)."""
    mol = embed_3d(mol)
    if mol is None:
        return None, 0.0, 0.0, "3D embedding failed"

    props = MMFFGetMoleculeProperties(mol, mmffVariant=variant)
    if props is None:
        return None, 0.0, 0.0, f"{variant} properties unavailable (unsupported atom type?)"

    ff = MMFFGetMoleculeForceField(mol, props)
    if ff is None:
        return None, 0.0, 0.0, f"{variant} force field setup failed"

    e_before = ff.CalcEnergy()
    ff.Minimize(maxIts=max_iters, forceTol=rmsg)
    e_after = ff.CalcEnergy()
    return mol, e_before, e_after, "OK"


def optimize_uff(mol, max_iters: int = 2000,
                  rmsg: float = 0.001) -> Tuple[Optional[Any], float, float, str]:
    """Optimize with Universal Force Field."""
    mol = embed_3d(mol)
    if mol is None:
        return None, 0.0, 0.0, "3D embedding failed"

    ff = UFFGetMoleculeForceField(mol)
    if ff is None:
        return None, 0.0, 0.0, "UFF force field setup failed"

    e_before = ff.CalcEnergy()
    ff.Minimize(maxIts=max_iters, forceTol=rmsg)
    e_after = ff.CalcEnergy()
    return mol, e_before, e_after, "OK"


def optimize_etkdg(mol) -> Tuple[Optional[Any], float, float, str]:
    """ETKDGv3 embedding only – no force field minimisation.
    Energy reported as MMFF94 single-point to allow comparison."""
    mol = embed_3d(mol)
    if mol is None:
        return None, 0.0, 0.0, "3D embedding failed"

    # Report single-point MMFF energy on the embedded geometry.
    props = MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
    if props is not None:
        ff = MMFFGetMoleculeForceField(mol, props)
        if ff:
            e = ff.CalcEnergy()
            return mol, e, e, "OK (ETKDG embedding – MMFF SP energy reported, no minimisation)"

    return mol, 0.0, 0.0, "OK (ETKDG embedding – energy unavailable)"


def load_sdf_mol(sdf_path: Path) -> Optional[Any]:
    """Load first molecule from SDF file."""
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
    for mol in suppl:
        if mol is not None:
            return mol
    # Try without H removal.
    suppl2 = Chem.SDMolSupplier(str(sdf_path), removeHs=True, sanitize=False)
    for mol in suppl2:
        if mol is not None:
            try:
                Chem.SanitizeMol(mol)
                return mol
            except Exception:
                return mol
    return None


def run_optimization(
    sdf_path: Path,
    method: str,
    max_iters: int,
    rmsg: float,
) -> Tuple[Optional[Any], float, float, str]:
    """Dispatch to the chosen optimization method."""
    mol = load_sdf_mol(sdf_path)
    if mol is None:
        return None, 0.0, 0.0, "Could not parse SDF file"

    method = method.lower()
    if method == "mmff94":
        return optimize_mmff(mol, variant="MMFF94", max_iters=max_iters, rmsg=rmsg)
    if method == "mmff94s":
        return optimize_mmff(mol, variant="MMFF94s", max_iters=max_iters, rmsg=rmsg)
    if method == "uff":
        return optimize_uff(mol, max_iters=max_iters, rmsg=rmsg)
    if method == "etkdg":
        return optimize_etkdg(mol)

    return None, 0.0, 0.0, f"Unknown method: {method}"


# ──────────────────────────────────────────────────────────────────────────────

def read_renaming_log(log_csv: Path) -> Dict[str, Dict[str, str]]:
    """Return mapping numbered_file -> {name, smiles} from renaming_log.csv."""
    mapping: Dict[str, Dict[str, str]] = {}
    if not log_csv.exists():
        return mapping
    with log_csv.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            nf = row.get("numbered_file", "").strip()
            if nf:
                mapping[nf] = {
                    "name": row.get("input_name", "").strip(),
                    "smiles": row.get("CanonicalSMILES", "").strip(),
                }
    return mapping


def run_pipeline(
    input_dir: Path,
    output_dir: Path,
    method: Optional[str],
    auto: bool,
    max_iters: int,
    rmsg: float,
) -> int:
    sdf_out = output_dir / "optimized_sdf"
    mol2_out = output_dir / "optimized_mol2"
    sdf_out.mkdir(parents=True, exist_ok=True)
    mol2_out.mkdir(parents=True, exist_ok=True)

    log_file = output_dir / "optimization.log"
    results_csv = output_dir / "optimization_results.csv"
    logger = setup_logger(log_file)

    logger.info("=== Optimization Pipeline ===")
    logger.info("Input dir  : %s", input_dir)
    logger.info("Auto mode  : %s", "ON" if auto else "OFF")
    if not auto:
        logger.info("Method     : %s", method.upper() if method else "none")
    logger.info("Max steps  : %d", max_iters)
    logger.info("Force tol  : %g", rmsg)

    renaming_log_csv = input_dir / "renaming_log.csv"
    smiles_map = read_renaming_log(renaming_log_csv)

    # Collect numbered SDF files (1.sdf, 2.sdf, ...)
    sdf_files = sorted(
        [f for f in input_dir.glob("*.sdf") if re.match(r"^\d+\.sdf$", f.name)],
        key=lambda p: int(p.stem),
    )

    if not sdf_files:
        logger.error("No numbered SDF files found in: %s", input_dir)
        return 1

    logger.info("Found %d structures to optimize", len(sdf_files))

    rows: List[Dict[str, Any]] = []
    done = 0
    errors = 0

    for sdf_path in sdf_files:
        num = int(sdf_path.stem)
        info_dict = smiles_map.get(sdf_path.name, {})
        compound_name = info_dict.get("name", sdf_path.stem)
        smiles = info_dict.get("smiles", "")
        
        # Auto-detect method based on SMILES if enabled.
        active_method = method
        if auto:
            active_method = sniff_method_from_smiles(smiles)

        logger.info("[%d/%d] %s – %s [method: %s]", num, len(sdf_files), sdf_path.name, compound_name, active_method.upper())

        mol, e_before, e_after, info = run_optimization(sdf_path, active_method, max_iters, rmsg)

        if mol is None or info.startswith("Could not") or "failed" in info:
            logger.error("  ERROR: %s", info)
            errors += 1
            rows.append({
                "number": num,
                "numbered_file": sdf_path.name,
                "compound_name": compound_name,
                "method": active_method.upper(),
                "energy_before_kcal_mol": "",
                "energy_after_kcal_mol": "",
                "energy_change_kcal_mol": "",
                "sdf_written": "NO",
                "mol2_written": "NO",
                "status": "ERROR",
                "info": info,
            })
            continue

        delta = e_after - e_before

        # Write optimized SDF.
        out_sdf = sdf_out / sdf_path.name
        try:
            with SDWriter(str(out_sdf)) as w:
                w.write(mol)
            sdf_ok = True
        except Exception as exc:
            sdf_ok = False
            logger.warning("  SDF write failed: %s", exc)

        # Convert to mol2 via obabel.
        out_mol2 = mol2_out / f"{sdf_path.stem}.mol2"
        mol2_ok, mol2_err = sdf_to_mol2(out_sdf if sdf_ok else sdf_path, out_mol2)
        if not mol2_ok:
            logger.warning("  mol2 conversion failed: %s", mol2_err)

        status = "DONE" if (sdf_ok or mol2_ok) else "ERROR"
        if status == "DONE":
            done += 1
            logger.info(
                "  %s | E_before=%.4f  E_after=%.4f  ΔE=%.4f kcal/mol | SDF=%s mol2=%s",
                status, e_before, e_after, delta,
                "✓" if sdf_ok else "✗", "✓" if mol2_ok else "✗",
            )
        else:
            errors += 1
            logger.error("  ERROR: SDF and mol2 both failed to write")

        rows.append({
            "number": num,
            "numbered_file": sdf_path.name,
            "compound_name": compound_name,
            "method": active_method.upper(),
            "energy_before_kcal_mol": round(e_before, 6),
            "energy_after_kcal_mol": round(e_after, 6),
            "energy_change_kcal_mol": round(delta, 6),
            "sdf_written": "YES" if sdf_ok else "NO",
            "mol2_written": "YES" if mol2_ok else "NO",
            "status": status,
            "info": info,
        })

    write_csv(
        results_csv,
        rows,
        [
            "number", "numbered_file", "compound_name", "method",
            "energy_before_kcal_mol", "energy_after_kcal_mol", "energy_change_kcal_mol",
            "sdf_written", "mol2_written", "status", "info",
        ],
    )

    logger.info("Optimization complete. DONE=%d | ERROR=%d | TOTAL=%d", done, errors, len(sdf_files))

    print("\n=== OPTIMIZATION SUMMARY ===")
    if auto:
        print(f"Mode            : AUTO (per-compound method selection)")
    else:
        print(f"Method          : {method.upper()}")
    print(f"Input dir       : {input_dir}")
    print(f"Output dir      : {output_dir}")
    print(f"Done            : {done}")
    print(f"Errors          : {errors}")
    print(f"Optimized SDF   : {sdf_out}")
    print(f"Optimized mol2  : {mol2_out}")
    print(f"Results CSV     : {results_csv}")
    print(f"Log file        : {log_file}")

    return 0 if errors == 0 else 2


# ──────────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Optimize 3D molecular structures for docking / DFT.\n\n"
            "Methods:\n"
            "  mmff94   MMFF94 force field – best for drug-like molecules & docking\n"
            "  mmff94s  MMFF94s – better for aromatic/conjugated systems\n"
            "  uff      Universal Force Field – works on any element (good fallback)\n"
            "  etkdg    ETKDGv3 embedding only – best starting geometry for DFT\n"
            "\nAuto mode:\n"
            "  --auto   Automatically selects method per compound based on SMILES:\n"
            "           - Aromatic rings detected → MMFF94s\n"
            "           - Unusual atoms detected → UFF\n"
            "           - Otherwise → MMFF94\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input-dir",
        default="output/numbered_sdf",
        help="Folder with numbered SDF files (from extract.py). Default: output/numbered_sdf",
    )
    parser.add_argument(
        "--output-dir",
        default="output/optimized",
        help="Output folder. Default: output/optimized",
    )
    parser.add_argument(
        "--auto",
        action="store_true",
        help="Auto-detect optimization method per compound based on SMILES. Overrides --method.",
    )
    parser.add_argument(
        "--method",
        default="mmff94",
        choices=["mmff94", "mmff94s", "uff", "etkdg"],
        help=(
            "Optimization method (ignored if --auto is set):\n"
            "  mmff94   – MMFF94 (recommended for docking)\n"
            "  mmff94s  – MMFF94s (recommended for aromatics)\n"
            "  uff      – Universal Force Field (works on everything)\n"
            "  etkdg    – 3D embedding only, no minimisation (for DFT input)\n"
            "Default: mmff94"
        ),
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=2000,
        help="Max minimization steps (ignored for etkdg). Default: 2000",
    )
    parser.add_argument(
        "--force-tol",
        type=float,
        default=0.001,
        help="Force convergence tolerance in kcal/mol/Å (ignored for etkdg). Default: 0.001",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    input_path = Path(args.input_dir)
    output_path = Path(args.output_dir)

    if not input_path.exists():
        print(f"Input directory not found: {input_path}")
        raise SystemExit(1)

    # Use auto mode if --auto is set, otherwise use explicit --method.
    method_to_use = None if args.auto else args.method

    raise SystemExit(
        run_pipeline(input_path, output_path, method_to_use, args.auto, args.steps, args.force_tol)
    )
