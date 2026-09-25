#!/usr/bin/env python3
"""
combine.py – Combine optimized ligands for virtual screening
=============================================================
Copyright (c) Abbas Aliyu. All rights reserved.
Developed by Abbas Aliyu Abbas

Merges all optimized structures into a single file for AutoDock/YASARA docking runs.
"""

import argparse
import csv
import logging
import re
import subprocess
from pathlib import Path
from rdkit import Chem


# ─────────────────────────────────────────────────────────────
# SETUP
# ─────────────────────────────────────────────────────────────

def setup_logger(log_path: Path):
    """Setup dual file + console logging."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    
    # File handler
    fh = logging.FileHandler(log_path)
    fh.setLevel(logging.DEBUG)
    
    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    
    # Format
    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)-5s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger


# ─────────────────────────────────────────────────────────────
# CORE FUNCTIONS
# ─────────────────────────────────────────────────────────────

def read_renaming_log(log_path: Path) -> dict:
    """
    Read renaming_log.csv and return dict mapping numbered_file → {name, cid, smiles}.
    File format: index,numbered_file,original_file,cid,input_name,structure_type,source,CanonicalSMILES,IsomericSMILES
    """
    ligand_map = {}
    if not log_path.exists():
        return ligand_map
    
    with open(log_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            numbered_file = row.get('numbered_file', '')
            ligand_map[numbered_file] = {
                'name': row.get('input_name', ''),
                'cid': row.get('cid', ''),
                'structure_type': row.get('structure_type', 'unknown'),
                'canonical_smiles': row.get('CanonicalSMILES', ''),
                'isomeric_smiles': row.get('IsomericSMILES', ''),
            }
    
    return ligand_map


def find_sdf_files(input_dir: Path, sort_numeric: bool = True) -> list:
    """Find all SDF files in directory, optionally sort numerically."""
    sdf_files = list(input_dir.glob("*.sdf"))
    
    if sort_numeric:
        # Sort by leading number: 1.sdf, 2.sdf, ..., 33.sdf
        sdf_files.sort(
            key=lambda p: int(re.match(r'(\d+)', p.name).group(1)) 
            if re.match(r'(\d+)', p.name) else 0
        )
    else:
        sdf_files.sort()
    
    return sdf_files


def combine_sdf_files(sdf_files: list, output_file: Path, logger: logging.Logger) -> int:
    """
    Combine multiple SDF files into one.
    Each molecule retains its name and conformer info.
    Returns number of molecules combined.
    """
    combined_writer = Chem.SDWriter(str(output_file))
    count = 0
    
    for sdf_file in sdf_files:
        try:
            supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=False)
            for mol in supplier:
                if mol is not None:
                    combined_writer.write(mol)
                    count += 1
        except Exception as e:
            logger.warning(f"Error reading {sdf_file.name}: {e}")
    
    combined_writer.close()
    return count


def combine_mol2_via_obabel(sdf_file: Path, output_file: Path, logger: logging.Logger) -> bool:
    """Convert combined SDF to MOL2 using OpenBabel."""
    try:
        cmd = [
            "obabel",
            "-isdf", str(sdf_file),
            "-omol2",
            "-O", str(output_file)
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        
        if result.returncode != 0:
            logger.warning(f"OpenBabel conversion failed: {result.stderr}")
            return False
        
        return True
    except FileNotFoundError:
        logger.warning("OpenBabel not found. Skipping MOL2 conversion.")
        return False


def create_manifest(sdf_files: list, output_dir: Path, ligand_map: dict, logger: logging.Logger):
    """
    Create manifest CSV with ligand info for reference during docking.
    Useful for tracking which compound is which in YASARA.
    """
    manifest_path = output_dir / "ligand_manifest.csv"
    
    with open(manifest_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'index', 'numbered_file', 'compound_name', 'cid', 
            'structure_type', 'canonical_smiles', 'isomeric_smiles'
        ])
        writer.writeheader()
        
        for idx, sdf_file in enumerate(sdf_files, start=1):
            info = ligand_map.get(sdf_file.name, {})
            writer.writerow({
                'index': idx,
                'numbered_file': sdf_file.name,
                'compound_name': info.get('name', 'unknown'),
                'cid': info.get('cid', ''),
                'structure_type': info.get('structure_type', ''),
                'canonical_smiles': info.get('canonical_smiles', ''),
                'isomeric_smiles': info.get('isomeric_smiles', ''),
            })
    
    logger.info(f"Manifest saved: {manifest_path}")


# ─────────────────────────────────────────────────────────────
# MAIN PIPELINE
# ─────────────────────────────────────────────────────────────

def run_pipeline(input_dir: Path, output_dir: Path, renaming_log: Path,
                 convert_mol2: bool, logger: logging.Logger) -> int:
    """Main combining pipeline."""
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find SDF files
    sdf_files = find_sdf_files(input_dir)
    
    if not sdf_files:
        logger.error(f"No SDF files found in {input_dir}")
        return 1
    
    logger.info(f"Found {len(sdf_files)} SDF files to combine")
    
    # Read ligand metadata
    ligand_map = read_renaming_log(renaming_log)
    
    # Combine SDF files
    combined_sdf = output_dir / "combined_ligands.sdf"
    logger.info(f"Combining {len(sdf_files)} molecules...")
    
    mol_count = combine_sdf_files(sdf_files, combined_sdf, logger)
    
    if mol_count == 0:
        logger.error("No molecules were combined successfully")
        return 1
    
    logger.info(f"✓ Combined {mol_count} molecules → {combined_sdf}")
    
    # Optional MOL2 conversion from combined SDF
    combined_mol2 = output_dir / "combined_ligands.mol2"
    if convert_mol2:
        logger.info("Converting to MOL2 format...")
        if combine_mol2_via_obabel(combined_sdf, combined_mol2, logger):
            logger.info(f"✓ Created MOL2 → {combined_mol2}")
        else:
            logger.warning("MOL2 conversion skipped (OpenBabel conversion failed)")
    
    # Create manifest
    create_manifest(sdf_files, output_dir, ligand_map, logger)
    
    # Summary
    logger.info("=" * 70)
    logger.info("COMBINATION COMPLETE")
    logger.info("=" * 70)
    logger.info(f"Combined file     : {combined_sdf}")
    logger.info(f"MOL2 file         : {combined_mol2 if convert_mol2 else 'not requested'}")
    logger.info(f"Molecules included: {mol_count}")
    logger.info(f"Manifest file     : {output_dir / 'ligand_manifest.csv'}")
    logger.info("")
    logger.info("Next step: Use combined_ligands.sdf (and combined_ligands.mol2 if generated)")
    logger.info("Reference ligand_manifest.csv to identify compounds during docking")
    
    return 0


# ─────────────────────────────────────────────────────────────
# CLI ENTRY POINT
# ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Combine optimized ligands into single file for virtual screening",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Combine MMFF94 optimized structures
  python combine.py --input-dir output/optimized_sdf --output-dir output/combined

  # Combine also generate MOL2 format
  python combine.py --input-dir output/optimized_sdf --output-dir output/combined --mol2

  # Use custom renaming log location
  python combine.py --input-dir output/optimized_sdf --renaming-log output/custom/renaming_log.csv
        """
    )
    
    parser.add_argument(
        "--input-dir",
        type=Path,
        default="output/optimized_sdf",
        help="Directory with optimized SDF files (default: output/optimized_sdf)"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default="output/combined",
        help="Output directory for combined ligand files (default: output/combined)"
    )
    parser.add_argument(
        "--renaming-log",
        type=Path,
        default="output/numbered_sdf/renaming_log.csv",
        help="Path to renaming_log.csv for ligand metadata (default: output/numbered_sdf/renaming_log.csv)"
    )
    parser.add_argument(
        "--mol2",
        action="store_true",
        help="Also generate combined_ligands.mol2 in output directory"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    log_path = args.output_dir / "combine.log"
    logger = setup_logger(log_path)
    
    logger.info("=" * 70)
    logger.info("LIGAND COMBINATION FOR VIRTUAL SCREENING")
    logger.info("=" * 70)
    logger.info(f"Input directory  : {args.input_dir}")
    logger.info(f"Output directory : {args.output_dir}")
    logger.info(f"Renaming log     : {args.renaming_log}")
    logger.info(f"MOL2 conversion  : {'ON' if args.mol2 else 'OFF'}")
    logger.info("")
    
    # Run pipeline
    returncode = run_pipeline(
        args.input_dir,
        args.output_dir,
        args.renaming_log,
        args.mol2,
        logger
    )
    
    logger.info("")
    return returncode


if __name__ == "__main__":
    exit(main())
