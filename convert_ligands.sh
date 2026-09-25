#!/usr/bin/env bash
###############################################################################
# convert_ligands.sh
# Convert individually optimized MOL2 ligands (from optimize.py) into
# individual PDBQT files, as expected by screen_sortase.sh (ligands/N.pdbqt)
###############################################################################
"""
=======================================
Copyright (c) Abbas Aliyu. All rights reserved.
Developed by Abbas Aliyu Abbas
"""

set -euo pipefail

MOL2_DIR="output/optimized/optimized_mol2"   # per-ligand mol2 files from optimize.py
LIGANDS_DIR="ligands"                        # where screen_sortase.sh expects N.pdbqt

mkdir -p "$LIGANDS_DIR"

if ! command -v obabel &> /dev/null; then
    echo "ERROR: obabel not found! Install OpenBabel first."
    exit 1
fi

if [ ! -d "$MOL2_DIR" ]; then
    echo "ERROR: $MOL2_DIR not found. Check that optimize.py has been run."
    exit 1
fi

count=0
fail=0
failed_list=()

for f in "$MOL2_DIR"/*.mol2; do
    [ -e "$f" ] || { echo "No mol2 files found in $MOL2_DIR"; exit 1; }
    n=$(basename "$f" .mol2)
    if obabel "$f" -opdbqt -O "$LIGANDS_DIR/${n}.pdbqt" 2>> "$LIGANDS_DIR/convert.log"; then
        count=$((count+1))
    else
        echo "WARNING: conversion failed for $n"
        fail=$((fail+1))
        failed_list+=("$n")
    fi
done

echo "----------------------------------------"
echo "Converted $count ligand(s) to PDBQT in $LIGANDS_DIR/"
if [ "$fail" -gt 0 ]; then
    echo "WARNING: $fail ligand(s) failed conversion: ${failed_list[*]}"
    echo "See $LIGANDS_DIR/convert.log for details."
fi
