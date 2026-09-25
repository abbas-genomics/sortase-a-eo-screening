#!/usr/bin/env bash
###############################################################################
# prepare_complexes_2kid.sh
# Build receptor-ligand complex PDBs from Vina docking results
# (Sortase A / 2KID, chain A only, substrate peptide removed)
###############################################################################
set -euo pipefail

# ============================================================
RECEPTOR_PDBQT="protein/receptor.pdbqt"        # receptor used in screen_sortase.sh
RESULTS_DIR="virtual_screening/results"        # where screen_sortase.sh wrote ligand_<i>_output.pdbqt
OUT_DIR="virtual_screening/complexes"          # single flat output folder
NUM_LIGANDS=100
RECEPTOR_TAG="2KID"                             # used in output filenames
# ============================================================

mkdir -p "$OUT_DIR"

# --- Convert receptor once, reused for every complex -----------------------
RECEPTOR_PDB="$OUT_DIR/receptor_${RECEPTOR_TAG}.pdb"
obabel -ipdbqt "$RECEPTOR_PDBQT" -opdb -O "$RECEPTOR_PDB"

# --- Start the combined complex with the receptor ---------------------------
COMBINED_PDB="$OUT_DIR/combined_ligands_${RECEPTOR_TAG}_complex.pdb"
grep -v "^END" "$RECEPTOR_PDB" > "$COMBINED_PDB"
echo "TER" >> "$COMBINED_PDB"

MISSING=()

# --- Process each ligand -----------------------------------------------------
for i in $(seq 1 "$NUM_LIGANDS"); do
    VINA_OUT="${RESULTS_DIR}/ligand_${i}_output.pdbqt"

    if [ ! -f "$VINA_OUT" ]; then
        echo "WARNING: $VINA_OUT not found — skipping ligand $i"
        MISSING+=("$i")
        continue
    fi

    # Vina writes poses in order of lowest binding energy first, so model 1
    # is always the best pose. -f 1 -l 1 tells obabel to convert ONLY model 1
    # (no vina_split needed at all).
    BEST_POSE_PDB="$OUT_DIR/ligand${i}_best_pose.pdb"
    obabel -ipdbqt "$VINA_OUT" -opdb -O "$BEST_POSE_PDB" -f 1 -l 1

    # --- Per-ligand complex: receptor + this ligand's best pose ------------
    COMPLEX_PDB="$OUT_DIR/ligand${i}_${RECEPTOR_TAG}_complex.pdb"
    {
        grep -v "^END" "$RECEPTOR_PDB"
        echo "TER"
        grep -v "^END\|^HEADER\|^COMPND" "$BEST_POSE_PDB"
        echo "END"
    } > "$COMPLEX_PDB"
    echo "Created: $COMPLEX_PDB"

    # --- Append this ligand's best pose into the combined complex ----------
    grep -v "^END\|^HEADER\|^COMPND" "$BEST_POSE_PDB" >> "$COMBINED_PDB"
    echo "TER" >> "$COMBINED_PDB"
done

echo "END" >> "$COMBINED_PDB"
echo "Created combined complex: $COMBINED_PDB"

if [ "${#MISSING[@]}" -gt 0 ]; then
    echo ""
    echo "NOTE: ${#MISSING[@]} ligand(s) had no docking output and were skipped: ${MISSING[*]}"
fi
