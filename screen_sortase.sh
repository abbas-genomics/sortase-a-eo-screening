#!/bin/bash
# ============================================================
# Virtual Screening Script for Sortase A (2KID)
# Docking ligands 1-100 against the same receptor
# ============================================================

NUM_LIGANDS=100
RECEPTOR_TAG="2KID"
RECEPTOR_PDBQT="protein/receptor.pdbqt"

# Grid box — catalytic triad centroid (His120/Cys184/Arg197), 22 A per side
CENTER_X=-8.081
CENTER_Y=-3.701
CENTER_Z=3.811
SIZE_X=22
SIZE_Y=22
SIZE_Z=22

EXHAUSTIVENESS=8
NUM_MODES=9
ENERGY_RANGE=3
CPU=4

mkdir -p virtual_screening/results
mkdir -p virtual_screening/logs

echo "=========================================="
echo "  Virtual Screening for Sortase A (${RECEPTOR_TAG})"
echo "=========================================="
echo "Date: $(date)"
echo "Working Directory: $(pwd)"
echo ""

# Check if vina is available
if ! command -v vina &> /dev/null; then
    echo "ERROR: vina not found! Please install AutoDock Vina."
    echo "Try: sudo apt install autodock-vina"
    exit 1
fi

# Check if receptor file exists
if [ ! -f "$RECEPTOR_PDBQT" ]; then
    echo "ERROR: $RECEPTOR_PDBQT not found!"
    exit 1
fi

SUCCESS=0
FAILED=0
FAILED_LIST=()

# Loop through ligands 1 to NUM_LIGANDS
for i in $(seq 1 "$NUM_LIGANDS"); do
    echo "----------------------------------------"
    echo "Processing Ligand $i / $NUM_LIGANDS..."

    LIGAND_PDBQT="ligands/${i}.pdbqt"
    if [ ! -f "$LIGAND_PDBQT" ]; then
        echo "WARNING: $LIGAND_PDBQT not found! Skipping..."
        FAILED=$((FAILED+1))
        FAILED_LIST+=("$i")
        continue
    fi

    # Create temporary config file
    cat > temp_config.txt << EOF
# Configuration for Ligand $i vs ${RECEPTOR_TAG}
receptor = $RECEPTOR_PDBQT
ligand = $LIGAND_PDBQT

center_x = $CENTER_X
center_y = $CENTER_Y
center_z = $CENTER_Z

size_x = $SIZE_X
size_y = $SIZE_Y
size_z = $SIZE_Z

out = virtual_screening/results/ligand_${i}_output.pdbqt

exhaustiveness = $EXHAUSTIVENESS
num_modes = $NUM_MODES
energy_range = $ENERGY_RANGE
cpu = $CPU
EOF

    LOG_FILE="virtual_screening/logs/ligand_${i}_log.txt"
    echo "Running: vina --config temp_config.txt"
    if vina --config temp_config.txt > "$LOG_FILE" 2>&1; then
        echo "✓ Ligand $i completed successfully"
        SUCCESS=$((SUCCESS+1))
    else
        echo "✗ ERROR: Ligand $i failed! Check $LOG_FILE"
        FAILED=$((FAILED+1))
        FAILED_LIST+=("$i")
    fi

done

# Clean up
rm -f temp_config.txt

echo "----------------------------------------"
echo ""
echo "=========================================="
echo "  Virtual Screening Complete!"
echo "=========================================="
echo "Succeeded: $SUCCESS / $NUM_LIGANDS"
echo "Failed:    $FAILED / $NUM_LIGANDS"
if [ "$FAILED" -gt 0 ]; then
    echo "Failed ligand IDs: ${FAILED_LIST[*]}"
    echo "(Re-run these individually once you've diagnosed the cause.)"
fi
echo "Results: virtual_screening/results/"
echo "Logs:    virtual_screening/logs/"
echo "=========================================="
