#!/usr/bin/env bash
# Standalone workflow script for running vasp_init from any location
#
# This script can be run from any directory. You just need to:
#   1) Have vasp_init package installed (pip install vasp_init or pip install -e /path/to/vasp_init)
#   2) Set the file paths below to your input files
#   3) Run: bash run_workflow_custom.sh
#
# The script will create output files in the current working directory.

set -euo pipefail

# ============================================================================
# USER CONFIGURATION - Edit these paths to match your files
# ============================================================================

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/../data"

# Input files (use absolute paths or relative to where you run this script)
POSCAR_IN="${DATA_DIR}/Framework_0_initial.vasp"
PDB_IN="${DATA_DIR}/Movie_LTA5A_crop_1.1.1_473.000000_0.000000_component_sodium_0.pdb"
DEF_IN="${DATA_DIR}/NH3_TraPPE.def"

# Output files (will be created in current directory)
OUT_IONS="./POSCAR_with_ions.vasp"
OUT_FINAL="./POSCAR_with_ions_NH3.vasp"

# NH3 placement coordinates (Cartesian Å)
X1=6.0560
Y1=2.9715
Z1=0.0000
X2=6.2220
Y2=9.8465
Z2=0.0000
PLACE=midpoint   # midpoint|first|second

# Offset settings
OFFSET_FROM_MIDPOINT=0.0  # Offset along connecting vector (Å)
OFFSET_DIRECTION=+        # + (toward point 2) or - (toward point 1)

# Custom axis-aligned offsets (applied after base placement)
OFFSET_X=0.0
OFFSET_Y=0.0
OFFSET_Z=18.41625  # Example: 24.555 * 0.75 for 75% of cell height

# Selective dynamics flags
ION_FLAGS="TTT"           # T=movable, F=fixed (e.g., TTT, FFT, TFT)
NH3_FLAGS="TTT"
FRAMEWORK_FLAGS="FFF"     # FFF = fix framework atoms

# Wrapping control (set to "--no-wrap" to disable wrapping, or "" to enable)
NO_WRAP_IONS=""           # "" = wrap, "--no-wrap" = don't wrap
NO_WRAP_NH3=""

# ============================================================================
# END USER CONFIGURATION
# ============================================================================

# Verify input files exist
[[ -f "${POSCAR_IN}" ]] || { echo "ERROR: Missing POSCAR: ${POSCAR_IN}" >&2; exit 1; }
[[ -f "${PDB_IN}" ]]    || { echo "ERROR: Missing PDB: ${PDB_IN}" >&2; exit 1; }
[[ -f "${DEF_IN}" ]]    || { echo "ERROR: Missing DEF: ${DEF_IN}" >&2; exit 1; }

# Verify vasp_init is installed
if ! python3 -c "import vasp_init" 2>/dev/null; then
    echo "ERROR: vasp_init package not found. Please install it:" >&2
    echo "  pip install vasp_init" >&2
    echo "  or: pip install -e /path/to/vasp_init" >&2
    exit 1
fi

echo "============================================================================"
echo "VASP Initialization Workflow"
echo "============================================================================"
echo "Working directory: $(pwd)"
echo "Input POSCAR: ${POSCAR_IN}"
echo "Input PDB: ${PDB_IN}"
echo "Input DEF: ${DEF_IN}"
echo "Output (ions): ${OUT_IONS}"
echo "Output (final): ${OUT_FINAL}"
echo "============================================================================"

echo ""
echo "[1/2] Merging ions from PDB into POSCAR..."
python3 -m vasp_init.cli.add_ions \
  --poscar "${POSCAR_IN}" \
  --pdb "${PDB_IN}" \
  --out "${OUT_IONS}" \
  --model-index -1 \
  --ion-flags "${ION_FLAGS}" \
  --framework-flags "${FRAMEWORK_FLAGS}" \
  ${NO_WRAP_IONS}

echo ""
echo "[2/2] Adding NH3 molecule..."
python3 -m vasp_init.cli.add_nh3 \
  --poscar "${OUT_IONS}" \
  --def "${DEF_IN}" \
  --x1 "${X1}" --y1 "${Y1}" --z1 "${Z1}" \
  --x2 "${X2}" --y2 "${Y2}" --z2 "${Z2}" \
  --place "${PLACE}" \
  --offset-from-midpoint "${OFFSET_FROM_MIDPOINT}" \
  --offset-direction "${OFFSET_DIRECTION}" \
  --offset-x "${OFFSET_X}" \
  --offset-y "${OFFSET_Y}" \
  --offset-z "${OFFSET_Z}" \
  --flags "${NH3_FLAGS}" \
  ${NO_WRAP_NH3} \
  --out "${OUT_FINAL}"

echo ""
echo "============================================================================"
echo "Done! Outputs:"
echo "  Ions merged: ${OUT_IONS}"
echo "  Final with NH3: ${OUT_FINAL}"
echo "============================================================================"
