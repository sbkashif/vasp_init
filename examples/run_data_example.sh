#!/usr/bin/env zsh
# Example run script using files in ./data
#
# Prereqs:
#   - pip install -e .
#   - Python 3 available as `python3`
#
# This will:
#   1) Merge ions from the provided PDB into the framework POSCAR
#   2) Add an NH3 molecule placed at the midpoint between two demo points
#
# Outputs:
#   - data/POSCAR_with_ions.vasp
#   - data/POSCAR_with_ions_NH3.vasp

set -euo pipefail

# Resolve repo root (directory containing this script is ./examples)
SCRIPT_DIR="${0:a:h}"
ROOT_DIR="${SCRIPT_DIR:h}"
DATA_DIR="${ROOT_DIR}/data"

POSCAR_IN="${DATA_DIR}/Framework_0_initial.vasp"
PDB_IN="${DATA_DIR}/Movie_LTA5A_crop_1.1.1_473.000000_0.000000_component_sodium_0.pdb"
DEF_IN="${DATA_DIR}/NH3_TraPPE.def"

OUT_IONS="${DATA_DIR}/POSCAR_with_ions.vasp"
OUT_FINAL="${DATA_DIR}/POSCAR_with_ions_NH3.vasp"


# Demo Cartesian points for placing NH3 (Å)
# Adjust if you need a different location
# Example: Place NH3 at the midpoint plus 1.5 Å offset in the C direction (from point 1 to point 2)
X1=6.0560
Y1=2.9715
Z1=0.0000
X2=6.2220
Y2=9.8465
Z2=0.0000
PLACE=midpoint   # midpoint|first|second
OFFSET_FROM_MIDPOINT=$((12.277/2.0))  # Ångström offset from midpoint (zsh arithmetic expansion)
OFFSET_DIRECTION=+        # + (toward point 2), - (toward point 1)

# Example: To use custom axis-aligned offsets (e.g., +1.0 Å in x, -2.0 Å in y, +0.5 Å in z), uncomment and set:
# OFFSET_X=1.0
# OFFSET_Y=-2.0
# OFFSET_Z=0.5
# and add these to the python command below:
#   --offset-x "${OFFSET_X}" --offset-y "${OFFSET_Y}" --offset-z "${OFFSET_Z}"

# Optional: selective dynamics flags
ION_FLAGS=${ION_FLAGS:-TTT}     # e.g., TTT, FFT, TFT
NH3_FLAGS=${NH3_FLAGS:-TTT}
FRAMEWORK_FLAGS=${FRAMEWORK_FLAGS:-FFF}  # e.g., FFF to fix framework atoms

# Optional: wrapping control
NO_WRAP_IONS=${NO_WRAP_IONS:-0}   # set to 1 to disable wrapping ions
NO_WRAP_NH3=${NO_WRAP_NH3:-0}     # set to 1 to disable wrapping NH3

[[ -f "${POSCAR_IN}" ]] || { echo "Missing POSCAR: ${POSCAR_IN}" >&2; exit 1; }
[[ -f "${PDB_IN}" ]]    || { echo "Missing PDB: ${PDB_IN}" >&2; exit 1; }
[[ -f "${DEF_IN}" ]]    || { echo "Missing DEF: ${DEF_IN}" >&2; exit 1; }

echo "[1/2] Merging ions into POSCAR..."
python3 "${ROOT_DIR}/examples/add_ions.py" \
  --poscar "${POSCAR_IN}" \
  --pdb "${PDB_IN}" \
  --out "${OUT_IONS}" \
  --model-index -1 \
  --ion-flags "${ION_FLAGS}" \
  --framework-flags "${FRAMEWORK_FLAGS}" \
  $( (( NO_WRAP_IONS == 1 )) && echo "--no-wrap" )


echo "[2/2] Adding NH3 molecule..."
PYTHON_NH3_CMD=(
  python3 "${ROOT_DIR}/examples/add_nh3.py"
  --poscar "${OUT_IONS}"
  --def "${DEF_IN}"
  --x1 "${X1}" --y1 "${Y1}" --z1 "${Z1}"
  --x2 "${X2}" --y2 "${Y2}" --z2 "${Z2}"
  --place "${PLACE}"
  --offset-from-midpoint "${OFFSET_FROM_MIDPOINT}"
  --offset-direction "${OFFSET_DIRECTION}"
  --flags "${NH3_FLAGS}"
  $( (( NO_WRAP_NH3 == 1 )) && echo "--no-wrap" )
  --out "${OUT_FINAL}"
)
# Uncomment the next lines to use custom axis-aligned offsets:
# PYTHON_NH3_CMD+=(--offset-x "${OFFSET_X}" --offset-y "${OFFSET_Y}" --offset-z "${OFFSET_Z}")

"${PYTHON_NH3_CMD[@]}"

echo "Done. Outputs:"
echo "  Ions merged: ${OUT_IONS}"
echo "  Final with NH3: ${OUT_FINAL}"
