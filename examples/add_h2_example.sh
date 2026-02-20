#!/bin/bash
# Example script for adding H2 molecules to a POSCAR using vasp_init

# Configuration
DATA_DIR="../data"
POSCAR_FILE="${DATA_DIR}/Framework_0_initial.vasp"
H2_DEF_FILE="${DATA_DIR}/H2_TraPPE.def"

# Output files (in data directory, following NH3 example pattern)
OUT_H2="${DATA_DIR}/POSCAR_with_H2.vasp"
OUT_H2_OFFSET="${DATA_DIR}/POSCAR_with_H2_offset.vasp"

# H2 placement parameters
X1=2.0
Y1=0.0
Z1=0.0
X2=8.0
Y2=0.0
Z2=0.0

# Add H2 molecule at midpoint between two coordinates
echo "Adding H2 molecule between coordinates (${X1}, ${Y1}, ${Z1}) and (${X2}, ${Y2}, ${Z2})"

python -m vasp_init.cli.add_h2 \
    --poscar "${POSCAR_FILE}" \
    --def "${H2_DEF_FILE}" \
    --x1 "${X1}" --y1 "${Y1}" --z1 "${Z1}" \
    --x2 "${X2}" --y2 "${Y2}" --z2 "${Z2}" \
    --place midpoint \
    --out "${OUT_H2}"

echo "H2 molecule added successfully!"
echo "Output file: ${OUT_H2}"

# Optional: Add H2 with custom offset
echo ""
echo "Adding H2 with custom z-offset of 3.0 Å"

python -m vasp_init.cli.add_h2 \
    --poscar "${POSCAR_FILE}" \
    --def "${H2_DEF_FILE}" \
    --x1 "${X1}" --y1 "${Y1}" --z1 "${Z1}" \
    --x2 "${X2}" --y2 "${Y2}" --z2 "${Z2}" \
    --place midpoint \
    --offset-z 3.0 \
    --out "${OUT_H2_OFFSET}"

echo "H2 molecule with offset added successfully!"
echo "Output file: ${OUT_H2_OFFSET}"