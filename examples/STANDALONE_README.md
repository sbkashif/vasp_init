# Standalone Workflow Script

This script (`run_workflow_custom.sh`) is designed to be run from **any directory** on your cluster or local machine, as long as you have the `vasp_init` package installed.

## Prerequisites

1. **Install vasp_init package:**
   ```bash
   pip install vasp_init
   # or for development:
   pip install -e /path/to/vasp_init
   ```

2. **Have your input files ready:**
   - POSCAR or CONTCAR file (VASP format)
   - PDB file with ions from RASPA or other source
   - NH3 .def file (TraPPE force field format)

## Usage

1. **Copy the script to your working directory:**
   ```bash
   cp /path/to/vasp_init/examples/run_workflow_custom.sh .
   ```

2. **Edit the configuration section** at the top of the script to set your file paths:
   ```bash
   # Input files (use absolute paths or relative to where you run this script)
   POSCAR_IN="./Framework_0_initial.vasp"
   PDB_IN="./Movie_LTA5A_crop_1.1.1_473.000000_0.000000_component_sodium_0.pdb"
   DEF_IN="./NH3_TraPPE.def"
   
   # NH3 placement coordinates (Cartesian Å)
   X1=6.0560
   Y1=2.9715
   Z1=0.0000
   X2=6.2220
   Y2=9.8465
   Z2=0.0000
   
   # Custom axis-aligned offsets
   OFFSET_X=0.0
   OFFSET_Y=0.0
   OFFSET_Z=18.41625  # Example: 24.555 * 0.75 for 75% of cell height
   ```

3. **Run the script:**
   ```bash
   bash run_workflow_custom.sh
   ```

## What it does

1. **Merges ions** from the PDB file into the POSCAR framework
2. **Adds NH3 molecule** at specified location with custom offsets
3. **Creates output files** in the current directory:
   - `POSCAR_with_ions.vasp` - framework + ions
   - `POSCAR_with_ions_NH3.vasp` - framework + ions + NH3

## Configuration Options

### Placement
- `X1, Y1, Z1`: First point (Cartesian Å)
- `X2, Y2, Z2`: Second point (Cartesian Å)
- `PLACE`: Where to place NH3 (`midpoint`, `first`, or `second`)

### Offsets
- `OFFSET_FROM_MIDPOINT`: Offset along connecting vector (Å)
- `OFFSET_DIRECTION`: `+` (toward point 2) or `-` (toward point 1)
- `OFFSET_X, OFFSET_Y, OFFSET_Z`: Custom axis-aligned offsets (Å)

### Selective Dynamics
- `ION_FLAGS`: Flags for ions (e.g., `TTT` = movable, `FFF` = fixed)
- `NH3_FLAGS`: Flags for NH3 atoms
- `FRAMEWORK_FLAGS`: Flags for framework atoms (e.g., `FFF` = fixed)

### Wrapping
- `NO_WRAP_IONS`: Set to `"--no-wrap"` to disable wrapping ions
- `NO_WRAP_NH3`: Set to `"--no-wrap"` to disable wrapping NH3

## Example for Cluster Use

```bash
# On your cluster, in your job directory:
cd /scratch/username/job_001

# Copy your input files here
cp /data/poscar/Framework_0_initial.vasp .
cp /data/pdb/ions.pdb .
cp /data/molecules/NH3_TraPPE.def .

# Copy and edit the workflow script
cp /path/to/vasp_init/examples/run_workflow_custom.sh .
nano run_workflow_custom.sh  # Edit file paths and coordinates

# Run it
bash run_workflow_custom.sh

# Your output files are now in the current directory
ls -lh POSCAR_with_ions*.vasp
```

## Troubleshooting

**Error: "vasp_init package not found"**
- Make sure you've installed the package: `pip install vasp_init`
- Or install in development mode: `pip install -e /path/to/vasp_init`

**Error: "Missing POSCAR: ./Framework_0_initial.vasp"**
- The script looks for files relative to where you run it
- Either copy files to current directory or use absolute paths in the configuration

**Error: Command not found**
- Make sure the script has execute permissions: `chmod +x run_workflow_custom.sh`
- Or run with bash: `bash run_workflow_custom.sh`
