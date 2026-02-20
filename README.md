# vasp_init

Utilities to merge ions and molecules into the `Framework` files from RASPA. The starting file will be the framework file written out in the `vasp` format from a RASPA simulatins. Then, the molecule and ions may be added to it which are generally written out in PDB format by RASPA.

Terminology for the purpose of this package:
 - **Framework**: Something like Zeolite or MOF or a membrane
 - **Ions**: Typically added to neutralize the framework and are needed at the same precise locations in VASP input as in the last frame of GCMC simulation in RASPA.
 - **Molecules**: Typically the system needed to be studied for interaction with the `Framework`. If a loaded (with molecules) zeolite needs to be simulated, then ion script may be called to insert the molecules in POSCAR format. If it about studying the interaction of molecule at a specific position in the framework, there is a specific module to do that. 

## Features
- Robust POSCAR reader/writer (Direct/Cartesian, Selective dynamics)
- Add ions from the last frame of a PDB (wrap into cell, keep symbols/counts)
- Add NH3 molecule between two Cartesian points with rigid Hs from a TraPPE .def file
- Add H2 molecule between two Cartesian points with rigid geometry from a TraPPE .def file
- Central `VaspWorkflow` class + CLI entry points

## Install

### From source (editable)

Clone the repository and install in editable mode:

```zsh
git clone https://github.com/<your-username>/vasp_init.git
cd vasp_init
pip install -U pip setuptools wheel  # Ensure modern versions (pip>=23, setuptools>=68)
pip install -e .
```

This installs the package and three console scripts:
- `vasp_init_add_ions` — merge ions from PDB into POSCAR
- `vasp_init_add_nh3` — place NH3 molecule between two atoms
- `vasp_init_add_h2` — place H2 molecule between two atoms

### From PyPI (once published)

```zsh
pip install vasp_init
```

## CLI usage

Add ions from PDB:
```zsh
vasp_init_add_ions \
  --poscar path/to/Framework_0_initial.vasp \
  --pdb path/to/raspa_output.pdb \
  --out path/to/POSCAR_with_ions
```
Options: `--model-index`, `--ion-flags TTT|FFF|...`, `--no-wrap`, `--out-coords Direct|Cartesian`.

Add NH3 between two points:
```zsh
vasp_init_add_nh3 \
  --poscar path/to/Framework_0_initial.vasp \
  --def path/to/NH3_TraPPE.def \
  --x1 1.0 --y1 2.0 --z1 3.0 \
  --x2 8.0 --y2 9.0 --z2 10.0 \
  --place midpoint \
  --out path/to/POSCAR_with_NH3
```
Options: `--place midpoint|first|second`, `--flags TTT|FFF|...`, `--no-wrap`, `--out-coords`.

Add H2 between two points:
```zsh
vasp_init_add_h2 \
  --poscar path/to/Framework_0_initial.vasp \
  --def path/to/H2_TraPPE.def \
  --x1 1.0 --y1 2.0 --z1 3.0 \
  --x2 8.0 --y2 9.0 --z2 10.0 \
  --place midpoint \
  --out path/to/POSCAR_with_H2
```
Options: `--place midpoint|first|second`, `--flags TTT|FFF|...`, `--no-wrap`, `--out-coords`.

## Python API

```python
from vasp_init import VaspWorkflow

wf = VaspWorkflow()
# Ions
wf.add_ions_from_pdb(
    poscar_path="Framework_0_initial.vasp",
    pdb_path="raspa_output.pdb",
    out_path="POSCAR_with_ions",
    model_index=-1,
    ion_flags=(True, True, True),
)
# Ammonia
wf.add_ammonia_between(
    poscar_path="Framework_0_initial.vasp",
    def_path="NH3_TraPPE.def",
    out_path="POSCAR_with_NH3",
    x1=1.0,y1=2.0,z1=3.0,
    x2=8.0,y2=9.0,z2=10.0,
    place='midpoint',
)
# Hydrogen
wf.add_hydrogen_between(
    poscar_path="Framework_0_initial.vasp", 
    def_path="H2_TraPPE.def",
    out_path="POSCAR_with_H2",
    x1=1.0,y1=2.0,z1=3.0,
    x2=8.0,y2=9.0,z2=10.0,
    place='midpoint',
)
```

### Direct function usage

```python
from vasp_init import read_poscar, add_ammonia_to_poscar, add_hydrogen_to_poscar

# Read POSCAR
p = read_poscar("POSCAR") 

# Add NH3 at midpoint with custom z-offset
p_nh3 = add_ammonia_to_poscar(p, "NH3_TraPPE.def", 
                              (2.0, 0.0, 0.0), (8.0, 0.0, 0.0),
                              place='midpoint', offset_z=3.0)

# Add H2 molecule
p_h2 = add_hydrogen_to_poscar(p, "H2_TraPPE.def",
                              (2.0, 0.0, 0.0), (8.0, 0.0, 0.0), 
                              place='midpoint', offset_z=3.0)
```
## Units and conventions
- POSCAR symbols line is preserved if present; if missing, counts are updated but symbols remain omitted.
- Coordinates in PDB and .def are treated as Cartesian Å.
- Wrapping places added atoms within the primary cell by default; disable with `--no-wrap`.

## Notes
### Repo layout and packaging

Project layout (PEP 517/518 with `pyproject.toml`):
- `vasp_init/` – package source (geometry/io/molecules/workflow/cli)
- `tests/` – pytest test suite
- `.github/workflows/ci.yml` – GitHub Actions to run tests on push/PR
- `pyproject.toml` – build metadata, entry points, dependencies
- `LICENSE`, `README.md`

GitHub Actions (`.github/workflows/ci.yml`) will run tests on push/PR.

## Running tests locally
```zsh
pip install -e .
pip install pytest
pytest -q
```

