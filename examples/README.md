# Examples

These example scripts show how to use the library programmatically to:
- Merge ions from the last frame (or a selected MODEL) of a PDB into a POSCAR
- Add a rigid NH3 molecule between two Cartesian points defined in Å
- Run both steps in a single script

All scripts assume you’ve installed the package locally in editable mode and run from the repo root:

```zsh
pip install -e .
```

Then run the scripts with Python. Replace placeholder paths with your files.

## 1) Merge ions into a POSCAR

```zsh
python examples/add_ions.py \
  --poscar path/to/Framework_0_initial.vasp \
  --pdb path/to/raspa_output.pdb \
  --out path/to/POSCAR_with_ions \
  --model-index -1 \
  --ion-flags TTT \
  --framework-flags FFF     # optional: fix framework atoms (F F F = fixed)
# add --no-wrap to avoid wrapping into the primary cell
```

- MODEL index -1 means "last model"; use 0-based index otherwise.
- `--ion-flags` sets selective dynamics flags for ions (e.g., TTT = movable, FFF = fixed)
- `--framework-flags` sets selective dynamics flags for existing framework atoms (e.g., FFF = fixed)
- If either flag option is provided, "Selective dynamics" will be enabled in the output POSCAR

## 2) Add NH3 between two points

```zsh
python examples/add_nh3.py \
  --poscar path/to/POSCAR_or_CONTCAR \
  --def path/to/NH3_TraPPE.def \
  --x1 2.0 --y1 0.0 --z1 0.0 \
  --x2 8.0 --y2 0.0 --z2 0.0 \
  --place midpoint \
  --out path/to/POSCAR_with_NH3 \
  --offset-from-midpoint 1.0 --offset-direction +   # optional \
  --flags TTT \
  --framework-flags FFF                              # optional: fix framework atoms
# add --no-wrap to avoid wrapping into the primary cell
```

- Coordinates are in Cartesian Å.
- `--place` can be `midpoint`, `first`, or `second`.
- `--offset-from-midpoint` shifts the N location along the line from point 1→2.
- `--flags` sets selective dynamics flags for NH3 atoms (e.g., TTT = movable)
- `--framework-flags` sets selective dynamics flags for framework atoms (e.g., FFF = fixed)
- If either flag option is provided, "Selective dynamics" will be enabled in the output POSCAR

## 3) Full workflow: ions then NH3

```zsh
python examples/full_workflow.py \
  --poscar path/to/Framework_0_initial.vasp \
  --pdb path/to/raspa_output.pdb \
  --def path/to/NH3_TraPPE.def \
  --x1 2.0 --y1 0.0 --z1 0.0 \
  --x2 8.0 --y2 0.0 --z2 0.0 \
  --place midpoint \
  --out path/to/POSCAR_final \
  --out-ions path/to/POSCAR_with_ions   # optional intermediate output \
  --ion-flags TTT \
  --framework-flags FFF \
  --flags TTT
# add --no-wrap-ions / --no-wrap-nh3 to avoid wrapping into the primary cell per-step
```

- `--framework-flags` is applied during the ion merge step and persists through the NH3 addition
- `--ion-flags` applies to the ions, `--flags` applies to NH3 atoms

## Notes
- The scripts don't ship sample datasets; use your own POSCAR/PDB/DEF inputs.
- If your input POSCAR lacks a symbols line, the code will still update counts, but symbols remain omitted by design.
- Coordinates from PDB/DEF are treated as Cartesian Å. Wrapping to [0,1) fractional is on by default and can be disabled with `--no-wrap*`.
- **Selective dynamics flags**:
  - Providing `--ion-flags`, `--flags`, or `--framework-flags` automatically enables "Selective dynamics" in the output
  - Use `TTT` for movable atoms (True/True/True), `FFF` for fixed atoms (False/False/False)
  - Common pattern: `--framework-flags FFF` to fix the framework, `--ion-flags TTT` to allow ions to move
