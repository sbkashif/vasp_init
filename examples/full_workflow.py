#!/usr/bin/env python3
"""Example script: merge ions from PDB and then add NH3 to a POSCAR.

This runs both steps in sequence and writes the final POSCAR. Optionally
writes an intermediate file after ion merge.

Usage:
  python examples/full_workflow.py \
    --poscar path/to/POSCAR \
    --pdb path/to/ions.pdb [--model-index -1] [--ion-flags TTT] [--no-wrap-ions] \
    --def path/to/NH3.def \
    --x1 2.0 --y1 0.0 --z1 0.0 \
    --x2 8.0 --y2 0.0 --z2 0.0 \
    --place midpoint [--offset-from-midpoint 0.0 --offset-direction +] [--flags TTT] [--no-wrap-nh3] \
    --out path/to/POSCAR_final [--out-ions path/to/POSCAR_with_ions]

Requires the package installed (e.g., `pip install -e .`).
"""
from __future__ import annotations

import argparse
from typing import Tuple

from vasp_init.io import read_poscar, write_poscar, read_pdb_last_frame, merge_ions_into_poscar
from vasp_init.molecules import add_ammonia_to_poscar


def _parse_tff(s: str) -> Tuple[bool, bool, bool]:
    s = s.strip().upper()
    if len(s) == 3 and all(ch in 'TF' for ch in s):
        return tuple(ch == 'T' for ch in s)  # type: ignore[return-value]
    raise argparse.ArgumentTypeError("Flags must be like TTT, TFT, FFT, etc.")


ess = argparse.ArgumentParser(description="Merge ions then add NH3 in one go")
# Ions
ess.add_argument("--poscar", required=True, help="Path to input POSCAR/CONTCAR")
ess.add_argument("--pdb", required=True, help="Path to PDB file for ions")
ess.add_argument("--model-index", type=int, default=-1)
ess.add_argument("--ion-flags", type=_parse_tff, default=None)
ess.add_argument("--framework-flags", type=_parse_tff, default=None, help="Selective dynamics flags for framework atoms, e.g., FFF")
ess.add_argument("--no-wrap-ions", action="store_true")
# NH3
ess.add_argument("--def", dest="def_path", required=True, help="Path to TraPPE NH3 .def file")
ess.add_argument("--x1", type=float, required=True)
ess.add_argument("--y1", type=float, required=True)
ess.add_argument("--z1", type=float, required=True)
ess.add_argument("--x2", type=float, required=True)
ess.add_argument("--y2", type=float, required=True)
ess.add_argument("--z2", type=float, required=True)
ess.add_argument("--place", choices=["midpoint", "first", "second"], default="midpoint")
ess.add_argument("--offset-from-midpoint", type=float, default=0.0)
ess.add_argument("--offset-direction", choices=["+", "-", "plus", "minus"], default="+")
ess.add_argument("--flags", type=_parse_tff, default=None, help="Selective dynamics flags for NH3 atoms")
ess.add_argument("--no-wrap-nh3", action="store_true")
# Output
ess.add_argument("--out", required=True, help="Final POSCAR path (after both steps)")
ess.add_argument("--out-ions", default=None, help="Optional path to write POSCAR after ion merge")


def main() -> None:
    args = ess.parse_args()

    # Read initial and merge ions
    p = read_poscar(args.poscar)
    ions = read_pdb_last_frame(args.pdb, model_index=args.model_index)
    p_ions = merge_ions_into_poscar(p, ions, wrap=not args.no_wrap_ions, ion_flags=args.ion_flags, framework_flags=args.framework_flags)

    if args.out_ions:
        write_poscar(p_ions, args.out_ions)

    # Add NH3
    coord1 = (args.x1, args.y1, args.z1)
    coord2 = (args.x2, args.y2, args.z2)
    p_final = add_ammonia_to_poscar(
        p_ions,
        args.def_path,
        coord1,
        coord2,
        place=args.place,
        wrap=not args.no_wrap_nh3,
        flags=args.flags,
        framework_flags=None,  # framework flags already applied in merge step
        offset_from_midpoint=args.offset_from_midpoint,
        offset_direction=args.offset_direction,
    )

    write_poscar(p_final, args.out)


if __name__ == "__main__":
    main()
