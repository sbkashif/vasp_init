#!/usr/bin/env python3
"""CLI tool: merge ions from a PDB into a POSCAR.

Usage:
  python -m vasp_init.cli.add_ions \
    --poscar path/to/POSCAR \
    --pdb path/to/ions.pdb \
    --out path/to/POSCAR_with_ions \
    [--model-index -1] [--ion-flags TTT] [--no-wrap]
"""
from __future__ import annotations

import argparse
from typing import Tuple

from vasp_init.io import read_poscar, write_poscar, read_pdb_last_frame, merge_ions_into_poscar


def _parse_tff(s: str) -> Tuple[bool, bool, bool]:
    s = s.strip().upper()
    if len(s) == 3 and all(ch in 'TF' for ch in s):
        return tuple(ch == 'T' for ch in s)  # type: ignore[return-value]
    raise argparse.ArgumentTypeError("--ion-flags must be like TTT, TFT, FFT, etc.")


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge ions from PDB into POSCAR")
    parser.add_argument("--poscar", required=True, help="Path to input POSCAR/CONTCAR")
    parser.add_argument("--pdb", required=True, help="Path to PDB file (ions) — last frame or MODEL used")
    parser.add_argument("--out", required=True, help="Output POSCAR path")
    parser.add_argument("--model-index", type=int, default=-1, help="MODEL index in multi-model PDB (-1 = last)")
    parser.add_argument("--ion-flags", type=_parse_tff, default=None, help="Selective dynamics flags for ions, e.g., TTT or FFT")
    parser.add_argument("--framework-flags", type=_parse_tff, default=None, help="Selective dynamics flags for framework atoms, e.g., FFF or TTT")
    parser.add_argument("--no-wrap", action="store_true", help="Do not wrap ions into [0,1) fractional cell")
    
    args = parser.parse_args()

    p = read_poscar(args.poscar)
    ions = read_pdb_last_frame(args.pdb, model_index=args.model_index)
    merged = merge_ions_into_poscar(p, ions, wrap=not args.no_wrap, ion_flags=args.ion_flags, framework_flags=args.framework_flags)
    write_poscar(merged, args.out)


if __name__ == "__main__":
    main()
