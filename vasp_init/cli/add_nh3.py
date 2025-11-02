#!/usr/bin/env python3
"""CLI tool: add an NH3 molecule to a POSCAR using a TraPPE .def file.

Usage:
  python -m vasp_init.cli.add_nh3 \
    --poscar path/to/POSCAR \
    --def path/to/NH3.def \
    --x1 2.0 --y1 0.0 --z1 0.0 \
    --x2 8.0 --y2 0.0 --z2 0.0 \
    --place midpoint \
    --out path/to/POSCAR_with_NH3 \
    [--offset-from-midpoint 1.0 --offset-direction +] [--flags TTT] [--no-wrap]
"""
from __future__ import annotations

import argparse
from typing import Tuple

from vasp_init.io import read_poscar, write_poscar
from vasp_init.molecules import add_ammonia_to_poscar


def _parse_tff(s: str) -> Tuple[bool, bool, bool]:
    s = s.strip().upper()
    if len(s) == 3 and all(ch in 'TF' for ch in s):
        return tuple(ch == 'T' for ch in s)  # type: ignore[return-value]
    raise argparse.ArgumentTypeError("--flags must be like TTT, TFT, FFT, etc.")


def main() -> None:
    parser = argparse.ArgumentParser(description="Add NH3 between two Cartesian points in a POSCAR")
    parser.add_argument("--poscar", required=True, help="Path to input POSCAR/CONTCAR")
    parser.add_argument("--def", dest="def_path", required=True, help="Path to TraPPE NH3 .def file")
    parser.add_argument("--x1", type=float, required=True)
    parser.add_argument("--y1", type=float, required=True)
    parser.add_argument("--z1", type=float, required=True)
    parser.add_argument("--x2", type=float, required=True)
    parser.add_argument("--y2", type=float, required=True)
    parser.add_argument("--z2", type=float, required=True)
    parser.add_argument("--place", choices=["midpoint", "first", "second"], default="midpoint")
    parser.add_argument("--offset-from-midpoint", type=float, default=0.0, help="Optional distance (Å) to offset from midpoint along the line 1→2")
    parser.add_argument("--offset-direction", choices=["+", "-", "plus", "minus"], default="+", help="Direction of offset from midpoint")
    parser.add_argument("--flags", type=_parse_tff, default=None, help="Selective dynamics flags for added NH3 atoms, e.g., TTT")
    parser.add_argument("--framework-flags", type=_parse_tff, default=None, help="Selective dynamics flags for framework atoms, e.g., FFF or TTT")
    parser.add_argument("--no-wrap", action="store_true", help="Do not wrap into [0,1) fractional cell")
    parser.add_argument("--out", required=True, help="Output POSCAR path")
    parser.add_argument("--offset-x", type=float, default=0.0, help="Custom offset in x direction (Å) after midpoint/first/second placement")
    parser.add_argument("--offset-y", type=float, default=0.0, help="Custom offset in y direction (Å) after midpoint/first/second placement")
    parser.add_argument("--offset-z", type=float, default=0.0, help="Custom offset in z direction (Å) after midpoint/first/second placement")
    
    args = parser.parse_args()
    p = read_poscar(args.poscar)

    coord1 = (args.x1, args.y1, args.z1)
    coord2 = (args.x2, args.y2, args.z2)
    p2 = add_ammonia_to_poscar(
        p,
        args.def_path,
        coord1,
        coord2,
        place=args.place,
        wrap=not args.no_wrap,
        flags=args.flags,
        framework_flags=args.framework_flags,
        offset_from_midpoint=args.offset_from_midpoint,
        offset_direction=args.offset_direction,
        offset_x=args.offset_x,
        offset_y=args.offset_y,
        offset_z=args.offset_z,
    )

    write_poscar(p2, args.out)


if __name__ == "__main__":
    main()
