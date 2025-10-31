from __future__ import annotations

import argparse

from .workflow import VaspWorkflow
from .io import read_poscar
from .geometry import mat_vec


def main_add_ions():
    ap = argparse.ArgumentParser(description="Add ions from last PDB frame to a POSCAR/CONTCAR")
    ap.add_argument('--poscar', required=True, help='Input POSCAR/CONTCAR path (framework)')
    ap.add_argument('--pdb', required=True, help='Input PDB file with ions (from RASPA)')
    ap.add_argument('--out', required=True, help='Output POSCAR path')
    ap.add_argument('--model-index', type=int, default=-1, help='MODEL index in PDB (default: -1 last)')
    ap.add_argument('--ion-flags', default=None, help="Selective-dynamics flags for ions if POSCAR uses them: e.g., 'TTT' or 'FFF'")
    ap.add_argument('--no-wrap', action='store_true', help='Do not wrap ions into the primary cell (default wraps)')
    ap.add_argument('--out-coords', choices=['Direct', 'Cartesian'], default=None, help='Force output coordinate type (default: same as POSCAR)')

    args = ap.parse_args()

    wf = VaspWorkflow()

    ion_flags = None
    if args.ion_flags is not None:
        s = args.ion_flags.strip().upper()
        if len(s) != 3 or any(ch not in 'TF' for ch in s):
            raise ValueError("--ion-flags must be a 3-char combo of T/F like TTT or FFF")
        ion_flags = (s[0]=='T', s[1]=='T', s[2]=='T')

    wf.add_ions_from_pdb(
        poscar_path=args.poscar,
        pdb_path=args.pdb,
        out_path=args.out,
        model_index=args.model_index,
        ion_flags=ion_flags,
        wrap=(not args.no_wrap),
        out_coords=args.out_coords,
    )
    print(f"Wrote updated POSCAR with ions: {args.out}")


essages = """
usage examples:
  vasp-init-add-ions --poscar POSCAR --pdb ions.pdb --out POSCAR_out
  vasp-init-add-nh3 --poscar POSCAR --def NH3.def --x1 ... --y1 ... --z1 ... --x2 ... --y2 ... --z2 ... --out POSCAR_out
"""

def main_add_nh3():
    ap = argparse.ArgumentParser(description="Add an NH3 molecule between two Cartesian coordinates to a POSCAR")
    ap.add_argument('--poscar', required=True, help='Input POSCAR/CONTCAR path (framework)')
    ap.add_argument('--def', dest='deffile', required=True, help='TraPPE .def file with ammonia geometry')
    ap.add_argument('--out', required=True, help='Output POSCAR path')
    # Either provide explicit Cartesian coordinates OR two atom indices (1-based)
    ap.add_argument('--x1', type=float, required=False)
    ap.add_argument('--y1', type=float, required=False)
    ap.add_argument('--z1', type=float, required=False)
    ap.add_argument('--x2', type=float, required=False)
    ap.add_argument('--y2', type=float, required=False)
    ap.add_argument('--z2', type=float, required=False)
    ap.add_argument('--idx1', type=int, required=False, help='1-based atom index for first point (from POSCAR order)')
    ap.add_argument('--idx2', type=int, required=False, help='1-based atom index for second point (from POSCAR order)')
    ap.add_argument('--place', choices=['midpoint', 'first', 'second'], default='midpoint')
    ap.add_argument('--offset', type=float, default=0.0, help='Distance in Å to move from the midpoint along the line (only for place=midpoint)')
    ap.add_argument('--dir', dest='direction', choices=['+','-','plus','minus'], default='+', help='Direction from midpoint along the line ("+" toward idx2/x2, "-" toward idx1/x1)')
    ap.add_argument('--flags', default=None, help="Selective-dynamics flags for added atoms if POSCAR uses them: e.g., 'TTT' or 'FFF'")
    ap.add_argument('--no-wrap', action='store_true', help='Do not wrap molecule into the primary cell (default wraps)')
    ap.add_argument('--out-coords', choices=['Direct', 'Cartesian'], default=None, help='Force output coordinate type (default: same as POSCAR)')

    args = ap.parse_args()

    wf = VaspWorkflow()

    flags = None
    if args.flags is not None:
        s = args.flags.strip().upper()
        if len(s) != 3 or any(ch not in 'TF' for ch in s):
            raise ValueError("--flags must be a 3-char combo of T/F like TTT or FFF")
        flags = (s[0]=='T', s[1]=='T', s[2]=='T')

    # Determine coordinates: either from indices or explicit xyz
    if args.idx1 is not None and args.idx2 is not None:
        p = read_poscar(args.poscar)
        L = p.lattice_cart()
        n_atoms = len(p.frac_coords)
        for name, idx in [('idx1', args.idx1), ('idx2', args.idx2)]:
            if idx < 1 or idx > n_atoms:
                raise IndexError(f"{name}={idx} out of range (1..{n_atoms})")
        f1 = p.frac_coords[args.idx1 - 1]
        f2 = p.frac_coords[args.idx2 - 1]
        x1, y1, z1 = mat_vec(L, f1)
        x2, y2, z2 = mat_vec(L, f2)
    else:
        missing = [v for v in ['x1','y1','z1','x2','y2','z2'] if getattr(args, v) is None]
        if missing:
            raise ValueError("Provide either --idx1/--idx2 or all of --x1 --y1 --z1 --x2 --y2 --z2")
        x1, y1, z1 = args.x1, args.y1, args.z1
        x2, y2, z2 = args.x2, args.y2, args.z2

    wf.add_ammonia_between(
        poscar_path=args.poscar,
        def_path=args.deffile,
        out_path=args.out,
        x1=x1, y1=y1, z1=z1,
        x2=x2, y2=y2, z2=z2,
        place=args.place,
        flags=flags,
        wrap=(not args.no_wrap),
        out_coords=args.out_coords,
        offset_from_midpoint=args.offset,
        offset_direction=args.direction,
    )
    print(f"Wrote updated POSCAR with NH3: {args.out}")
