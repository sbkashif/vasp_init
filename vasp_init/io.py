from __future__ import annotations

from typing import List, Tuple, Optional, Dict
from collections import defaultdict

from .geometry import _normalize_element, det3, mat_inv3, mat_vec, vec_mod1


class Poscar:
    def __init__(self):
        self.comment: str = ""
        self.scale: float = 1.0
        self.lattice: List[List[float]] = [[0.0]*3 for _ in range(3)]  # 3x3
        self.symbols: List[str] = []  # optional; may be empty
        self.counts: List[int] = []
        self.has_selective: bool = False
        self.coord_type: str = "Direct"  # or "Cartesian"
        self.frac_coords: List[List[float]] = []  # fractional coords for all atoms
        self.flags: Optional[List[Tuple[bool,bool,bool]]] = None  # if selective

    def lattice_cart(self) -> List[List[float]]:
        # Apply scale; handle negative scale (volume mode)
        L = [[self.lattice[r][c] for c in range(3)] for r in range(3)]
        s = self.scale
        if s > 0:
            for r in range(3):
                for c in range(3):
                    L[r][c] *= s
            return L
        # Negative scale: abs(scale) is target volume
        target_vol = abs(s)
        vol = det3(L)
        if vol <= 0:
            raise ValueError("Invalid lattice volume in POSCAR")
        factor = (target_vol / vol) ** (1.0/3.0)
        for r in range(3):
            for c in range(3):
                L[r][c] *= factor
        return L

    def total_atoms(self) -> int:
        return sum(self.counts)


def _is_int_tokens(tokens: List[str]) -> bool:
    try:
        _ = [int(t) for t in tokens]
        return True
    except Exception:
        return False


def read_poscar(path: str) -> Poscar:
    with open(path, 'r', encoding='utf-8') as f:
        lines = [ln.rstrip('\n') for ln in f]
    if len(lines) < 8:
        raise ValueError("POSCAR seems too short")

    p = Poscar()
    idx = 0
    p.comment = lines[idx].strip(); idx += 1

    p.scale = float(lines[idx].strip()); idx += 1

    # lattice 3x3
    lat = []
    for _ in range(3):
        toks = lines[idx].split()
        if len(toks) < 3:
            raise ValueError("Invalid lattice line in POSCAR")
        lat.append([float(toks[0]), float(toks[1]), float(toks[2])])
        idx += 1
    p.lattice = lat

    # symbols or counts
    toks = lines[idx].split()
    if not toks:
        raise ValueError("Unexpected blank line in POSCAR before counts")

    if _is_int_tokens(toks):
        # no symbols line
        p.symbols = []
        p.counts = [int(t) for t in toks]
        idx += 1
    else:
        p.symbols = [_normalize_element(t) for t in toks]
        idx += 1
        toks2 = lines[idx].split()
        if not _is_int_tokens(toks2):
            raise ValueError("Expected counts line after symbols in POSCAR")
        p.counts = [int(t) for t in toks2]
        if len(p.symbols) and len(p.symbols) != len(p.counts):
            if len(p.symbols) > len(p.counts):
                p.symbols = p.symbols[:len(p.counts)]
            else:
                raise ValueError("Mismatch between symbols and counts lengths")
        idx += 1

    # optional Selective dynamics line
    line = lines[idx].strip()
    if line.lower().startswith('s'):
        p.has_selective = True
        idx += 1
        line = lines[idx].strip()

    # coordinate type
    if line.lower().startswith('d'):
        p.coord_type = 'Direct'
    elif line.lower().startswith('c') or line.lower().startswith('k'):
        p.coord_type = 'Cartesian'
    else:
        raise ValueError(f"Unknown coordinate type line: {line}")
    idx += 1

    # Read coordinates; number = sum(counts)
    nat = p.total_atoms()
    frac_coords: List[List[float]] = []
    flags: List[Tuple[bool,bool,bool]] = []

    L = p.lattice_cart()
    Linv = mat_inv3(L)

    for _ in range(nat):
        if idx >= len(lines):
            raise ValueError("Unexpected end of POSCAR while reading coordinates")
        toks = lines[idx].split()
        idx += 1
        if len(toks) < 3:
            raise ValueError("Coordinate line has fewer than 3 values")
        x, y, z = float(toks[0]), float(toks[1]), float(toks[2])
        if p.coord_type == 'Direct':
            fcoord = [x, y, z]
        else:  # Cartesian
            fcoord = mat_vec(Linv, [x, y, z])
        frac_coords.append(fcoord)
        if p.has_selective and len(toks) >= 6:
            flags.append(tuple(t.upper().startswith('T') for t in toks[3:6]))

    p.frac_coords = frac_coords
    if p.has_selective:
        p.flags = flags if flags else [(True, True, True)] * nat
    else:
        p.flags = None

    return p


def write_poscar(p: Poscar, path: str, out_coord_type: Optional[str] = None) -> None:
    out_coord_type = out_coord_type or p.coord_type
    L = p.lattice_cart()

    with open(path, 'w', encoding='utf-8') as f:
        f.write(f"{p.comment}\n")
        f.write(f"{p.scale:.16f}\n")
        for r in range(3):
            f.write("  " + "  ".join(f"{L[r][c]: .16f}" for c in range(3)) + "\n")
        # symbols + counts
        if p.symbols:
            f.write(" ".join(p.symbols) + " \n")
        f.write(" ".join(str(c) for c in p.counts) + " \n")

        if p.has_selective:
            f.write("Selective dynamics\n")
        f.write(out_coord_type + "\n")

        # Coordinates
        if out_coord_type.lower().startswith('d'):
            coords = p.frac_coords
        else:
            coords = [mat_vec(L, fc) for fc in p.frac_coords]

        for i, c in enumerate(coords):
            if p.has_selective and p.flags is not None and i < len(p.flags):
                fl = p.flags[i]
                f.write(f"{c[0]: .16f}  {c[1]: .16f}  {c[2]: .16f}   {'T' if fl[0] else 'F'}  {'T' if fl[1] else 'F'}  {'T' if fl[2] else 'F'}\n")
            else:
                f.write(f"{c[0]: .16f}  {c[1]: .16f}  {c[2]: .16f}\n")


class PdbAtom:
    def __init__(self, element: str, x: float, y: float, z: float):
        self.element = _normalize_element(element)
        self.xyz = [x, y, z]


def _element_from_pdb_line(line: str) -> str:
    # Columns 77-78 are element symbol (right-justified)
    if len(line) >= 78:
        elem = line[76:78].strip()
        if elem:
            return elem
    # Fallback: derive from atom name (columns 13-16)
    name = line[12:16].strip()
    if not name:
        return "X"
    if len(name) >= 2 and name[0].isalpha() and name[1].islower():
        return name[:2]
    for ch in name:
        if ch.isalpha():
            return ch
    return name[0]


def read_pdb_last_frame(path: str, model_index: int = -1) -> List[PdbAtom]:
    atoms_frames: List[List[PdbAtom]] = []
    current: List[PdbAtom] = []
    in_model = False
    have_model = False

    with open(path, 'r', encoding='utf-8') as f:
        for raw in f:
            line = raw.rstrip('\n')
            rec = line[:6].strip().upper()
            if rec == 'MODEL':
                if in_model:
                    atoms_frames.append(current)
                    current = []
                in_model = True
                have_model = True
                continue
            if rec == 'ENDMDL':
                if in_model:
                    atoms_frames.append(current)
                    current = []
                    in_model = False
                continue
            if rec in ('ATOM', 'HETATM'):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    toks = line.split()
                    coords: List[float] = []
                    for t in toks:
                        try:
                            coords.append(float(t))
                        except Exception:
                            pass
                    if len(coords) < 3:
                        continue
                    x, y, z = coords[:3]
                elem = _element_from_pdb_line(line)
                current.append(PdbAtom(elem, x, y, z))

    if have_model:
        if in_model:
            atoms_frames.append(current)
        if not atoms_frames:
            return []
        idx = model_index if model_index >= 0 else len(atoms_frames) - 1
        if idx < 0 or idx >= len(atoms_frames):
            raise IndexError(f"MODEL index {model_index} out of range (n_models={len(atoms_frames)})")
        return atoms_frames[idx]

    return current


def merge_ions_into_poscar(p: Poscar, ions: List[PdbAtom], wrap: bool = True,
                           ion_flags: Optional[Tuple[bool,bool,bool]] = None,
                           framework_flags: Optional[Tuple[bool,bool,bool]] = None) -> Poscar:
    """Merge ions into a POSCAR.
    
    Args:
        p: Input POSCAR
        ions: List of ion atoms to add
        wrap: Whether to wrap fractional coordinates to [0,1)
        ion_flags: Selective dynamics flags for added ions (enables selective dynamics if provided)
        framework_flags: Selective dynamics flags for existing framework atoms (enables selective dynamics if provided)
    
    Returns:
        New POSCAR with ions merged
    """
    if not ions:
        return p

    L = p.lattice_cart()
    Linv = mat_inv3(L)

    ion_frac: List[List[float]] = []
    ion_symbols: List[str] = []
    for atom in ions:
        f = mat_vec(Linv, atom.xyz)
        if wrap:
            f = vec_mod1(f)
        ion_frac.append(f)
        ion_symbols.append(_normalize_element(atom.element))

    sym_to_coords: Dict[str, List[List[float]]] = defaultdict(list)
    for sym, fc in zip(ion_symbols, ion_frac):
        sym_to_coords[sym].append(fc)

    base_symbols = list(p.symbols) if p.symbols else []
    base_counts = list(p.counts)

    new_symbols = base_symbols[:]
    new_counts = base_counts[:] if base_counts else []

    seen_order: List[str] = []
    for sym in ion_symbols:
        if sym not in seen_order:
            seen_order.append(sym)

    for sym in seen_order:
        coords_list = sym_to_coords[sym]
        if new_symbols and sym in new_symbols:
            idx = new_symbols.index(sym)
            new_counts[idx] += len(coords_list)
        elif new_symbols:
            new_symbols.append(sym)
            new_counts.append(len(coords_list))
        else:
            # No symbols line originally: only counts are tracked
            new_counts.append(len(coords_list))

    new_frac = list(p.frac_coords)
    ordering: List[Tuple[str, List[List[float]]]] = []
    if new_symbols:
        for sym in new_symbols:
            if sym in sym_to_coords:
                ordering.append((sym, sym_to_coords[sym]))
    else:
        for sym in seen_order:
            ordering.append((sym, sym_to_coords[sym]))

    for _, coords_list in ordering:
        new_frac.extend(coords_list)

    # Determine if we need selective dynamics
    enable_selective = p.has_selective or ion_flags is not None or framework_flags is not None
    
    if enable_selective:
        # Build flags for framework atoms
        # If framework_flags specified, it overrides any existing flags
        if framework_flags is not None:
            new_flags = [framework_flags] * len(p.frac_coords)
        elif p.has_selective and p.flags is not None:
            new_flags = list(p.flags)
        else:
            new_flags = [(True, True, True)] * len(p.frac_coords)
        
        # Add flags for ions
        fl = ion_flags if ion_flags is not None else (True, True, True)
        for _ in range(len(new_frac) - len(new_flags)):
            new_flags.append(fl)
    else:
        new_flags = None

    out = Poscar()
    out.comment = p.comment
    out.scale = p.scale
    out.lattice = [row[:] for row in p.lattice]
    out.symbols = new_symbols
    out.counts = new_counts
    out.has_selective = enable_selective
    out.coord_type = p.coord_type
    out.frac_coords = new_frac
    out.flags = new_flags
    return out
