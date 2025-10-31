from __future__ import annotations

from typing import List, Tuple, Dict, Optional

from .geometry import mat_inv3, mat_vec, vec_mod1
from .io import Poscar


def parse_def_ammonia(def_path: str) -> List[Tuple[str, float, float, float]]:
    """Parse a TraPPE-style ammonia .def file and return [(name, x, y, z)].
    Keeps only atoms whose name starts with 'N_' or 'H_'. Coordinates are in Å.
    """
    atoms: List[Tuple[str, float, float, float]] = []
    with open(def_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    start_idx = None
    for i, line in enumerate(lines):
        if line.strip().lower().startswith('# atomic positions'):
            start_idx = i + 1
            break

    if start_idx is None:
        raise ValueError('Cannot find "# atomic positions" in def file')

    for line in lines[start_idx:]:
        s = line.strip()
        if not s:
            break
        if s.startswith('#'):
            break
        parts = s.split()
        if len(parts) < 5:
            continue
        name = parts[1]
        try:
            x = float(parts[2]); y = float(parts[3]); z = float(parts[4])
        except ValueError:
            continue
        if name.startswith('N_') or name.startswith('H_'):
            atoms.append((name, x, y, z))
    if not atoms:
        raise ValueError('No N/H atoms found in def file')
    return atoms


def add_ammonia_to_poscar(p: Poscar,
                           def_path: str,
                           coord1_cart: Tuple[float, float, float],
                           coord2_cart: Tuple[float, float, float],
                           place: str = 'midpoint',
                           wrap: bool = True,
                           flags: Optional[Tuple[bool, bool, bool]] = None,
                           offset_from_midpoint: float = 0.0,
                           offset_direction: str = '+') -> Poscar:
    """Add an NH3 molecule to a POSCAR, positioning the N atom and adding rigid Hs.

    - def_path: TraPPE .def file describing NH3 geometry (Å), with atoms N_*, H_*.
    - coord1_cart/coord2_cart: two Cartesian coordinates in Å. The N position is
      chosen at midpoint/first/second based on 'place'.
    - wrap: wrap fractional coords to [0,1).
    - flags: selective dynamics flags for added atoms when POSCAR uses them.

    Returns a new Poscar instance with atoms appended and counts updated.
    """
    atoms = parse_def_ammonia(def_path)
    n_atoms = [a for a in atoms if a[0].startswith('N_')]
    if not n_atoms:
        raise ValueError('No N atom found in def file')
    n_atom = n_atoms[0]
    n_ref = (n_atom[1], n_atom[2], n_atom[3])

    included = [a for a in atoms if a[0].startswith('N_') or a[0].startswith('H_')]
    rels: List[Tuple[str, float, float, float]] = []  # (sym, dx, dy, dz) in Å
    for name, x, y, z in included:
        sym = 'N' if name.startswith('N_') else 'H'
        rels.append((sym, x - n_ref[0], y - n_ref[1], z - n_ref[2]))

    if place == 'midpoint':
        # Base midpoint
        mx = (coord1_cart[0] + coord2_cart[0]) / 2.0
        my = (coord1_cart[1] + coord2_cart[1]) / 2.0
        mz = (coord1_cart[2] + coord2_cart[2]) / 2.0
        # Optional offset along the line (coord1 -> coord2)
        if abs(offset_from_midpoint) > 0:
            vx = coord2_cart[0] - coord1_cart[0]
            vy = coord2_cart[1] - coord1_cart[1]
            vz = coord2_cart[2] - coord1_cart[2]
            norm = (vx*vx + vy*vy + vz*vz) ** 0.5
            if norm == 0:
                raise ValueError("Cannot offset from midpoint: the two points are identical (zero-length vector)")
            ux, uy, uz = vx / norm, vy / norm, vz / norm
            sign = 1.0 if (offset_direction == '+' or offset_direction.lower() == 'plus') else -1.0
            cx = mx + sign * offset_from_midpoint * ux
            cy = my + sign * offset_from_midpoint * uy
            cz = mz + sign * offset_from_midpoint * uz
        else:
            cx, cy, cz = mx, my, mz
    elif place == 'first':
        cx, cy, cz = coord1_cart
    elif place == 'second':
        cx, cy, cz = coord2_cart
    else:
        raise ValueError("place must be one of: midpoint, first, second")

    L = p.lattice_cart()
    Linv = mat_inv3(L)

    # Build new fractional coordinates and track counts per species
    add_frac: List[List[float]] = []
    add_symbols: List[str] = []

    for sym, dx, dy, dz in rels:
        x = cx + dx
        y = cy + dy
        z = cz + dz
        f = mat_vec(Linv, [x, y, z])
        if wrap:
            f = vec_mod1(f)
        add_frac.append(f)
        add_symbols.append(sym)

    # Update symbols and counts
    new_symbols = list(p.symbols) if p.symbols else []
    new_counts = list(p.counts)

    # Count additions
    add_count: Dict[str, int] = {}
    for s in add_symbols:
        add_count[s] = add_count.get(s, 0) + 1

    for s, cnt in add_count.items():
        if new_symbols and s in new_symbols:
            idx = new_symbols.index(s)
            new_counts[idx] += cnt
        elif new_symbols:
            new_symbols.append(s)
            new_counts.append(cnt)
        else:
            # No symbols line originally: only adjust counts
            new_counts.append(cnt)

    # Append coords at the end (grouping by symbol order in new_symbols if present)
    new_frac = list(p.frac_coords)
    if new_symbols:
        # Group added atoms by species order
        for sym in new_symbols:
            for i, s in enumerate(add_symbols):
                if s == sym:
                    new_frac.append(add_frac[i])
    else:
        # No symbols: append in input order
        new_frac.extend(add_frac)

    # Flags
    if p.has_selective:
        new_flags = list(p.flags) if p.flags is not None else [(True, True, True)] * len(p.frac_coords)
        fl = flags if flags is not None else (True, True, True)
        for _ in range(len(new_frac) - len(new_flags)):
            new_flags.append(fl)
    else:
        new_flags = None

    from .io import Poscar as _Poscar  # avoid circular type issues at runtime
    out = _Poscar()
    out.comment = p.comment
    out.scale = p.scale
    out.lattice = [row[:] for row in p.lattice]
    out.symbols = new_symbols
    out.counts = new_counts
    out.has_selective = p.has_selective
    out.coord_type = p.coord_type
    out.frac_coords = new_frac
    out.flags = new_flags
    return out
