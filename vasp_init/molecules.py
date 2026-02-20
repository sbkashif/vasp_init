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
                           framework_flags: Optional[Tuple[bool, bool, bool]] = None,
                           offset_from_midpoint: float = 0.0,
                           offset_direction: str = '+',
                           offset_x: float = 0.0,
                           offset_y: float = 0.0,
                           offset_z: float = 0.0) -> Poscar:
    """Add an NH3 molecule to a POSCAR, positioning the N atom and adding rigid Hs.

    - def_path: TraPPE .def file describing NH3 geometry (Å), with atoms N_*, H_*.
    - coord1_cart/coord2_cart: two Cartesian coordinates in Å. The N position is
      chosen at midpoint/first/second based on 'place'.
    - wrap: wrap fractional coords to [0,1).
    - flags: selective dynamics flags for added NH3 atoms (enables selective dynamics if provided).
    - framework_flags: selective dynamics flags for existing framework atoms (enables selective dynamics if provided).

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
        # Apply custom offsets along axes
        cx += offset_x
        cy += offset_y
        cz += offset_z
    elif place == 'first':
        cx, cy, cz = coord1_cart
        cx += offset_x
        cy += offset_y
        cz += offset_z
    elif place == 'second':
        cx, cy, cz = coord2_cart
        cx += offset_x
        cy += offset_y
        cz += offset_z
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

    # Determine if we need selective dynamics
    enable_selective = p.has_selective or flags is not None or framework_flags is not None
    
    # Flags
    if enable_selective:
        # Build flags for framework atoms
        # If framework_flags specified, it overrides any existing flags
        if framework_flags is not None:
            new_flags = [framework_flags] * len(p.frac_coords)
        elif p.has_selective and p.flags is not None:
            new_flags = list(p.flags)
        else:
            new_flags = [(True, True, True)] * len(p.frac_coords)
        
        # Add flags for NH3 atoms
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
    out.has_selective = enable_selective
    out.coord_type = p.coord_type
    out.frac_coords = new_frac
    out.flags = new_flags
    return out


def parse_def_hydrogen(def_path: str) -> List[Tuple[str, float, float, float]]:
    """Parse a TraPPE-style hydrogen .def file and return [(name, x, y, z)].
    Keeps only atoms whose name starts with 'H_'. Coordinates are in Å.
    Ignores dummy atoms (M_*).
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
        # Only include H atoms, ignore dummy atoms (M_*)
        if name.startswith('H_'):
            atoms.append((name, x, y, z))
    if not atoms:
        raise ValueError('No H atoms found in def file')
    return atoms


def add_hydrogen_to_poscar(p: Poscar,
                          def_path: str,
                          coord1_cart: Tuple[float, float, float],
                          coord2_cart: Tuple[float, float, float],
                          place: str = 'midpoint',
                          wrap: bool = True,
                          flags: Optional[Tuple[bool, bool, bool]] = None,
                          framework_flags: Optional[Tuple[bool, bool, bool]] = None,
                          offset_from_midpoint: float = 0.0,
                          offset_direction: str = '+',
                          offset_x: float = 0.0,
                          offset_y: float = 0.0,
                          offset_z: float = 0.0) -> Poscar:
    """Add an H2 molecule to a POSCAR, positioning the molecule and adding rigid H atoms.

    - def_path: TraPPE .def file describing H2 geometry (Å), with atoms H_*.
    - coord1_cart/coord2_cart: two Cartesian coordinates in Å. The H2 center position is
      chosen at midpoint/first/second based on 'place'.
    - wrap: wrap fractional coords to [0,1).
    - flags: selective dynamics flags for added H2 atoms (enables selective dynamics if provided).
    - framework_flags: selective dynamics flags for existing framework atoms (enables selective dynamics if provided).

    Returns a new Poscar instance with atoms appended and counts updated.
    """
    atoms = parse_def_hydrogen(def_path)
    if not atoms:
        raise ValueError('No H atoms found in def file')
    
    # Calculate the center of the H2 molecule
    # For H2, this is the midpoint between the two H atoms
    if len(atoms) != 2:
        raise ValueError(f'Expected exactly 2 H atoms in H2 def file, found {len(atoms)}')
    
    h1_pos = (atoms[0][1], atoms[0][2], atoms[0][3])
    h2_pos = (atoms[1][1], atoms[1][2], atoms[1][3])
    molecule_center = ((h1_pos[0] + h2_pos[0]) / 2.0,
                      (h1_pos[1] + h2_pos[1]) / 2.0,
                      (h1_pos[2] + h2_pos[2]) / 2.0)

    included = [a for a in atoms if a[0].startswith('H_')]
    rels: List[Tuple[str, float, float, float]] = []  # (sym, dx, dy, dz) in Å
    for name, x, y, z in included:
        sym = 'H'  # All are H atoms
        rels.append((sym, x - molecule_center[0], y - molecule_center[1], z - molecule_center[2]))

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
        # Apply custom offsets along axes
        cx += offset_x
        cy += offset_y
        cz += offset_z
    elif place == 'first':
        cx, cy, cz = coord1_cart
        cx += offset_x
        cy += offset_y
        cz += offset_z
    elif place == 'second':
        cx, cy, cz = coord2_cart
        cx += offset_x
        cy += offset_y
        cz += offset_z
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

    # Determine if we need selective dynamics
    enable_selective = p.has_selective or flags is not None or framework_flags is not None
    
    # Flags
    if enable_selective:
        # Build flags for framework atoms
        # If framework_flags specified, it overrides any existing flags
        if framework_flags is not None:
            new_flags = [framework_flags] * len(p.frac_coords)
        elif p.has_selective and p.flags is not None:
            new_flags = list(p.flags)
        else:
            new_flags = [(True, True, True)] * len(p.frac_coords)
        
        # Add flags for H2 atoms
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
    out.has_selective = enable_selective
    out.coord_type = p.coord_type
    out.frac_coords = new_frac
    out.flags = new_flags
    return out
