from __future__ import annotations

from typing import Optional, Tuple

from .io import read_poscar, write_poscar, read_pdb_last_frame, merge_ions_into_poscar, Poscar
from .molecules import add_ammonia_to_poscar, add_hydrogen_to_poscar


class VaspWorkflow:
    """Central workflow facade for adding ions and molecules to POSCARs."""

    def add_ions_from_pdb(self,
                          poscar_path: str,
                          pdb_path: str,
                          out_path: str,
                          model_index: int = -1,
                          ion_flags: Optional[Tuple[bool, bool, bool]] = None,
                          framework_flags: Optional[Tuple[bool, bool, bool]] = None,
                          wrap: bool = True,
                          out_coords: Optional[str] = None) -> str:
        p = read_poscar(poscar_path)
        ions = read_pdb_last_frame(pdb_path, model_index=model_index)
        merged = merge_ions_into_poscar(p, ions, wrap=wrap, ion_flags=ion_flags, framework_flags=framework_flags)
        write_poscar(merged, out_path, out_coord_type=out_coords)
        return out_path

    def add_ammonia_between(self,
                            poscar_path: str,
                            def_path: str,
                            out_path: str,
                            x1: float, y1: float, z1: float,
                            x2: float, y2: float, z2: float,
                            place: str = 'midpoint',
                            flags: Optional[Tuple[bool, bool, bool]] = None,
                            framework_flags: Optional[Tuple[bool, bool, bool]] = None,
                            wrap: bool = True,
                            out_coords: Optional[str] = None,
                            offset_from_midpoint: float = 0.0,
                            offset_direction: str = '+') -> str:
        p = read_poscar(poscar_path)
        p2: Poscar = add_ammonia_to_poscar(
            p, def_path, (x1, y1, z1), (x2, y2, z2),
            place=place, wrap=wrap, flags=flags,
            framework_flags=framework_flags,
            offset_from_midpoint=offset_from_midpoint,
            offset_direction=offset_direction,
        )
        write_poscar(p2, out_path, out_coord_type=out_coords)
        return out_path

    def add_hydrogen_between(self,
                            poscar_path: str,
                            def_path: str,
                            out_path: str,
                            x1: float, y1: float, z1: float,
                            x2: float, y2: float, z2: float,
                            place: str = 'midpoint',
                            flags: Optional[Tuple[bool, bool, bool]] = None,
                            framework_flags: Optional[Tuple[bool, bool, bool]] = None,
                            wrap: bool = True,
                            out_coords: Optional[str] = None,
                            offset_from_midpoint: float = 0.0,
                            offset_direction: str = '+') -> str:
        p = read_poscar(poscar_path)
        p2: Poscar = add_hydrogen_to_poscar(
            p, def_path, (x1, y1, z1), (x2, y2, z2),
            place=place, wrap=wrap, flags=flags,
            framework_flags=framework_flags,
            offset_from_midpoint=offset_from_midpoint,
            offset_direction=offset_direction,
        )
        write_poscar(p2, out_path, out_coord_type=out_coords)
        return out_path
