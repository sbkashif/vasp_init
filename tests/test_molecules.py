"""Tests for adding NH3 to a POSCAR.

Assumes package is installed (e.g., `pip install -e .`) and tests are
run from the repository root with `pytest`.
"""

from vasp_init.io import read_poscar
from vasp_init.molecules import add_ammonia_to_poscar


def make_min_poscar(path):
    """make_min_poscar
    
    Create a minimal POSCAR file for testing.
    
    Args:
        path (str): Path to write the POSCAR file.
    """
    lines = [
        "# test cell",
        "1.0",
        "  10.0 0.0 0.0",
        "  0.0 10.0 0.0",
        "  0.0 0.0 10.0",
        "Si",
        "1",
        "Direct",
        " 0.000000 0.000000 0.000000",
    ]
    with open(path, 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")


def make_nh3_def(path):
    content = """# ammonia TraPPE
# atomic positions
1 N_1 0.000 0.000 0.000
2 H_1 0.940 0.000 0.000
3 H_2 -0.313 0.889 0.000
4 H_3 -0.313 -0.889 0.000
# end
"""
    with open(path, 'w', encoding='utf-8') as f:
        f.write(content)


def test_add_ammonia_midpoint_and_offset(tmp_path):
    poscar_p = tmp_path / 'POSCAR'
    def_p = tmp_path / 'NH3.def'
    make_min_poscar(str(poscar_p))
    make_nh3_def(str(def_p))

    p = read_poscar(str(poscar_p))

    # place at midpoint between (2,0,0) and (8,0,0)
    p2 = add_ammonia_to_poscar(p, str(def_p), (2.0, 0.0, 0.0), (8.0, 0.0, 0.0), place='midpoint', wrap=True)
    # expect Si + N + 3H counts (but symbols were ['Si'], so counts extended)
    assert len(p2.frac_coords) == len(p.frac_coords) + 4
    # With symbols, we didn't add symbols in a no-symbols situation; here symbols exist ['Si'], so appended
    assert 'N' in p2.symbols and 'H' in p2.symbols

    # offset 1.0 Å toward second point
    p3 = add_ammonia_to_poscar(p, str(def_p), (2.0, 0.0, 0.0), (8.0, 0.0, 0.0), place='midpoint', wrap=True, offset_from_midpoint=1.0, offset_direction='+')
    assert len(p3.frac_coords) == len(p.frac_coords) + 4
