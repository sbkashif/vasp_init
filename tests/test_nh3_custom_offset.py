"""Test custom offset functionality for NH3 placement."""
from vasp_init.io import read_poscar
from vasp_init.molecules import add_ammonia_to_poscar
import tempfile
import os

def make_poscar(path):
    lines = [
        "# test cell",
        "1.0",
        "  10.0 0.0 0.0",
        "  0.0 10.0 0.0",
        "  0.0 0.0 10.0",
        "Si",
        "2",
        "Direct",
        " 0.0 0.0 0.0",
        " 0.25 0.25 0.25",
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

def test_add_nh3_custom_offset():
    with tempfile.TemporaryDirectory() as tmpdir:
        poscar_p = os.path.join(tmpdir, 'POSCAR')
        def_p = os.path.join(tmpdir, 'NH3.def')
        make_poscar(poscar_p)
        make_nh3_def(def_p)
        p = read_poscar(poscar_p)
        # Place at midpoint, then offset +1.0 in x, -2.0 in y, +0.5 in z
        p2 = add_ammonia_to_poscar(
            p, def_p,
            (2.0, 0.0, 0.0), (8.0, 0.0, 0.0),
            place='midpoint',
            wrap=True,
            flags=None,
            framework_flags=None,
            offset_from_midpoint=0.0,
            offset_direction='+',
            offset_x=1.0,
            offset_y=-2.0,
            offset_z=0.5,
        )
        # The N atom should be at (5.0+1.0, 0.0-2.0, 0.0+0.5) = (6.0, -2.0, 0.5)
        # But with wrapping in a 10 Å cell, -2.0 becomes 8.0
        L = p2.lattice
        def cart_from_frac(frac):
            return [sum(frac[j]*L[i][j] for j in range(3)) for i in range(3)]
        nh3_frac = p2.frac_coords[-4:]
        nh3_cart = [cart_from_frac(f) for f in nh3_frac]
        # N atom is first in NH3
        n_cart = nh3_cart[0]
        assert abs(n_cart[0] - 6.0) < 1e-6
        assert abs(n_cart[1] - 8.0) < 1e-6  # wrapped from -2.0 to 8.0
        assert abs(n_cart[2] - 0.5) < 1e-6
        # The other H atoms are rigidly offset from N, so just check N is correct
