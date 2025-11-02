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
    """Test custom axis-aligned offsets for NH3 placement.
    
    This test places NH3 at the midpoint of two points in x,y plane (z=0),
    then applies a custom offset in the z direction to place it at half the cell height.
    This matches the real-world use case in the example script.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        poscar_p = os.path.join(tmpdir, 'POSCAR')
        def_p = os.path.join(tmpdir, 'NH3.def')
        make_poscar(poscar_p)
        make_nh3_def(def_p)
        p = read_poscar(poscar_p)
        
        # Two points in x,y plane at z=0
        # Midpoint: (5.0, 5.0, 0.0)
        # Custom offset: +5.0 in z to place at half cell height
        p2 = add_ammonia_to_poscar(
            p, def_p,
            (2.0, 3.0, 0.0), (8.0, 7.0, 0.0),
            place='midpoint',
            wrap=True,
            flags=None,
            framework_flags=None,
            offset_from_midpoint=0.0,
            offset_direction='+',
            offset_x=0.0,
            offset_y=0.0,
            offset_z=5.0,  # Place at half cell height
        )
        
        # The N atom should be at midpoint (5.0, 5.0, 0.0) + offset (0.0, 0.0, 5.0) = (5.0, 5.0, 5.0)
        L = p2.lattice
        def cart_from_frac(frac):
            return [sum(frac[j]*L[i][j] for j in range(3)) for i in range(3)]
        nh3_frac = p2.frac_coords[-4:]
        nh3_cart = [cart_from_frac(f) for f in nh3_frac]
        # N atom is first in NH3
        n_cart = nh3_cart[0]
        assert abs(n_cart[0] - 5.0) < 1e-5, f"Expected x=5.0, got {n_cart[0]}"
        assert abs(n_cart[1] - 5.0) < 1e-5, f"Expected y=5.0, got {n_cart[1]}"
        assert abs(n_cart[2] - 5.0) < 1e-5, f"Expected z=5.0, got {n_cart[2]}"
        
        # Verify H atoms maintain rigid geometry relative to N
        # (just check they exist and are offset from N)
        for i in range(1, 4):
            h_cart = nh3_cart[i]
            dist = sum((h_cart[j] - n_cart[j])**2 for j in range(3))**0.5
            assert 0.8 < dist < 1.1, f"H atom {i} distance from N should be ~0.94 Å, got {dist}"
