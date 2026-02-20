"""Test hydrogen molecule addition functionality."""
from vasp_init.io import read_poscar
from vasp_init.molecules import add_hydrogen_to_poscar
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

def make_h2_def(path):
    content = """# hydrogen TraPPE
# atomic positions
0 H_h2    0.0000      0.0000    0.0000
1 H_h2    0.7410      0.0000    0.0000
2 M_h2    0.3705     0.0000    0.0000
# end
"""
    with open(path, 'w', encoding='utf-8') as f:
        f.write(content)

def test_add_h2_basic():
    """Test basic H2 molecule addition at midpoint."""
    with tempfile.TemporaryDirectory() as tmpdir:
        poscar_p = os.path.join(tmpdir, 'POSCAR')
        def_p = os.path.join(tmpdir, 'H2.def')
        make_poscar(poscar_p)
        make_h2_def(def_p)
        p = read_poscar(poscar_p)
        
        # Add H2 at midpoint between two points
        p2 = add_hydrogen_to_poscar(
            p, def_p,
            (2.0, 3.0, 4.0), (8.0, 7.0, 6.0),
            place='midpoint',
            wrap=True,
        )
        
        # Check that we added 2 H atoms
        h_count = 0 if 'H' not in p.symbols else p.counts[p.symbols.index('H')]
        h_count_new = p2.counts[p2.symbols.index('H')]
        assert h_count_new == h_count + 2, f"Expected 2 more H atoms, got {h_count_new - h_count}"
        
        # Check that the H2 molecule is placed correctly
        L = p2.lattice
        def cart_from_frac(frac):
            return [sum(frac[j]*L[i][j] for j in range(3)) for i in range(3)]
        
        # Get the H atoms (last 2 atoms added)
        h2_frac = p2.frac_coords[-2:]
        h2_cart = [cart_from_frac(f) for f in h2_frac]
        
        # Check H-H distance is approximately 0.741 Å
        h1_cart, h2_cart_pos = h2_cart[0], h2_cart[1]
        dist = sum((h2_cart_pos[i] - h1_cart[i])**2 for i in range(3))**0.5
        assert abs(dist - 0.741) < 1e-3, f"H-H distance should be 0.741 Å, got {dist}"
        
        # Check that center is approximately at midpoint
        center_x = (h1_cart[0] + h2_cart_pos[0]) / 2.0
        center_y = (h1_cart[1] + h2_cart_pos[1]) / 2.0
        center_z = (h1_cart[2] + h2_cart_pos[2]) / 2.0
        expected_midpoint = (5.0, 5.0, 5.0)  # midpoint of (2,3,4) and (8,7,6)
        
        assert abs(center_x - expected_midpoint[0]) < 1e-3, f"Expected center x={expected_midpoint[0]}, got {center_x}"
        assert abs(center_y - expected_midpoint[1]) < 1e-3, f"Expected center y={expected_midpoint[1]}, got {center_y}"
        assert abs(center_z - expected_midpoint[2]) < 1e-3, f"Expected center z={expected_midpoint[2]}, got {center_z}"

def test_add_h2_custom_offset():
    """Test H2 addition with custom offsets."""
    with tempfile.TemporaryDirectory() as tmpdir:
        poscar_p = os.path.join(tmpdir, 'POSCAR')
        def_p = os.path.join(tmpdir, 'H2.def')
        make_poscar(poscar_p)
        make_h2_def(def_p)
        p = read_poscar(poscar_p)
        
        # Add H2 with custom z offset
        p2 = add_hydrogen_to_poscar(
            p, def_p,
            (0.0, 0.0, 0.0), (10.0, 10.0, 0.0),
            place='midpoint',
            wrap=True,
            offset_x=0.0,
            offset_y=0.0,
            offset_z=3.0,  # 3 Å offset in z
        )
        
        # Get H2 center position
        L = p2.lattice
        def cart_from_frac(frac):
            return [sum(frac[j]*L[i][j] for j in range(3)) for i in range(3)]
        
        h2_frac = p2.frac_coords[-2:]
        h2_cart = [cart_from_frac(f) for f in h2_frac]
        
        # Calculate center
        center_x = (h2_cart[0][0] + h2_cart[1][0]) / 2.0
        center_y = (h2_cart[0][1] + h2_cart[1][1]) / 2.0  
        center_z = (h2_cart[0][2] + h2_cart[1][2]) / 2.0
        
        # Should be at midpoint (5,5,0) + offset (0,0,3) = (5,5,3)
        assert abs(center_x - 5.0) < 1e-3, f"Expected center x=5.0, got {center_x}"
        assert abs(center_y - 5.0) < 1e-3, f"Expected center y=5.0, got {center_y}"
        assert abs(center_z - 3.0) < 1e-3, f"Expected center z=3.0, got {center_z}"