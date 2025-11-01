"""Tests for selective dynamics flag handling in merge and add operations.

Assumes package is installed (e.g., `pip install -e .`) and tests are
run from the repository root with `pytest`.
"""

from vasp_init.io import read_poscar, write_poscar, merge_ions_into_poscar, PdbAtom
from vasp_init.molecules import add_ammonia_to_poscar


def make_poscar_no_selective(path):
    """Create a minimal POSCAR without selective dynamics."""
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


def make_poscar_with_selective(path):
    """Create a minimal POSCAR with selective dynamics already enabled."""
    lines = [
        "# test cell with selective",
        "1.0",
        "  10.0 0.0 0.0",
        "  0.0 10.0 0.0",
        "  0.0 0.0 10.0",
        "Si",
        "2",
        "Selective dynamics",
        "Direct",
        " 0.0 0.0 0.0 F F F",
        " 0.25 0.25 0.25 T T T",
    ]
    with open(path, 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")


def make_nh3_def(path):
    """Create a minimal NH3 def file."""
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


# ============================================================================
# Tests for merge_ions_into_poscar
# ============================================================================

def test_merge_ions_no_flags_no_selective(tmp_path):
    """When input has no selective and no flags provided, output has no selective."""
    poscar_p = tmp_path / 'POSCAR'
    make_poscar_no_selective(str(poscar_p))
    
    p = read_poscar(str(poscar_p))
    assert not p.has_selective
    
    ions = [PdbAtom('Na', 1.0, 2.0, 3.0)]
    merged = merge_ions_into_poscar(p, ions, wrap=True, ion_flags=None, framework_flags=None)
    
    assert not merged.has_selective
    assert merged.flags is None
    assert len(merged.frac_coords) == 3  # 2 Si + 1 Na


def test_merge_ions_with_ion_flags_enables_selective(tmp_path):
    """Providing ion_flags enables selective dynamics even if input doesn't have it."""
    poscar_p = tmp_path / 'POSCAR'
    make_poscar_no_selective(str(poscar_p))
    
    p = read_poscar(str(poscar_p))
    assert not p.has_selective
    
    ions = [PdbAtom('Na', 1.0, 2.0, 3.0)]
    merged = merge_ions_into_poscar(p, ions, wrap=True, ion_flags=(False, False, False), framework_flags=None)
    
    assert merged.has_selective
    assert merged.flags is not None
    assert len(merged.flags) == 3
    # Framework atoms should get default TTT
    assert merged.flags[0] == (True, True, True)
    assert merged.flags[1] == (True, True, True)
    # Ion should get FFF
    assert merged.flags[2] == (False, False, False)


def test_merge_ions_with_framework_flags_enables_selective(tmp_path):
    """Providing framework_flags enables selective dynamics."""
    poscar_p = tmp_path / 'POSCAR'
    make_poscar_no_selective(str(poscar_p))
    
    p = read_poscar(str(poscar_p))
    
    ions = [PdbAtom('Na', 1.0, 2.0, 3.0)]
    merged = merge_ions_into_poscar(p, ions, wrap=True, ion_flags=None, framework_flags=(False, False, False))
    
    assert merged.has_selective
    assert merged.flags is not None
    # Framework atoms should get FFF
    assert merged.flags[0] == (False, False, False)
    assert merged.flags[1] == (False, False, False)
    # Ion should get default TTT
    assert merged.flags[2] == (True, True, True)


def test_merge_ions_both_flags_specified(tmp_path):
    """When both ion_flags and framework_flags are specified."""
    poscar_p = tmp_path / 'POSCAR'
    make_poscar_no_selective(str(poscar_p))
    
    p = read_poscar(str(poscar_p))
    
    ions = [PdbAtom('Na', 1.0, 2.0, 3.0), PdbAtom('Na', 4.0, 5.0, 6.0)]
    merged = merge_ions_into_poscar(
        p, ions, wrap=True, 
        ion_flags=(True, False, True),
        framework_flags=(False, False, False)
    )
    
    assert merged.has_selective
    assert merged.flags is not None
    assert len(merged.flags) == 4  # 2 Si + 2 Na
    # Framework atoms: FFF
    assert merged.flags[0] == (False, False, False)
    assert merged.flags[1] == (False, False, False)
    # Ions: TFT
    assert merged.flags[2] == (True, False, True)
    assert merged.flags[3] == (True, False, True)


def test_merge_ions_preserves_existing_selective(tmp_path):
    """When input already has selective dynamics, it's preserved and extended."""
    poscar_p = tmp_path / 'POSCAR'
    make_poscar_with_selective(str(poscar_p))
    
    p = read_poscar(str(poscar_p))
    assert p.has_selective
    assert p.flags[0] == (False, False, False)
    assert p.flags[1] == (True, True, True)
    
    ions = [PdbAtom('Na', 1.0, 2.0, 3.0)]
    merged = merge_ions_into_poscar(p, ions, wrap=True, ion_flags=(True, True, False), framework_flags=None)
    
    assert merged.has_selective
    assert len(merged.flags) == 3
    # Original flags preserved
    assert merged.flags[0] == (False, False, False)
    assert merged.flags[1] == (True, True, True)
    # Ion gets specified flags
    assert merged.flags[2] == (True, True, False)


def test_merge_ions_framework_flags_override_existing(tmp_path):
    """When framework_flags specified with existing selective, it overrides the original flags."""
    poscar_p = tmp_path / 'POSCAR'
    make_poscar_with_selective(str(poscar_p))
    
    p = read_poscar(str(poscar_p))
    
    ions = [PdbAtom('Na', 1.0, 2.0, 3.0)]
    merged = merge_ions_into_poscar(
        p, ions, wrap=True,
        ion_flags=(True, True, True),
        framework_flags=(False, True, False)  # Override original flags
    )
    
    assert merged.has_selective
    # Framework atoms should all have the new flags
    assert merged.flags[0] == (False, True, False)
    assert merged.flags[1] == (False, True, False)
    # Ion gets specified flags
    assert merged.flags[2] == (True, True, True)


# ============================================================================
# Tests for add_ammonia_to_poscar
# ============================================================================

def test_add_nh3_no_flags_no_selective(tmp_path):
    """When input has no selective and no flags provided, output has no selective."""
    poscar_p = tmp_path / 'POSCAR'
    def_p = tmp_path / 'NH3.def'
    make_poscar_no_selective(str(poscar_p))
    make_nh3_def(str(def_p))
    
    p = read_poscar(str(poscar_p))
    assert not p.has_selective
    
    p2 = add_ammonia_to_poscar(
        p, str(def_p),
        (2.0, 0.0, 0.0), (8.0, 0.0, 0.0),
        place='midpoint', wrap=True,
        flags=None, framework_flags=None
    )
    
    assert not p2.has_selective
    assert p2.flags is None
    assert len(p2.frac_coords) == 6  # 2 Si + 1 N + 3 H


def test_add_nh3_with_flags_enables_selective(tmp_path):
    """Providing flags enables selective dynamics for NH3 atoms."""
    poscar_p = tmp_path / 'POSCAR'
    def_p = tmp_path / 'NH3.def'
    make_poscar_no_selective(str(poscar_p))
    make_nh3_def(str(def_p))
    
    p = read_poscar(str(poscar_p))
    
    p2 = add_ammonia_to_poscar(
        p, str(def_p),
        (2.0, 0.0, 0.0), (8.0, 0.0, 0.0),
        place='midpoint', wrap=True,
        flags=(True, False, True),
        framework_flags=None
    )
    
    assert p2.has_selective
    assert p2.flags is not None
    assert len(p2.flags) == 6
    # Framework atoms: default TTT
    assert p2.flags[0] == (True, True, True)
    assert p2.flags[1] == (True, True, True)
    # NH3 atoms: TFT (1 N + 3 H)
    for i in range(2, 6):
        assert p2.flags[i] == (True, False, True)


def test_add_nh3_with_framework_flags_enables_selective(tmp_path):
    """Providing framework_flags enables selective dynamics."""
    poscar_p = tmp_path / 'POSCAR'
    def_p = tmp_path / 'NH3.def'
    make_poscar_no_selective(str(poscar_p))
    make_nh3_def(str(def_p))
    
    p = read_poscar(str(poscar_p))
    
    p2 = add_ammonia_to_poscar(
        p, str(def_p),
        (2.0, 0.0, 0.0), (8.0, 0.0, 0.0),
        place='midpoint', wrap=True,
        flags=None,
        framework_flags=(False, False, False)
    )
    
    assert p2.has_selective
    assert p2.flags is not None
    # Framework atoms: FFF
    assert p2.flags[0] == (False, False, False)
    assert p2.flags[1] == (False, False, False)
    # NH3 atoms: default TTT
    for i in range(2, 6):
        assert p2.flags[i] == (True, True, True)


def test_add_nh3_both_flags_specified(tmp_path):
    """When both flags and framework_flags are specified."""
    poscar_p = tmp_path / 'POSCAR'
    def_p = tmp_path / 'NH3.def'
    make_poscar_no_selective(str(poscar_p))
    make_nh3_def(str(def_p))
    
    p = read_poscar(str(poscar_p))
    
    p2 = add_ammonia_to_poscar(
        p, str(def_p),
        (2.0, 0.0, 0.0), (8.0, 0.0, 0.0),
        place='midpoint', wrap=True,
        flags=(True, True, False),
        framework_flags=(False, False, False)
    )
    
    assert p2.has_selective
    # Framework: FFF
    assert p2.flags[0] == (False, False, False)
    assert p2.flags[1] == (False, False, False)
    # NH3: TTF
    for i in range(2, 6):
        assert p2.flags[i] == (True, True, False)


def test_add_nh3_preserves_existing_selective(tmp_path):
    """When input already has selective dynamics, it's preserved."""
    poscar_p = tmp_path / 'POSCAR'
    def_p = tmp_path / 'NH3.def'
    make_poscar_with_selective(str(poscar_p))
    make_nh3_def(str(def_p))
    
    p = read_poscar(str(poscar_p))
    assert p.has_selective
    
    p2 = add_ammonia_to_poscar(
        p, str(def_p),
        (2.0, 0.0, 0.0), (8.0, 0.0, 0.0),
        place='midpoint', wrap=True,
        flags=(False, True, False),
        framework_flags=None
    )
    
    assert p2.has_selective
    # Original flags preserved
    assert p2.flags[0] == (False, False, False)
    assert p2.flags[1] == (True, True, True)
    # NH3 atoms get specified flags
    for i in range(2, 6):
        assert p2.flags[i] == (False, True, False)


def test_add_nh3_framework_flags_override_existing(tmp_path):
    """When framework_flags specified with existing selective, it overrides."""
    poscar_p = tmp_path / 'POSCAR'
    def_p = tmp_path / 'NH3.def'
    make_poscar_with_selective(str(poscar_p))
    make_nh3_def(str(def_p))
    
    p = read_poscar(str(poscar_p))
    
    p2 = add_ammonia_to_poscar(
        p, str(def_p),
        (2.0, 0.0, 0.0), (8.0, 0.0, 0.0),
        place='midpoint', wrap=True,
        flags=(True, False, True),
        framework_flags=(True, True, False)  # Override
    )
    
    assert p2.has_selective
    # Framework atoms should have new flags
    assert p2.flags[0] == (True, True, False)
    assert p2.flags[1] == (True, True, False)
    # NH3 atoms
    for i in range(2, 6):
        assert p2.flags[i] == (True, False, True)


# ============================================================================
# Integration test: ions then NH3
# ============================================================================

def test_workflow_ions_then_nh3_with_flags(tmp_path):
    """Test full workflow: merge ions with flags, then add NH3 with flags."""
    poscar_p = tmp_path / 'POSCAR'
    def_p = tmp_path / 'NH3.def'
    make_poscar_no_selective(str(poscar_p))
    make_nh3_def(str(def_p))
    
    # Step 1: merge ions
    p = read_poscar(str(poscar_p))
    ions = [PdbAtom('Na', 1.0, 2.0, 3.0)]
    p_ions = merge_ions_into_poscar(
        p, ions, wrap=True,
        ion_flags=(True, True, True),
        framework_flags=(False, False, False)
    )
    
    assert p_ions.has_selective
    assert p_ions.flags[0] == (False, False, False)  # Si
    assert p_ions.flags[1] == (False, False, False)  # Si
    assert p_ions.flags[2] == (True, True, True)     # Na
    
    # Step 2: add NH3 (framework_flags=None to preserve existing)
    p_final = add_ammonia_to_poscar(
        p_ions, str(def_p),
        (2.0, 0.0, 0.0), (8.0, 0.0, 0.0),
        place='midpoint', wrap=True,
        flags=(True, True, False),
        framework_flags=None  # Don't override existing
    )
    
    assert p_final.has_selective
    assert len(p_final.frac_coords) == 7  # 2 Si + 1 Na + 1 N + 3 H
    # Framework (Si) preserved
    assert p_final.flags[0] == (False, False, False)
    assert p_final.flags[1] == (False, False, False)
    # Ion (Na) preserved
    assert p_final.flags[2] == (True, True, True)
    # NH3 atoms get new flags
    for i in range(3, 7):
        assert p_final.flags[i] == (True, True, False)


def test_write_and_read_selective_roundtrip(tmp_path):
    """Test that selective dynamics survives write/read roundtrip."""
    poscar_p = tmp_path / 'POSCAR'
    poscar_out = tmp_path / 'POSCAR_out'
    make_poscar_no_selective(str(poscar_p))
    
    p = read_poscar(str(poscar_p))
    ions = [PdbAtom('Na', 1.0, 2.0, 3.0)]
    merged = merge_ions_into_poscar(
        p, ions, wrap=True,
        ion_flags=(False, True, False),
        framework_flags=(True, False, True)
    )
    
    # Write and read back
    write_poscar(merged, str(poscar_out))
    p_read = read_poscar(str(poscar_out))
    
    assert p_read.has_selective
    assert len(p_read.flags) == 3
    assert p_read.flags[0] == (True, False, True)
    assert p_read.flags[1] == (True, False, True)
    assert p_read.flags[2] == (False, True, False)
