"""Tests for POSCAR IO utilities.

Assumes package is installed (e.g., `pip install -e .`) and tests are
run from the repository root with `pytest`.
"""

from vasp_init.io import read_poscar, write_poscar, merge_ions_into_poscar, PdbAtom


def make_min_poscar(path, coord_type='Direct'):
    lines = [
        "# test cell",
        "1.0",
        "  10.0 0.0 0.0",
        "  0.0 10.0 0.0",
        "  0.0 0.0 10.0",
        "Si",
        "1",
        coord_type,
        " 0.000000 0.000000 0.000000",
    ]
    with open(path, 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")


def test_poscar_read_write_roundtrip(tmp_path):
    pth = tmp_path / 'POSCAR'
    make_min_poscar(str(pth), coord_type='Direct')
    p = read_poscar(str(pth))
    assert p.comment == '# test cell'
    assert p.coord_type == 'Direct'
    assert p.symbols == ['Si']
    assert p.counts == [1]
    # write and re-read
    out = tmp_path / 'POSCAR_out'
    write_poscar(p, str(out))
    p2 = read_poscar(str(out))
    assert p2.counts == [1]
    assert len(p2.frac_coords) == 1


def test_merge_ions_updates_counts_and_coords(tmp_path):
    pth = tmp_path / 'POSCAR'
    make_min_poscar(str(pth))
    p = read_poscar(str(pth))
    # two Na ions in Cartesian coords
    ions = [PdbAtom('Na', 1.0, 2.0, 3.0), PdbAtom('Na', 4.0, 5.0, 6.0)]
    merged = merge_ions_into_poscar(p, ions, wrap=True)
    # symbols should include Na and counts updated
    assert merged.symbols == ['Si', 'Na']
    assert merged.counts == [1, 2]
    # coords increased by 2
    assert len(merged.frac_coords) == 3

