import math
from pathlib import Path

from Bio.PDB.vectors import Vector, calc_dihedral
import pytest

from glyphipsi.angles import normalize_degrees
from glyphipsi.io import parse_structure
from glyphipsi.residues import extract_gly_phi_psi


def test_includes_complete_standard_glycine(tmp_path):
    pdb_path = write_pdb(tmp_path / "complete.pdb", base_atoms())

    rows = extract_rows(pdb_path)

    assert len(rows) == 1
    assert rows[0].residue_name == "GLY"
    assert rows[0].chain_id == "A"
    assert rows[0].residue_number == "2"


def test_extracted_phi_psi_match_biopython_conventional_atom_order(tmp_path):
    pdb_path = write_pdb(tmp_path / "complete.pdb", base_atoms())
    structure = parse_structure(pdb_path)
    rows = extract_gly_phi_psi(structure, str(pdb_path))
    residues = list(next(next(structure.get_models()).get_chains()).get_residues())

    expected_phi = biopython_dihedral(
        parsed_atom_coord(residues[0], "C"),
        parsed_atom_coord(residues[1], "N"),
        parsed_atom_coord(residues[1], "CA"),
        parsed_atom_coord(residues[1], "C"),
    )
    expected_psi = biopython_dihedral(
        parsed_atom_coord(residues[1], "N"),
        parsed_atom_coord(residues[1], "CA"),
        parsed_atom_coord(residues[1], "C"),
        parsed_atom_coord(residues[2], "N"),
    )

    assert len(rows) == 1
    assert rows[0].phi_deg == pytest.approx(expected_phi, abs=1e-12)
    assert rows[0].psi_deg == pytest.approx(expected_psi, abs=1e-12)


@pytest.mark.parametrize("atom_name", ["N", "CA", "C"])
def test_skips_glycine_missing_current_required_atom(tmp_path, atom_name):
    atoms = [atom for atom in base_atoms() if not (atom["resname"] == "GLY" and atom["name"] == atom_name)]
    pdb_path = write_pdb(tmp_path / f"missing_{atom_name}.pdb", atoms)

    assert extract_rows(pdb_path) == []


def test_skips_glycine_missing_previous_c(tmp_path):
    atoms = [atom for atom in base_atoms() if not (atom["resname"] == "ALA" and atom["name"] == "C")]
    pdb_path = write_pdb(tmp_path / "missing_previous_c.pdb", atoms)

    assert extract_rows(pdb_path) == []


def test_skips_glycine_missing_next_n(tmp_path):
    atoms = [atom for atom in base_atoms() if not (atom["resname"] == "SER" and atom["name"] == "N")]
    pdb_path = write_pdb(tmp_path / "missing_next_n.pdb", atoms)

    assert extract_rows(pdb_path) == []


def test_skips_first_chain_residue(tmp_path):
    atoms = [atom for atom in base_atoms() if atom["resseq"] in {2, 3}]
    for atom in atoms:
        if atom["resname"] == "GLY":
            atom["resseq"] = 1
        elif atom["resname"] == "SER":
            atom["resseq"] = 2
    pdb_path = write_pdb(tmp_path / "first_residue.pdb", atoms)

    assert extract_rows(pdb_path) == []


def test_skips_last_chain_residue(tmp_path):
    atoms = [atom for atom in base_atoms() if atom["resseq"] in {1, 2}]
    pdb_path = write_pdb(tmp_path / "last_residue.pdb", atoms)

    assert extract_rows(pdb_path) == []


def test_does_not_calculate_across_chain_boundary(tmp_path):
    atoms = [atom for atom in base_atoms() if atom["resseq"] in {1, 2}]
    atoms.extend(
        [
            atom_dict(7, "N", "SER", "B", 1, 4.5, 1.0, 1.0),
            atom_dict(8, "CA", "SER", "B", 1, 5.0, 2.0, 1.0),
            atom_dict(9, "C", "SER", "B", 1, 6.0, 2.0, 1.0),
        ]
    )
    pdb_path = write_pdb(tmp_path / "chain_boundary.pdb", atoms)

    assert extract_rows(pdb_path) == []


def test_skips_glycine_adjacent_to_simple_chain_break(tmp_path):
    atoms = base_atoms()
    for atom in atoms:
        if atom["resname"] == "ALA" and atom["name"] == "C":
            atom["x"] = -10.0
    pdb_path = write_pdb(tmp_path / "chain_break.pdb", atoms)

    assert extract_rows(pdb_path) == []


def test_skips_non_standard_previous_neighbour(tmp_path):
    atoms = base_atoms()
    for atom in atoms:
        if atom["resname"] == "ALA":
            atom["resname"] = "MSE"
            atom["record"] = "HETATM"
    pdb_path = write_pdb(tmp_path / "nonstandard_neighbour.pdb", atoms)

    assert extract_rows(pdb_path) == []


def test_skips_unresolved_alternate_conformation_on_required_atom(tmp_path):
    atoms = [atom for atom in base_atoms() if not (atom["resname"] == "GLY" and atom["name"] == "CA")]
    atoms.extend(
        [
            atom_dict(10, "CA", "GLY", "A", 2, 2.0, 1.0, 0.0, altloc="A", occupancy=0.5),
            atom_dict(11, "CA", "GLY", "A", 2, 2.1, 1.1, 0.1, altloc="B", occupancy=0.5),
        ]
    )
    pdb_path = write_pdb(tmp_path / "altloc.pdb", atoms)

    assert extract_rows(pdb_path) == []


def test_preserves_insertion_code(tmp_path):
    atoms = base_atoms()
    for atom in atoms:
        if atom["resname"] == "GLY":
            atom["icode"] = "A"
    pdb_path = write_pdb(tmp_path / "insertion_code.pdb", atoms)

    rows = extract_rows(pdb_path)

    assert len(rows) == 1
    assert rows[0].insertion_code == "A"


def test_parses_minimal_mmcif(tmp_path):
    cif_path = write_mmcif(tmp_path / "complete.cif", base_atoms())

    rows = extract_rows(cif_path)

    assert len(rows) == 1
    assert rows[0].residue_name == "GLY"


def extract_rows(path: Path):
    return extract_gly_phi_psi(parse_structure(path), str(path))


def biopython_dihedral(*points):
    return normalize_degrees(math.degrees(calc_dihedral(*(Vector(*point) for point in points))))


def parsed_atom_coord(residue, atom_name):
    coord = residue[atom_name].get_coord()
    return (float(coord[0]), float(coord[1]), float(coord[2]))


def base_atoms():
    return [
        atom_dict(1, "N", "ALA", "A", 1, -1.2, 0.0, 0.0),
        atom_dict(2, "CA", "ALA", "A", 1, -0.6, 0.5, 0.0),
        atom_dict(3, "C", "ALA", "A", 1, 0.0, 0.0, 0.0),
        atom_dict(4, "N", "GLY", "A", 2, 1.3, 0.0, 0.0),
        atom_dict(5, "CA", "GLY", "A", 2, 2.0, 1.0, 0.0),
        atom_dict(6, "C", "GLY", "A", 2, 3.2, 1.0, 1.0),
        atom_dict(7, "N", "SER", "A", 3, 4.5, 1.0, 1.0),
        atom_dict(8, "CA", "SER", "A", 3, 5.0, 2.0, 1.0),
        atom_dict(9, "C", "SER", "A", 3, 6.0, 2.0, 1.0),
    ]


def atom_dict(
    serial,
    name,
    resname,
    chain,
    resseq,
    x,
    y,
    z,
    *,
    altloc=" ",
    occupancy=1.0,
    record="ATOM",
    icode=" ",
):
    return {
        "serial": serial,
        "name": name,
        "resname": resname,
        "chain": chain,
        "resseq": resseq,
        "x": x,
        "y": y,
        "z": z,
        "altloc": altloc,
        "occupancy": occupancy,
        "record": record,
        "icode": icode,
    }


def write_pdb(path: Path, atoms):
    lines = [_pdb_atom_line(**atom) for atom in atoms]
    lines.extend(["TER\n", "END\n"])
    path.write_text("".join(lines), encoding="utf-8")
    return path


def _pdb_atom_line(
    serial,
    name,
    resname,
    chain,
    resseq,
    x,
    y,
    z,
    altloc=" ",
    occupancy=1.0,
    record="ATOM",
    icode=" ",
):
    element = name.strip()[0]
    return (
        f"{record:<6}{serial:5d} {name:>4}{altloc:1}{resname:>3} {chain:1}"
        f"{resseq:4d}{icode:1}   {x:8.3f}{y:8.3f}{z:8.3f}"
        f"{occupancy:6.2f}{20.0:6.2f}          {element:>2}\n"
    )


def write_mmcif(path: Path, atoms):
    header = """data_glyphipsi_test
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
"""
    lines = [header]
    for atom in atoms:
        lines.append(
            "{record} {serial} {element} {name} {altloc} {resname} {chain} 1 {resseq} {icode} "
            "{x:.3f} {y:.3f} {z:.3f} {occupancy:.2f} 20.00 {resseq} {resname} {chain} {name} 1\n".format(
                record=atom["record"],
                serial=atom["serial"],
                element=atom["name"].strip()[0],
                name=atom["name"],
                altloc="." if atom["altloc"] == " " else atom["altloc"],
                resname=atom["resname"],
                chain=atom["chain"],
                resseq=atom["resseq"],
                icode="?" if atom["icode"] == " " else atom["icode"],
                x=atom["x"],
                y=atom["y"],
                z=atom["z"],
                occupancy=atom["occupancy"],
            )
        )
    lines.append("#\n")
    path.write_text("".join(lines), encoding="utf-8")
    return path
