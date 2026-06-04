import csv
from pathlib import Path

import pytest

from glyphipsi.cli import main
from glyphipsi.residues import CSV_COLUMNS

from .test_parser_filters import base_atoms, write_pdb


def test_cli_writes_required_csv_columns_and_one_row(tmp_path):
    pdb_path = write_pdb(tmp_path / "input.pdb", base_atoms())
    out_path = tmp_path / "results" / "gly_phi_psi.csv"

    exit_code = main([str(pdb_path), "--out", str(out_path)])

    assert exit_code == 0
    with out_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert rows
    assert list(rows[0].keys()) == CSV_COLUMNS
    assert rows[0]["source_file"] == str(pdb_path)
    assert rows[0]["model_id"] == "1"
    assert rows[0]["chain_id"] == "A"
    assert rows[0]["residue_name"] == "GLY"
    assert rows[0]["residue_number"] == "2"
    assert rows[0]["insertion_code"] == ""
    assert float(rows[0]["phi_deg"]) == pytest.approx(float(rows[0]["phi_deg"]))
    assert float(rows[0]["psi_deg"]) == pytest.approx(float(rows[0]["psi_deg"]))


def test_cli_writes_empty_csv_with_header_when_no_rows_pass(tmp_path):
    atoms = [atom for atom in base_atoms() if atom["resname"] != "GLY"]
    pdb_path = write_pdb(tmp_path / "no_gly.pdb", atoms)
    out_path = tmp_path / "empty.csv"

    exit_code = main([str(pdb_path), "--out", str(out_path)])

    assert exit_code == 0
    assert out_path.read_text(encoding="utf-8").strip() == ",".join(CSV_COLUMNS)


def test_cli_accepts_ent_extension(tmp_path):
    ent_path = write_pdb(tmp_path / "input.ent", base_atoms())
    out_path = tmp_path / "out.csv"

    assert main([str(ent_path), "--out", str(out_path)]) == 0


def test_cli_rejects_missing_input_file(tmp_path):
    with pytest.raises(SystemExit) as exc_info:
        main([str(tmp_path / "missing.pdb"), "--out", str(tmp_path / "out.csv")])

    assert exc_info.value.code == 2


def test_cli_rejects_unsupported_extension(tmp_path):
    unsupported = Path(tmp_path / "input.txt")
    unsupported.write_text("not a structure", encoding="utf-8")

    with pytest.raises(SystemExit) as exc_info:
        main([str(unsupported), "--out", str(tmp_path / "out.csv")])

    assert exc_info.value.code == 2
