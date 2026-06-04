import csv
from pathlib import Path

from matplotlib import image as mpimg
import pytest

from glyphipsi.cli import main
from glyphipsi.residues import CSV_COLUMNS

from .test_parser_filters import base_atoms, mixed_atoms, write_pdb


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
    assert rows[0]["residue_group"] == "gly"
    assert rows[0]["residue_number"] == "2"
    assert rows[0]["insertion_code"] == ""
    assert float(rows[0]["phi_deg"]) == pytest.approx(float(rows[0]["phi_deg"]))
    assert float(rows[0]["psi_deg"]) == pytest.approx(float(rows[0]["psi_deg"]))


def test_cli_handles_multiple_inputs_and_preserves_source_file(tmp_path):
    first_path = write_pdb(tmp_path / "first.pdb", base_atoms())
    second_path = write_pdb(tmp_path / "second.pdb", base_atoms())
    out_path = tmp_path / "multi.csv"

    exit_code = main([str(first_path), str(second_path), "--out", str(out_path)])

    assert exit_code == 0
    with out_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert [row["source_file"] for row in rows] == [str(first_path), str(second_path)]


def test_cli_general_mode_writes_general_rows_only(tmp_path):
    pdb_path = write_pdb(tmp_path / "mixed.pdb", mixed_atoms())
    out_path = tmp_path / "general.csv"

    exit_code = main([str(pdb_path), "--out", str(out_path), "--residue-mode", "general"])

    assert exit_code == 0
    with out_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert [(row["residue_name"], row["residue_group"]) for row in rows] == [
        ("SER", "general"),
        ("VAL", "general"),
    ]


def test_cli_gly_vs_general_writes_group_labels(tmp_path):
    pdb_path = write_pdb(tmp_path / "mixed.pdb", mixed_atoms())
    out_path = tmp_path / "comparison.csv"

    exit_code = main([str(pdb_path), "--out", str(out_path), "--residue-mode", "gly-vs-general"])

    assert exit_code == 0
    with out_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert [(row["residue_name"], row["residue_group"]) for row in rows] == [
        ("SER", "general"),
        ("GLY", "gly"),
        ("VAL", "general"),
    ]


def test_cli_writes_plot_when_plot_supplied(tmp_path):
    pdb_path = write_pdb(tmp_path / "input.pdb", base_atoms())
    out_path = tmp_path / "gly_phi_psi.csv"
    plot_path = tmp_path / "plots" / "gly_ramachandran.png"

    exit_code = main([str(pdb_path), "--out", str(out_path), "--plot", str(plot_path)])

    assert exit_code == 0
    assert plot_path.is_file()
    assert plot_path.stat().st_size > 0
    assert mpimg.imread(plot_path).size > 0


def test_cli_writes_comparison_plot_when_requested(tmp_path):
    pdb_path = write_pdb(tmp_path / "mixed.pdb", mixed_atoms())
    out_path = tmp_path / "comparison.csv"
    plot_path = tmp_path / "comparison.png"

    exit_code = main(
        [
            str(pdb_path),
            "--out",
            str(out_path),
            "--plot",
            str(plot_path),
            "--residue-mode",
            "gly-vs-general",
        ]
    )

    assert exit_code == 0
    assert plot_path.is_file()
    assert plot_path.stat().st_size > 0
    assert mpimg.imread(plot_path).size > 0


def test_cli_writes_valid_empty_plot_when_no_rows_pass(tmp_path):
    atoms = [atom for atom in base_atoms() if atom["resname"] != "GLY"]
    pdb_path = write_pdb(tmp_path / "no_gly.pdb", atoms)
    out_path = tmp_path / "empty.csv"
    plot_path = tmp_path / "empty.png"

    exit_code = main([str(pdb_path), "--out", str(out_path), "--plot", str(plot_path)])

    assert exit_code == 0
    assert plot_path.is_file()
    assert plot_path.stat().st_size > 0
    assert mpimg.imread(plot_path).size > 0


def test_cli_writes_valid_empty_comparison_plot_when_no_rows_pass(tmp_path):
    atoms = [atom for atom in base_atoms() if atom["resname"] != "GLY"]
    pdb_path = write_pdb(tmp_path / "empty_comparison.pdb", atoms)
    out_path = tmp_path / "empty_comparison.csv"
    plot_path = tmp_path / "empty_comparison.png"

    exit_code = main(
        [
            str(pdb_path),
            "--out",
            str(out_path),
            "--plot",
            str(plot_path),
            "--residue-mode",
            "gly-vs-general",
        ]
    )

    assert exit_code == 0
    assert out_path.read_text(encoding="utf-8").strip() == ",".join(CSV_COLUMNS)
    assert plot_path.is_file()
    assert plot_path.stat().st_size > 0
    assert mpimg.imread(plot_path).size > 0


def test_plotting_does_not_change_csv_output(tmp_path):
    pdb_path = write_pdb(tmp_path / "input.pdb", base_atoms())
    csv_only_path = tmp_path / "csv_only.csv"
    csv_with_plot_path = tmp_path / "csv_with_plot.csv"
    plot_path = tmp_path / "plot.png"

    assert main([str(pdb_path), "--out", str(csv_only_path)]) == 0
    assert main([str(pdb_path), "--out", str(csv_with_plot_path), "--plot", str(plot_path)]) == 0

    assert csv_with_plot_path.read_text(encoding="utf-8") == csv_only_path.read_text(encoding="utf-8")


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
