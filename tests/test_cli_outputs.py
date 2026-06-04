import csv
from pathlib import Path

from matplotlib import image as mpimg
import pytest

from glyphipsi.cli import main
from glyphipsi.residues import CSV_COLUMNS

from .test_parser_filters import base_atoms, write_pdb


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_cli_writes_required_csv_columns_and_one_row(tmp_path):
    pdb_path = write_pdb(tmp_path / "input.pdb", base_atoms())
    out_path = tmp_path / "results" / "gly_phi_psi.csv"

    exit_code = main([str(pdb_path), "--out", str(out_path)])

    assert exit_code == 0
    rows = read_csv(out_path)

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


def test_cli_handles_multiple_inputs_and_preserves_source_file(tmp_path):
    first_path = write_pdb(tmp_path / "first.pdb", base_atoms())
    second_path = write_pdb(tmp_path / "second.pdb", base_atoms())
    out_path = tmp_path / "multi.csv"

    exit_code = main([str(first_path), str(second_path), "--out", str(out_path)])

    assert exit_code == 0
    rows = read_csv(out_path)

    assert [row["source_file"] for row in rows] == [str(first_path), str(second_path)]


def test_cli_expands_glob_inputs_for_batch_atlas(tmp_path):
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()
    first_path = write_pdb(atlas_dir / "a_first.pdb", base_atoms())
    second_path = write_pdb(atlas_dir / "b_second.pdb", base_atoms())
    out_path = tmp_path / "results" / "atlas.csv"

    exit_code = main([str(atlas_dir / "*.pdb"), "--out", str(out_path)])

    assert exit_code == 0
    rows = read_csv(out_path)
    assert [row["source_file"] for row in rows] == [str(first_path), str(second_path)]


def test_cli_writes_summary_for_multiple_valid_inputs(tmp_path):
    first_path = write_pdb(tmp_path / "first.pdb", base_atoms())
    second_path = write_pdb(tmp_path / "second.pdb", base_atoms())
    out_path = tmp_path / "atlas.csv"
    summary_path = tmp_path / "summary.csv"

    exit_code = main([str(first_path), str(second_path), "--out", str(out_path), "--summary", str(summary_path)])

    assert exit_code == 0
    assert read_csv(summary_path) == [
        {
            "source_file": str(first_path),
            "parsed_successfully": "true",
            "accepted_residue_count": "1",
            "skipped_or_error_count": "",
            "error_message": "",
        },
        {
            "source_file": str(second_path),
            "parsed_successfully": "true",
            "accepted_residue_count": "1",
            "skipped_or_error_count": "",
            "error_message": "",
        },
    ]


def test_cli_continues_after_failed_file_and_records_summary(tmp_path):
    valid_path = write_pdb(tmp_path / "valid.pdb", base_atoms())
    bad_path = tmp_path / "bad.cif"
    bad_path.write_text("not a valid mmCIF file\n", encoding="utf-8")
    out_path = tmp_path / "atlas.csv"
    summary_path = tmp_path / "summary.csv"

    exit_code = main([str(valid_path), str(bad_path), "--out", str(out_path), "--summary", str(summary_path)])

    assert exit_code == 0
    rows = read_csv(out_path)
    assert len(rows) == 1
    assert rows[0]["source_file"] == str(valid_path)

    summary_rows = read_csv(summary_path)
    assert summary_rows[0] == {
        "source_file": str(valid_path),
        "parsed_successfully": "true",
        "accepted_residue_count": "1",
        "skipped_or_error_count": "",
        "error_message": "",
    }
    assert summary_rows[1]["source_file"] == str(bad_path)
    assert summary_rows[1]["parsed_successfully"] == "false"
    assert summary_rows[1]["accepted_residue_count"] == "0"
    assert summary_rows[1]["skipped_or_error_count"] == "1"
    assert summary_rows[1]["error_message"]


def test_cli_strict_mode_stops_on_failed_file(tmp_path):
    valid_path = write_pdb(tmp_path / "valid.pdb", base_atoms())
    bad_path = tmp_path / "bad.cif"
    bad_path.write_text("not a valid mmCIF file\n", encoding="utf-8")
    out_path = tmp_path / "atlas.csv"
    summary_path = tmp_path / "summary.csv"

    exit_code = main(
        [str(valid_path), str(bad_path), "--out", str(out_path), "--summary", str(summary_path), "--strict"]
    )

    assert exit_code == 1
    assert not out_path.exists()
    assert [row["parsed_successfully"] for row in read_csv(summary_path)] == ["true", "false"]


def test_cli_aggregate_csv_preserves_row_identifiers(tmp_path):
    first_path = write_pdb(tmp_path / "first.pdb", base_atoms())
    second_path = write_pdb(tmp_path / "second.pdb", base_atoms())
    out_path = tmp_path / "atlas.csv"

    exit_code = main([str(first_path), str(second_path), "--out", str(out_path)])

    assert exit_code == 0
    rows = read_csv(out_path)
    assert [
        (
            row["source_file"],
            row["model_id"],
            row["chain_id"],
            row["residue_name"],
            row["residue_number"],
            row["insertion_code"],
            bool(row["phi_deg"]),
            bool(row["psi_deg"]),
        )
        for row in rows
    ] == [
        (str(first_path), "1", "A", "GLY", "2", "", True, True),
        (str(second_path), "1", "A", "GLY", "2", "", True, True),
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


def test_cli_writes_aggregate_plot_for_multiple_inputs(tmp_path):
    first_path = write_pdb(tmp_path / "first.pdb", base_atoms())
    second_path = write_pdb(tmp_path / "second.pdb", base_atoms())
    out_path = tmp_path / "atlas.csv"
    plot_path = tmp_path / "atlas.png"

    exit_code = main([str(first_path), str(second_path), "--out", str(out_path), "--plot", str(plot_path)])

    assert exit_code == 0
    assert len(read_csv(out_path)) == 2
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


def test_cli_writes_empty_aggregate_outputs_gracefully(tmp_path):
    atoms = [atom for atom in base_atoms() if atom["resname"] != "GLY"]
    first_path = write_pdb(tmp_path / "first_no_gly.pdb", atoms)
    second_path = write_pdb(tmp_path / "second_no_gly.pdb", atoms)
    out_path = tmp_path / "empty_atlas.csv"
    plot_path = tmp_path / "empty_atlas.png"
    summary_path = tmp_path / "empty_summary.csv"

    exit_code = main(
        [
            str(first_path),
            str(second_path),
            "--out",
            str(out_path),
            "--plot",
            str(plot_path),
            "--summary",
            str(summary_path),
        ]
    )

    assert exit_code == 0
    assert out_path.read_text(encoding="utf-8").strip() == ",".join(CSV_COLUMNS)
    assert [row["accepted_residue_count"] for row in read_csv(summary_path)] == ["0", "0"]
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
