"""Command-line interface for CSV extraction and optional Phase 2 plotting."""

from __future__ import annotations

import argparse
import csv
import glob
import statistics
import sys
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path

from .io import StructureParseError, is_supported_structure_path, parse_structure
from .plotting import plot_phi_psi
from .residues import CSV_COLUMNS, GlyPhiPsiRow, extract_gly_phi_psi

SUMMARY_COLUMNS = [
    "source_file",
    "parsed_successfully",
    "accepted_residue_count",
    "skipped_or_error_count",
    "error_message",
]

STATS_COLUMNS = [
    "total_input_files",
    "files_parsed_successfully",
    "files_failed",
    "total_accepted_residues",
    "accepted_glycine_residues",
    "unique_source_files_represented",
    "unique_chains_represented",
    "phi_min",
    "phi_max",
    "psi_min",
    "psi_max",
    "mean_phi",
    "mean_psi",
    "median_phi",
    "median_psi",
]


@dataclass(frozen=True)
class SourceSummaryRow:
    source_file: str
    parsed_successfully: bool
    accepted_residue_count: int
    skipped_or_error_count: str
    error_message: str


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="glyphipsi",
        description="Extract glycine phi/psi angles from local PDB and mmCIF files and write CSV output.",
    )
    parser.add_argument("inputs", nargs="+", help="Local .pdb, .ent, .cif, or .mmcif files. Globs are accepted.")
    parser.add_argument("--out", required=True, help="Output CSV path.")
    parser.add_argument("--plot", help="Optional PNG path for a glycine phi/psi scatter plot.")
    parser.add_argument("--summary", help="Optional per-source-file processing summary CSV path.")
    parser.add_argument("--stats", help="Optional aggregate dataset statistics output path, .csv or text.")
    parser.add_argument(
        "--model-policy",
        choices=["first"],
        default="first",
        help="Model handling policy. Phase 1 supports only 'first'.",
    )
    parser.add_argument(
        "--altloc-policy",
        choices=["skip"],
        default="skip",
        help="Alternate conformation policy. Phase 1 skips unresolved required atom alternates.",
    )
    parser.add_argument(
        "--break-max-c-n-distance",
        type=float,
        default=1.8,
        help="Maximum C-N distance in Angstroms for the simple peptide-continuity check.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail on the first file-level parse error instead of continuing with other inputs.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.break_max_c_n_distance <= 0:
        parser.error("--break-max-c-n-distance must be greater than zero.")

    input_paths = _expand_inputs(args.inputs)
    if not input_paths:
        parser.error("No input files matched.")

    if len(input_paths) == 1:
        if not input_paths[0].is_file():
            parser.error(f"Input file does not exist: {input_paths[0]}")
        if not is_supported_structure_path(input_paths[0]):
            parser.error(f"Unsupported input file extension: {input_paths[0]}")

    all_rows: list[GlyPhiPsiRow] = []
    parsed_files = 0
    file_errors = 0
    summary_rows: list[SourceSummaryRow] = []

    for path in input_paths:
        if not path.is_file():
            file_errors += 1
            message = f"Input file does not exist: {path}"
            summary_rows.append(_error_summary(path, message))
            print(f"warning: {message}", file=sys.stderr)
            if args.strict:
                _write_summary_if_requested(args.summary, summary_rows)
                _write_stats_if_requested(args.stats, input_paths, summary_rows, all_rows)
                return 1
            continue

        if not is_supported_structure_path(path):
            file_errors += 1
            message = f"Unsupported input file extension: {path}"
            summary_rows.append(_error_summary(path, message))
            print(f"warning: {message}", file=sys.stderr)
            if args.strict:
                _write_summary_if_requested(args.summary, summary_rows)
                _write_stats_if_requested(args.stats, input_paths, summary_rows, all_rows)
                return 1
            continue

        try:
            structure = parse_structure(path)
        except StructureParseError as exc:
            message = str(exc)
            file_errors += 1
            summary_rows.append(_error_summary(path, message))
            print(f"warning: {message}", file=sys.stderr)
            if args.strict:
                _write_summary_if_requested(args.summary, summary_rows)
                _write_stats_if_requested(args.stats, input_paths, summary_rows, all_rows)
                return 1
            continue

        parsed_files += 1
        file_rows = extract_gly_phi_psi(
            structure,
            str(path),
            model_policy=args.model_policy,
            altloc_policy=args.altloc_policy,
            break_max_c_n_distance=args.break_max_c_n_distance,
        )
        all_rows.extend(file_rows)
        summary_rows.append(
            SourceSummaryRow(
                source_file=str(path),
                parsed_successfully=True,
                accepted_residue_count=len(file_rows),
                skipped_or_error_count="",
                error_message="",
            )
        )

    if parsed_files == 0 and file_errors:
        _write_summary_if_requested(args.summary, summary_rows)
        _write_stats_if_requested(args.stats, input_paths, summary_rows, all_rows)
        return 1

    _write_csv(Path(args.out), all_rows)
    _write_summary_if_requested(args.summary, summary_rows)
    _write_stats_if_requested(args.stats, input_paths, summary_rows, all_rows)
    if args.plot:
        plot_phi_psi(all_rows, Path(args.plot))
    return 0


def _expand_inputs(inputs: list[str]) -> list[Path]:
    paths: list[Path] = []
    for raw_input in inputs:
        matches = glob.glob(raw_input, recursive=True) if _contains_glob(raw_input) else []
        if matches:
            paths.extend(Path(match) for match in sorted(matches))
        else:
            paths.append(Path(raw_input))
    return paths


def _contains_glob(path: str) -> bool:
    return any(character in path for character in "*?[]")


def _write_csv(path: Path, rows: list[GlyPhiPsiRow]) -> None:
    if path.parent != Path("."):
        path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for row in rows:
            data = asdict(row)
            data["phi_deg"] = f"{row.phi_deg:.6f}"
            data["psi_deg"] = f"{row.psi_deg:.6f}"
            writer.writerow(data)


def _write_summary_if_requested(path: str | None, rows: list[SourceSummaryRow]) -> None:
    if path:
        _write_summary_csv(Path(path), rows)


def _write_summary_csv(path: Path, rows: list[SourceSummaryRow]) -> None:
    if path.parent != Path("."):
        path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS)
        writer.writeheader()
        for row in rows:
            data = asdict(row)
            data["parsed_successfully"] = "true" if row.parsed_successfully else "false"
            writer.writerow(data)


def _write_stats_if_requested(
    path: str | None,
    input_paths: list[Path],
    summary_rows: list[SourceSummaryRow],
    rows: list[GlyPhiPsiRow],
) -> None:
    if path:
        _write_stats(Path(path), _build_stats(input_paths, summary_rows, rows))


def _build_stats(
    input_paths: list[Path],
    summary_rows: list[SourceSummaryRow],
    rows: list[GlyPhiPsiRow],
) -> dict[str, str]:
    phi_values = [row.phi_deg for row in rows]
    psi_values = [row.psi_deg for row in rows]
    stats = {
        "total_input_files": str(len(input_paths)),
        "files_parsed_successfully": str(sum(1 for row in summary_rows if row.parsed_successfully)),
        "files_failed": str(sum(1 for row in summary_rows if not row.parsed_successfully)),
        "total_accepted_residues": str(len(rows)),
        "accepted_glycine_residues": str(sum(1 for row in rows if row.residue_name == "GLY")),
        "unique_source_files_represented": str(len({row.source_file for row in rows})),
        "unique_chains_represented": str(len({(row.source_file, row.model_id, row.chain_id) for row in rows})),
        "phi_min": _format_float(min(phi_values)) if phi_values else "",
        "phi_max": _format_float(max(phi_values)) if phi_values else "",
        "psi_min": _format_float(min(psi_values)) if psi_values else "",
        "psi_max": _format_float(max(psi_values)) if psi_values else "",
        "mean_phi": _format_float(statistics.fmean(phi_values)) if phi_values else "",
        "mean_psi": _format_float(statistics.fmean(psi_values)) if psi_values else "",
        "median_phi": _format_float(statistics.median(phi_values)) if phi_values else "",
        "median_psi": _format_float(statistics.median(psi_values)) if psi_values else "",
    }

    residue_groups = Counter(
        group
        for row in rows
        for group in [getattr(row, "residue_group", None)]
        if group
    )
    for residue_group in sorted(residue_groups):
        stats[f"residue_group_{residue_group}_count"] = str(residue_groups[residue_group])

    return stats


def _write_stats(path: Path, stats: dict[str, str]) -> None:
    if path.parent != Path("."):
        path.parent.mkdir(parents=True, exist_ok=True)

    if path.suffix.lower() == ".csv":
        fieldnames = STATS_COLUMNS + [key for key in stats if key not in STATS_COLUMNS]
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(stats)
        return

    with path.open("w", encoding="utf-8") as handle:
        for key in STATS_COLUMNS:
            handle.write(f"{key}: {stats.get(key, '')}\n")
        for key in stats:
            if key not in STATS_COLUMNS:
                handle.write(f"{key}: {stats[key]}\n")


def _format_float(value: float) -> str:
    return f"{value:.6f}"


def _error_summary(path: Path, message: str) -> SourceSummaryRow:
    return SourceSummaryRow(
        source_file=str(path),
        parsed_successfully=False,
        accepted_residue_count=0,
        skipped_or_error_count="1",
        error_message=message,
    )


if __name__ == "__main__":
    raise SystemExit(main())
