"""Command-line interface for CSV extraction and optional Phase 2 plotting."""

from __future__ import annotations

import argparse
import csv
import glob
import sys
from dataclasses import asdict
from pathlib import Path

from .io import StructureParseError, is_supported_structure_path, parse_structure
from .plotting import plot_phi_psi
from .residues import CSV_COLUMNS, RESIDUE_MODES, GlyPhiPsiRow, extract_gly_phi_psi


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="glyphipsi",
        description="Extract phi/psi angles from local PDB and mmCIF files and write CSV output.",
    )
    parser.add_argument("inputs", nargs="+", help="Local .pdb, .ent, .cif, or .mmcif files. Globs are accepted.")
    parser.add_argument("--out", required=True, help="Output CSV path.")
    parser.add_argument("--plot", help="Optional PNG path for a phi/psi scatter plot.")
    parser.add_argument(
        "--residue-mode",
        choices=RESIDUE_MODES,
        default="gly",
        help="Residue selection mode. Defaults to glycine-only output.",
    )
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

    missing_paths = [path for path in input_paths if not path.is_file()]
    if missing_paths:
        parser.error("Input file does not exist: " + ", ".join(str(path) for path in missing_paths))

    unsupported_paths = [path for path in input_paths if not is_supported_structure_path(path)]
    if unsupported_paths:
        parser.error("Unsupported input file extension: " + ", ".join(str(path) for path in unsupported_paths))

    all_rows: list[GlyPhiPsiRow] = []
    parsed_files = 0
    parse_errors: list[str] = []

    for path in input_paths:
        try:
            structure = parse_structure(path)
        except StructureParseError as exc:
            message = str(exc)
            parse_errors.append(message)
            print(f"warning: {message}", file=sys.stderr)
            if args.strict:
                return 1
            continue

        parsed_files += 1
        all_rows.extend(
            extract_gly_phi_psi(
                structure,
                str(path),
                residue_mode=args.residue_mode,
                model_policy=args.model_policy,
                altloc_policy=args.altloc_policy,
                break_max_c_n_distance=args.break_max_c_n_distance,
            )
        )

    if parsed_files == 0 and parse_errors:
        return 1

    _write_csv(Path(args.out), all_rows)
    if args.plot:
        plot_phi_psi(all_rows, Path(args.plot), residue_mode=args.residue_mode)
    return 0


def _expand_inputs(inputs: list[str]) -> list[Path]:
    paths: list[Path] = []
    for raw_input in inputs:
        matches = glob.glob(raw_input) if _contains_glob(raw_input) else []
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


if __name__ == "__main__":
    raise SystemExit(main())
