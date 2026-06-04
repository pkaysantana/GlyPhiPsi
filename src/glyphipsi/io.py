"""Structure parsing helpers."""

from __future__ import annotations

from pathlib import Path

from Bio.PDB import MMCIFParser, PDBParser

SUPPORTED_SUFFIXES = {".pdb", ".ent", ".cif", ".mmcif"}


class StructureParseError(RuntimeError):
    """Raised when a structure file cannot be parsed."""


def is_supported_structure_path(path: Path) -> bool:
    return path.suffix.lower() in SUPPORTED_SUFFIXES


def parse_structure(path: Path):
    """Parse a local PDB or mmCIF file with Biopython."""
    suffix = path.suffix.lower()
    if suffix not in SUPPORTED_SUFFIXES:
        raise StructureParseError(f"Unsupported structure file extension: {path}")

    try:
        if suffix in {".pdb", ".ent"}:
            parser = PDBParser(QUIET=True, PERMISSIVE=True)
        else:
            parser = MMCIFParser(QUIET=True)
        return parser.get_structure(path.stem, str(path))
    except Exception as exc:  # pragma: no cover - parser-specific exception types vary.
        raise StructureParseError(f"Failed to parse {path}: {exc}") from exc
