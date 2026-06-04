"""Residue filtering and glycine phi/psi extraction."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from Bio.PDB.Polypeptide import is_aa

from .angles import dihedral_degrees, distance

CSV_COLUMNS = [
    "source_file",
    "model_id",
    "chain_id",
    "residue_name",
    "residue_number",
    "insertion_code",
    "phi_deg",
    "psi_deg",
]

STANDARD_AMINO_ACIDS = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}

BACKBONE_ATOMS = {"N", "CA", "C"}
WATER_NAMES = {"HOH", "H2O", "WAT", "DOD"}


@dataclass(frozen=True)
class GlyPhiPsiRow:
    source_file: str
    model_id: str
    chain_id: str
    residue_name: str
    residue_number: str
    insertion_code: str
    phi_deg: float
    psi_deg: float


class AtomSelectionError(RuntimeError):
    """Raised when a required atom is missing or ambiguous."""


def extract_gly_phi_psi(
    structure,
    source_file: str,
    *,
    model_policy: str = "first",
    altloc_policy: str = "skip",
    break_max_c_n_distance: float | None = 1.8,
) -> list[GlyPhiPsiRow]:
    """Extract accepted glycine phi/psi rows from a parsed structure."""
    if model_policy != "first":
        raise ValueError("Phase 1 only supports model_policy='first'.")
    if altloc_policy != "skip":
        raise ValueError("Phase 1 only supports altloc_policy='skip'.")

    model = _first_model(structure)
    if model is None:
        return []

    rows: list[GlyPhiPsiRow] = []
    model_id = _model_identifier(model)

    for chain in model.get_chains():
        chain_residues = list(_chain_residue_candidates(chain))
        for index, residue in enumerate(chain_residues):
            if not is_standard_glycine(residue):
                continue
            if index == 0 or index == len(chain_residues) - 1:
                continue

            previous_residue = chain_residues[index - 1]
            next_residue = chain_residues[index + 1]
            if not is_standard_amino_acid(previous_residue) or not is_standard_amino_acid(next_residue):
                continue

            try:
                previous_c = _required_atom_coord(previous_residue, "C")
                current_n = _required_atom_coord(residue, "N")
                current_ca = _required_atom_coord(residue, "CA")
                current_c = _required_atom_coord(residue, "C")
                next_n = _required_atom_coord(next_residue, "N")
            except AtomSelectionError:
                continue

            if _has_simple_chain_break(previous_c, current_n, current_c, next_n, break_max_c_n_distance):
                continue

            try:
                phi = dihedral_degrees(previous_c, current_n, current_ca, current_c)
                psi = dihedral_degrees(current_n, current_ca, current_c, next_n)
            except ValueError:
                continue

            residue_number, insertion_code = _residue_identifier(residue)
            rows.append(
                GlyPhiPsiRow(
                    source_file=source_file,
                    model_id=model_id,
                    chain_id=str(chain.id).strip() or "",
                    residue_name=_residue_name(residue),
                    residue_number=residue_number,
                    insertion_code=insertion_code,
                    phi_deg=phi,
                    psi_deg=psi,
                )
            )

    return rows


def is_standard_glycine(residue) -> bool:
    return is_standard_amino_acid(residue) and _residue_name(residue) == "GLY"


def is_standard_amino_acid(residue) -> bool:
    hetfield = residue.id[0].strip()
    return not hetfield and _residue_name(residue) in STANDARD_AMINO_ACIDS


def _chain_residue_candidates(chain) -> Iterable:
    for residue in chain.get_residues():
        if _is_chain_residue_candidate(residue):
            yield residue


def _is_chain_residue_candidate(residue) -> bool:
    resname = _residue_name(residue)
    if resname in WATER_NAMES or residue.id[0].strip().upper().startswith("W"):
        return False
    try:
        if is_aa(residue, standard=False):
            return True
    except Exception:
        pass
    if any(residue.has_id(atom_name) for atom_name in BACKBONE_ATOMS):
        return True
    return not residue.id[0].strip()


def _required_atom_coord(residue, atom_name: str) -> tuple[float, float, float]:
    if residue.is_disordered() == 2:
        raise AtomSelectionError(f"Residue has unresolved alternate conformations: {residue}")
    if not residue.has_id(atom_name):
        raise AtomSelectionError(f"Missing required atom {atom_name}: {residue}")

    atom = residue[atom_name]
    if atom.is_disordered():
        raise AtomSelectionError(f"Atom has alternate conformations: {atom_name}")
    if atom.get_altloc() not in (" ", ""):
        raise AtomSelectionError(f"Atom has alternate location identifier: {atom_name}")

    coord = atom.get_coord()
    return (float(coord[0]), float(coord[1]), float(coord[2]))


def _has_simple_chain_break(
    previous_c: tuple[float, float, float],
    current_n: tuple[float, float, float],
    current_c: tuple[float, float, float],
    next_n: tuple[float, float, float],
    max_c_n_distance: float | None,
) -> bool:
    if max_c_n_distance is None:
        return False
    return distance(previous_c, current_n) > max_c_n_distance or distance(current_c, next_n) > max_c_n_distance


def _first_model(structure):
    for model in structure.get_models():
        return model
    return None


def _model_identifier(model) -> str:
    serial_num = getattr(model, "serial_num", None)
    if serial_num not in (None, 0):
        return str(serial_num)
    if isinstance(model.id, int):
        return str(model.id + 1)
    return str(model.id)


def _residue_identifier(residue) -> tuple[str, str]:
    _, residue_number, insertion_code = residue.id
    insertion_code = "" if insertion_code in (" ", "?") else str(insertion_code).strip()
    return str(residue_number), insertion_code


def _residue_name(residue) -> str:
    return residue.get_resname().strip().upper()
