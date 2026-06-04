from pathlib import Path

import pytest

from glyphipsi.io import parse_structure
from glyphipsi.residues import extract_gly_phi_psi


def test_tracked_real_pdb_fixture_extracts_glycine_angles():
    pdb_path = Path("data/examples/3IWX.pdb")

    rows = extract_gly_phi_psi(parse_structure(pdb_path), str(pdb_path))

    assert len(rows) == 12
    assert {row.residue_name for row in rows} == {"GLY"}
    assert {row.residue_group for row in rows} == {"gly"}
    assert {row.chain_id for row in rows} == {"A", "B"}
    assert rows[0].chain_id == "A"
    assert rows[0].residue_number == "13"
    assert rows[0].phi_deg == pytest.approx(-55.244632771080006)
    assert rows[0].psi_deg == pytest.approx(-42.31592473720343)
    assert all(-180.0 <= row.phi_deg <= 180.0 for row in rows)
    assert all(-180.0 <= row.psi_deg <= 180.0 for row in rows)
