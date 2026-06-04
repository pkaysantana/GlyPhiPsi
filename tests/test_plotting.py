import pytest
from matplotlib import pyplot as plt

from glyphipsi.plotting import _build_phi_psi_figure
from glyphipsi.residues import GlyPhiPsiRow


def test_plot_uses_required_axes_labels_limits_and_square_figure():
    figure = _build_phi_psi_figure(
        [
            GlyPhiPsiRow(
                source_file="example.pdb",
                model_id="1",
                chain_id="A",
                residue_name="GLY",
                residue_number="2",
                insertion_code="",
                phi_deg=-75.0,
                psi_deg=145.0,
            )
        ]
    )

    try:
        axis = figure.axes[0]
        width, height = figure.get_size_inches()

        assert axis.get_xlabel() == "phi (degrees)"
        assert axis.get_ylabel() == "psi (degrees)"
        assert axis.get_xlim() == pytest.approx((-180.0, 180.0))
        assert axis.get_ylim() == pytest.approx((-180.0, 180.0))
        assert axis.get_aspect() == pytest.approx(1.0)
        assert width == pytest.approx(height)
    finally:
        plt.close(figure)


def test_empty_plot_uses_required_axes_labels_and_limits():
    figure = _build_phi_psi_figure([])

    try:
        axis = figure.axes[0]

        assert axis.get_xlabel() == "phi (degrees)"
        assert axis.get_ylabel() == "psi (degrees)"
        assert axis.get_xlim() == pytest.approx((-180.0, 180.0))
        assert axis.get_ylim() == pytest.approx((-180.0, 180.0))
    finally:
        plt.close(figure)
