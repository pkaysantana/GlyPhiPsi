"""Phase 2 glycine phi/psi scatter plotting."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt  # noqa: E402

from .residues import GlyPhiPsiRow


def plot_phi_psi(rows: Sequence[GlyPhiPsiRow], path: Path) -> None:
    """Save a glycine phi/psi scatter plot."""
    if path.parent != Path("."):
        path.parent.mkdir(parents=True, exist_ok=True)

    figure = _build_phi_psi_figure(rows)
    try:
        figure.savefig(path, dpi=150)
    finally:
        plt.close(figure)


def _build_phi_psi_figure(rows: Sequence[GlyPhiPsiRow]):
    figure, axis = plt.subplots(figsize=(6, 6))

    phi_values = [row.phi_deg for row in rows]
    psi_values = [row.psi_deg for row in rows]
    if phi_values and psi_values:
        axis.scatter(phi_values, psi_values, s=18, alpha=0.75, color="#1f77b4", edgecolors="none")

    axis.set_xlim(-180, 180)
    axis.set_ylim(-180, 180)
    axis.set_xlabel("phi (degrees)")
    axis.set_ylabel("psi (degrees)")
    axis.set_aspect("equal", adjustable="box")
    axis.grid(True, linewidth=0.5, alpha=0.35)
    figure.tight_layout()
    return figure
