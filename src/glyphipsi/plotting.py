"""Phi/psi scatter plotting."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt  # noqa: E402

from .residues import GlyPhiPsiRow, ResidueMode


GROUP_STYLES = {
    "gly": {"label": "gly", "color": "#1f77b4", "marker": "o", "size": 22, "alpha": 0.8},
    "general": {"label": "general", "color": "#6b7280", "marker": ".", "size": 14, "alpha": 0.45},
}


def plot_phi_psi(rows: Sequence[GlyPhiPsiRow], path: Path, *, residue_mode: ResidueMode = "gly") -> None:
    """Save a phi/psi scatter plot."""
    if path.parent != Path("."):
        path.parent.mkdir(parents=True, exist_ok=True)

    figure = _build_phi_psi_figure(rows, residue_mode=residue_mode)
    try:
        figure.savefig(path, dpi=150)
    finally:
        plt.close(figure)


def _build_phi_psi_figure(rows: Sequence[GlyPhiPsiRow], *, residue_mode: ResidueMode = "gly"):
    figure, axis = plt.subplots(figsize=(6, 6))

    if residue_mode == "gly-vs-general":
        _scatter_group(axis, rows, "general")
        _scatter_group(axis, rows, "gly")
        if rows:
            axis.legend(frameon=False, loc="upper right")
    else:
        group = "general" if residue_mode == "general" else "gly"
        _scatter_group(axis, rows, group, show_label=False)

    axis.set_xlim(-180, 180)
    axis.set_ylim(-180, 180)
    axis.set_xlabel("phi (degrees)")
    axis.set_ylabel("psi (degrees)")
    axis.set_aspect("equal", adjustable="box")
    axis.grid(True, linewidth=0.5, alpha=0.35)
    figure.tight_layout()
    return figure


def _scatter_group(axis, rows: Sequence[GlyPhiPsiRow], residue_group: str, *, show_label: bool = True) -> None:
    group_rows = [row for row in rows if row.residue_group == residue_group]
    if not group_rows:
        return

    style = GROUP_STYLES[residue_group]
    axis.scatter(
        [row.phi_deg for row in group_rows],
        [row.psi_deg for row in group_rows],
        s=style["size"],
        alpha=style["alpha"],
        color=style["color"],
        marker=style["marker"],
        edgecolors="none",
        label=style["label"] if show_label else None,
    )
