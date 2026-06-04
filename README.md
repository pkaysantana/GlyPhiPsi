# GlyPhiPsi

GlyPhiPsi is a Python command-line tool for glycine-specific Ramachandran analysis.

It parses local PDB or mmCIF protein structure files, extracts backbone phi/psi dihedral angles for standard glycine residues, writes a clean CSV table, and can generate a glycine-specific Ramachandran scatter plot. It can also extract non-glycine, non-proline standard residues for comparison when requested.

## Why Glycine?

Most amino acids have side chains that restrict the backbone conformations they can adopt. Glycine is unusual because its side chain is only a hydrogen atom, making it much less sterically restricted than other residues.

This means glycine can occupy regions of phi/psi conformational space that are less accessible to most other amino acids. GlyPhiPsi focuses on extracting and visualising those glycine-specific backbone angles.

## Current Features

- Parses local `.pdb`, `.ent`, `.cif`, and `.mmcif` files.
- Identifies standard `GLY` residues by default.
- Supports optional residue modes: `gly`, `general`, and `gly-vs-general`.
- Calculates backbone:
  - phi = C(i-1), N(i), CA(i), C(i)
  - psi = N(i), CA(i), C(i), N(i+1)
- Skips residues where phi/psi cannot be calculated deterministically.
- Writes a clean CSV output.
- Generates a glycine-only Ramachandran scatter plot.
- Handles empty outputs gracefully.
- Includes automated tests for geometry, parsing, filtering, CLI behaviour, and plotting.

## Non-Goals

GlyPhiPsi does not currently:

- classify residues as favoured, allowed, or outliers;
- perform full Ramachandran validation;
- infer missing atoms or missing residues;
- repair malformed structures;
- generate density maps or contour maps;
- assess global protein structure quality;
- automatically download structures from the PDB.

The current output is coordinate-derived geometry only.

## Installation

```bash
uv sync --extra dev
```

## Example Usage

Extract glycine phi/psi angles to CSV:

```bash
uv run glyphipsi data/examples/3IWX.pdb --out results/gly_phi_psi.csv
```

Generate a glycine-specific Ramachandran plot:

```bash
uv run glyphipsi data/examples/3IWX.pdb \
  --out results/gly_phi_psi.csv \
  --plot results/gly_ramachandran.png
```

Generate a glycine-vs-general comparison:

```bash
uv run glyphipsi data/examples/3IWX.pdb \
  --out results/gly_vs_general.csv \
  --plot results/gly_vs_general.png \
  --residue-mode gly-vs-general
```

## CSV Output

The output CSV contains one row per accepted glycine residue.

Required columns include:

- `source_file`
- `model_id`
- `chain_id`
- `residue_name`
- `residue_group`
- `residue_number`
- `insertion_code`
- `phi_deg`
- `psi_deg`

Example:

```csv
source_file,model_id,chain_id,residue_name,residue_group,residue_number,insertion_code,phi_deg,psi_deg
example.cif,1,A,GLY,gly,42,,-76.214,148.903
```

## Scientific Notes

Phi and psi are backbone dihedral angles.

For residue `i`:

- phi uses atoms C(i-1), N(i), CA(i), C(i)
- psi uses atoms N(i), CA(i), C(i), N(i+1)

The first residue in a chain cannot have phi calculated, and the last residue in a chain cannot have psi calculated. GlyPhiPsi skips residues where required atoms or neighbours are missing.

The `general` group currently includes standard amino acids except `GLY` and `PRO`. Proline is excluded from this comparison mode for now so it can be handled separately in future work.

## Example Result

A first glycine-specific Ramachandran plot from structure `3IWX` shows glycine residues occupying multiple phi/psi regions, including positive-phi conformational space.

This is consistent with glycine's reduced steric restriction, but a single structure is too small to define general glycine conformational preferences.

## Development

Run tests:

```bash
uv run --extra dev python -m pytest
```

## Roadmap

See [ROADMAP.md](ROADMAP.md).

## License

Apache-2.0.
