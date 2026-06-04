# GlyPhiPsi

GlyPhiPsi is a Python command-line tool for glycine-specific Ramachandran analysis.

It parses local PDB or mmCIF protein structure files, extracts backbone phi/psi dihedral angles for standard glycine residues, writes a clean CSV table, and can generate a glycine-specific Ramachandran scatter plot. It also supports local batch processing to build an aggregated glycine phi/psi atlas from many structures.

## Current Features

- Parses local `.pdb`, `.ent`, `.cif`, and `.mmcif` files.
- Accepts one file, many files, or glob patterns.
- Identifies standard `GLY` residues.
- Calculates backbone:
  - phi = C(i-1), N(i), CA(i), C(i)
  - psi = N(i), CA(i), C(i), N(i+1)
- Skips residues where phi/psi cannot be calculated deterministically.
- Writes a clean aggregated CSV output.
- Generates glycine-only Ramachandran scatter plots.
- Can write a per-source-file batch summary CSV.
- Continues past file-level parse errors unless `--strict` is supplied.
- Handles empty outputs gracefully.

## Non-Goals

GlyPhiPsi does not currently:

- classify residues as favored, allowed, or outliers;
- perform full Ramachandran validation;
- infer missing atoms or missing residues;
- repair malformed structures;
- generate density maps or contour maps;
- assess global protein structure quality;
- automatically download structures from the PDB;
- package or publish to PyPI.

The current output is coordinate-derived geometry only.

## Installation

```bash
uv sync --extra dev
```

## Single-Structure Usage

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

## Batch Atlas Usage

Process many local structures with a glob pattern:

```bash
uv run glyphipsi data/atlas/*.cif \
  --out results/atlas/gly_phi_psi_atlas.csv \
  --plot results/atlas/gly_atlas.png
```

Write a per-source-file summary:

```bash
uv run glyphipsi data/atlas/*.cif \
  --out results/atlas/gly_phi_psi_atlas.csv \
  --plot results/atlas/gly_atlas.png \
  --summary results/atlas/summary.csv
```

Use strict mode to stop on the first file-level error:

```bash
uv run glyphipsi data/atlas/*.cif \
  --out results/atlas/gly_phi_psi_atlas.csv \
  --summary results/atlas/summary.csv \
  --strict
```

## CSV Output

The main output CSV contains one row per accepted glycine residue.

Required columns include:

- `source_file`
- `model_id`
- `chain_id`
- `residue_name`
- `residue_number`
- `insertion_code`
- `phi_deg`
- `psi_deg`

Example:

```csv
source_file,model_id,chain_id,residue_name,residue_number,insertion_code,phi_deg,psi_deg
example.cif,1,A,GLY,42,,-76.214,148.903
```

## Summary Output

When `--summary` is supplied, GlyPhiPsi writes one summary row per input source.

Summary columns:

- `source_file`
- `parsed_successfully`
- `accepted_residue_count`
- `skipped_or_error_count`
- `error_message`

For successfully parsed files, `accepted_residue_count` is exact. Per-residue skipped counts are not tracked yet, so `skipped_or_error_count` is blank for successful files. For file-level failures, `parsed_successfully` is `false`, `accepted_residue_count` is `0`, and `skipped_or_error_count` is `1`.

## Scientific Notes

Phi and psi are backbone dihedral angles.

For residue `i`:

- phi uses atoms C(i-1), N(i), CA(i), C(i)
- psi uses atoms N(i), CA(i), C(i), N(i+1)

The first residue in a chain cannot have phi calculated, and the last residue in a chain cannot have psi calculated. GlyPhiPsi skips residues where required atoms or neighbours are missing, where simple peptide continuity checks indicate a chain break, or where required atoms have unresolved alternate conformations.

## Development

Run tests:

```bash
uv run --extra dev python -m pytest
```

## Roadmap

See [ROADMAP.md](ROADMAP.md).

## License

Apache-2.0.
