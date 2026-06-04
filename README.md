# GlyPhiPsi

Python project for extracting glycine phi/psi backbone dihedral angles from local PDB and mmCIF files.

See [PROJECT_SPEC.md](PROJECT_SPEC.md) for the initial project specification.

## Phase 1 Scope

GlyPhiPsi reads local `.pdb`, `.ent`, `.cif`, and `.mmcif` protein structure files, filters for standard `GLY` residues, calculates backbone phi and psi angles, writes a CSV table, and can optionally save a Phase 2 scatter plot.

The dihedral definitions are:

- `phi = C(i-1), N(i), CA(i), C(i)`
- `psi = N(i), CA(i), C(i), N(i+1)`

Angles are reported in degrees normalized to `-180` to `+180`.

## Usage

```bash
glyphipsi structure.pdb structure.cif --out results/gly_phi_psi.csv
glyphipsi structure.pdb --out results/gly_phi_psi.csv --plot results/gly_ramachandran.png
```

## Plotting

When `--plot` is supplied, GlyPhiPsi saves a scatter plot with `phi_deg` on the x-axis and `psi_deg` on the y-axis. Both axes span `-180` to `+180` degrees and are labeled `phi (degrees)` and `psi (degrees)`.

## Limitations

GlyPhiPsi intentionally does not implement Ramachandran region classification, density maps, contour maps, favored/allowed/outlier labels, or biological interpretation.

Residues are skipped when required atoms are missing, the glycine is first or last in a chain, neighbouring residues are non-standard, a simple C-N continuity check detects a chain break, or required atoms have unresolved alternate conformations. Multiple models use the first model only.
