# GlyPhiPsi

Python project for extracting glycine phi/psi backbone dihedral angles from local PDB and mmCIF files.

See [PROJECT_SPEC.md](PROJECT_SPEC.md) for the initial project specification.

## Phase 1 Scope

Phase 1 reads local `.pdb`, `.ent`, `.cif`, and `.mmcif` protein structure files, filters for standard `GLY` residues, calculates backbone phi and psi angles, and writes a CSV table.

The dihedral definitions are:

- `phi = C(i-1), N(i), CA(i), C(i)`
- `psi = N(i), CA(i), C(i), N(i+1)`

Angles are reported in degrees normalized to `-180` to `+180`.

## Usage

```bash
glyphipsi structure.pdb structure.cif --out results/gly_phi_psi.csv
```

## Phase 1 Limitations

Phase 1 writes CSV output only. It intentionally does not implement plotting, Ramachandran region classification, density maps, favored/allowed/outlier labels, or biological interpretation.

Residues are skipped when required atoms are missing, the glycine is first or last in a chain, neighbouring residues are non-standard, a simple C-N continuity check detects a chain break, or required atoms have unresolved alternate conformations. Multiple models use the first model only.
