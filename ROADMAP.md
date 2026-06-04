# GlyPhiPsi Roadmap

## v0.1 - Glycine Phi/Psi CSV Extractor

Status: complete.

Implemented:

- local PDB/mmCIF parsing;
- standard glycine filtering;
- phi/psi dihedral calculation;
- CSV output;
- skip rules for missing atoms, chain boundaries, first/last residues, non-standard neighbours, and unresolved alternate conformations;
- automated tests.

## v0.1.1 - Scientific Correctness Hardening

Status: complete.

Implemented:

- dihedral angle convention audit;
- additional tests against independent or trusted angle calculations;
- stronger parser/filtering tests;
- README scope clarification.

## v0.2 - Glycine Ramachandran Plotting

Status: complete.

Implemented:

- glycine-only Ramachandran scatter plot;
- fixed phi/psi axes from -180 to +180 degrees;
- correct axis labels;
- empty dataset handling;
- CLI `--plot` support;
- plotting tests.

## v0.3 - Glycine vs General Residue Comparison

Status: complete.

Implemented:

- `--residue-mode gly`;
- `--residue-mode general`;
- `--residue-mode gly-vs-general`;
- `residue_group` CSV labels;
- comparison plotting that distinguishes glycine and general residues.

The `general` group currently means standard amino acids excluding `GLY` and `PRO`.

## Possible Future Work

- Add optional plot styling controls.
- Add separate skipped-residue reporting.
- Add deterministic alternate-conformation selection policies.
- Handle proline and pre-proline as explicit separate groups.
- Add larger example datasets and documentation.
- Package for PyPI.

The project does not currently aim to perform full Ramachandran validation, density estimation, or structural quality assessment.
