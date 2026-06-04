# GlyPhiPsi Roadmap

## v0.1 - Glycine Phi/Psi CSV Extractor

Status: complete.

- Local PDB/mmCIF parsing.
- Standard glycine filtering.
- Phi/psi dihedral calculation.
- CSV output.
- Core skip rules and tests.

## v0.2 - Glycine Ramachandran Plotting

Status: complete.

- Glycine-only Ramachandran scatter plot.
- Fixed phi/psi axes from -180 to +180 degrees.
- Empty dataset handling.
- CLI `--plot` support.

## v0.4 - Local Batch Structure Atlas

Status: complete.

- Multiple local input files.
- Glob pattern inputs.
- Aggregated glycine phi/psi CSV output.
- Aggregated glycine scatter plots.
- Optional per-source-file summary CSV.
- Continue-on-error batch behavior unless `--strict` is supplied.

## Possible Future Work

- Add separate skipped-residue reporting with per-residue skip reasons.
- Add deterministic alternate-conformation selection policies.
- Add optional plot styling controls.
- Add larger example datasets and documentation.

The project does not currently aim to perform full Ramachandran validation, density estimation, contour-map generation, or structural quality assessment.
