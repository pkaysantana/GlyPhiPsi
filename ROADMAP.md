# GlyPhiPsi Roadmap

## v0.1 — Glycine Phi/Psi CSV Extractor

Status: complete.

Implemented:

- local PDB/mmCIF parsing;
- standard glycine filtering;
- phi/psi dihedral calculation;
- CSV output;
- skip rules for missing atoms, chain boundaries, first/last residues, non-standard neighbours, and unresolved alternate conformations;
- automated tests.

## v0.1.1 — Scientific Correctness Hardening

Status: complete.

Implemented:

- dihedral angle convention audit;
- additional tests against independent or trusted angle calculations;
- stronger parser/filtering tests;
- README scope clarification.

## v0.2 — Glycine Ramachandran Plotting

Status: complete.

Implemented:

- glycine-only Ramachandran scatter plot;
- fixed phi/psi axes from -180 to +180 degrees;
- correct axis labels;
- empty dataset handling;
- CLI `--plot` support;
- plotting tests.

## v0.3 — Glycine vs General Residue Comparison

Planned.

Goals:

- extract phi/psi angles for glycine and non-glycine standard residues;
- support comparison plotting;
- produce separate or combined plots for:
  - glycine residues;
  - general non-glycine, non-proline residues;
  - optionally proline and pre-proline later;
- avoid validation labels or density claims.

Scientific aim:

Show visually that glycine occupies a broader region of phi/psi space than most other residues, especially in positive-phi regions.

## v0.4 — Batch Structure Processing

Planned.

Goals:

- process many local structures in one command;
- aggregate glycine phi/psi values across a larger dataset;
- preserve source structure metadata;
- add summary counts by file, chain, and residue type;
- handle failed files gracefully.

Scientific aim:

Move from single-structure examples to a larger empirical glycine point cloud.

## v0.5 — Empirical Glycine Density Map

Planned.

Goals:

- estimate glycine phi/psi density from aggregated structures;
- generate heatmaps or contour plots;
- keep density visualisation separate from validation classification.

Scientific aim:

Begin approximating glycine-specific conformational preferences from observed structural data.

## v0.6 — Prototype Glycine-Specific Validator

Planned.

Goals:

- classify glycine residues against empirical regions;
- report possible outliers;
- label residues by source file, chain, and residue number;
- clearly distinguish prototype classification from professional validation tools.

Scientific aim:

Create a learning-focused prototype of residue-type-specific Ramachandran validation.

## Long-Term Ideas

- Add proline-specific and pre-proline-specific handling.
- Compare against published Ramachandran distributions.
- Add interactive plots.
- Add Jupyter examples.
- Add optional PDB download support.
- Package for PyPI.