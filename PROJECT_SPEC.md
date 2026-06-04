# GlyPhiPsi Project Specification

## Purpose

GlyPhiPsi is a Python command-line project for glycine-specific Ramachandran analysis. It parses local protein structure files in PDB or mmCIF format, extracts backbone phi and psi dihedral angles for standard glycine residues, writes a clean CSV table, and generates a scatter plot with both axes spanning -180 to +180 degrees.

The initial version intentionally avoids advanced structural validation, statistical classification, favored/outlier labeling, density estimation, and scientific interpretation beyond reporting calculated backbone angles.

## Core Requirements

1. Accept one or more local PDB or mmCIF files.
2. Parse protein chains and residue order from each input file.
3. Identify standard glycine residues only.
4. Calculate:
   - phi: C(i-1), N(i), CA(i), C(i)
   - psi: N(i), CA(i), C(i), N(i+1)
5. Skip residues when phi or psi cannot be calculated deterministically.
6. Write one CSV row per accepted glycine residue.
7. Generate a scatter plot of psi versus phi over the full -180 to +180 degree range.
8. Preserve enough identifiers in the CSV for each plotted point to be traced back to the source file, model, chain, and residue.

## Non-Goals For Initial Version

- No advanced structure validation.
- No Ramachandran region classification.
- No glycine-specific contour maps, density maps, or favored/disallowed labels.
- No biological claims about residue quality, refinement quality, folding correctness, or structural reliability.
- No correction of malformed structures.
- No inference of missing atoms or missing residues.
- No support for trajectories or multi-frame molecular dynamics files.
- No automatic download from the PDB or other remote services.

## Expected Inputs

### Supported File Types

- `.pdb`
- `.ent`
- `.cif`
- `.mmcif`

All inputs are local filesystem paths.

### Command-Line Interface

Proposed command:

```bash
glyphipsi input1.pdb input2.cif --out results/gly_phi_psi.csv --plot results/gly_ramachandran.png
```

Optional proposed arguments:

```bash
glyphipsi structures/*.cif \
  --out results/gly_phi_psi.csv \
  --plot results/gly_ramachandran.png \
  --model-policy first \
  --altloc-policy highest-occupancy \
  --break-max-c-n-distance 1.8
```

### Input Assumptions

- Files contain atomic coordinates for protein structures.
- Residues are processed in chain order as provided by the parser.
- Standard glycine is identified by residue name `GLY`.
- Protein chain boundaries are respected; phi and psi are never calculated across different chains.
- The initial version does not attempt to repair residue numbering, insertion codes, missing residues, or atom naming irregularities.

## Expected Outputs

### CSV Output

The CSV must contain a header row and one row per accepted glycine residue.

Required columns:

| Column | Description |
| --- | --- |
| `source_file` | Input file path or basename, consistently documented. |
| `model_id` | Model identifier used for calculation. |
| `chain_id` | Chain identifier. |
| `residue_name` | Residue name, expected to be `GLY` for accepted rows. |
| `residue_number` | Author or parser residue sequence number. |
| `insertion_code` | Insertion code if present, otherwise empty. |
| `phi_deg` | Phi angle in degrees, normalized to -180 to +180. |
| `psi_deg` | Psi angle in degrees, normalized to -180 to +180. |

Optional useful columns:

| Column | Description |
| --- | --- |
| `altloc_policy` | Alternate conformation policy used. |
| `prev_residue_number` | Previous residue used for phi. |
| `next_residue_number` | Next residue used for psi. |
| `skip_reason` | Only if a separate skipped-residue report is implemented. Do not mix skipped rows into the primary clean CSV unless explicitly requested. |

Example CSV:

```csv
source_file,model_id,chain_id,residue_name,residue_number,insertion_code,phi_deg,psi_deg
example.cif,1,A,GLY,42,,-76.214,148.903
```

### Plot Output

Required plot behavior:

- Scatter plot with phi on the x-axis and psi on the y-axis.
- X-axis range: -180 to +180 degrees.
- Y-axis range: -180 to +180 degrees.
- Axis labels:
  - `phi (degrees)`
  - `psi (degrees)`
- Equal axis scaling or visually square plotting area is preferred.
- Grid lines are allowed if they do not obscure points.
- Output format should support at least `.png`.

Optional plot behavior:

- `.svg` and `.pdf` export.
- Point alpha for dense datasets.
- Color by source file or chain, only when requested.

## Scientific Assumptions

- Backbone dihedral definitions follow the conventional atom order:
  - phi = dihedral angle C(i-1), N(i), CA(i), C(i)
  - psi = dihedral angle N(i), CA(i), C(i), N(i+1)
- Angles are reported in degrees and normalized into the closed-open or closed interval documented by the implementation. The recommended display and CSV convention is -180 to +180 degrees.
- Only standard glycine residues named `GLY` are included.
- The calculation reports coordinate-derived geometry only. It does not determine whether a conformation is favored, allowed, unusual, erroneous, or biologically meaningful.
- Chain continuity is a prerequisite for calculation. A glycine residue requires a valid previous residue in the same chain for phi and a valid next residue in the same chain for psi.
- Hydrogen atoms are not needed for phi or psi calculation.
- Multiple models are handled by a documented policy. The recommended initial policy is `first`, meaning only the first model is processed unless the user requests otherwise.

## Residue Inclusion Rules

A glycine residue is included only when all conditions below are satisfied:

1. Residue name is exactly standard `GLY`.
2. Residue belongs to a protein chain parsed from the structure file.
3. The selected model contains the residue.
4. Previous and next residues exist in the same chain.
5. Previous, current, and next residues are standard amino acid residues.
6. Required atoms are present after alternate-location handling:
   - Previous residue: `C`
   - Current glycine: `N`, `CA`, `C`
   - Next residue: `N`
7. The residue is not at the first or last calculable position of the chain.
8. The residue is not adjacent to a detected chain break.
9. Alternate conformations are absent or resolved by a deterministic policy.

## Required Skip Rules

The parser must skip and not calculate angles for:

- Glycine residues with any missing required atom.
- First residue in a chain, because phi cannot be calculated.
- Last residue in a chain, because psi cannot be calculated.
- Residues next to a chain break.
- Residues with unresolved alternate conformations.
- Non-standard residues.
- Residues whose previous or next neighbor is non-standard when that neighbor is needed for phi or psi.
- Residues from unsupported records or malformed inputs that the parser cannot interpret deterministically.

Skipped residues should not appear in the primary clean CSV. A separate skipped-residue report can be added later.

## Alternate Conformation Policy

The implementation must not silently choose arbitrary coordinates when alternate conformations are present.

Recommended initial deterministic policy:

1. Prefer blank alternate location identifier if present.
2. Otherwise choose the alternate location with highest occupancy for each required atom.
3. If occupancy is tied or missing in a way that prevents deterministic selection, skip the residue.
4. If required atoms for one residue would be selected from conflicting alternate identifiers and the implementation does not explicitly support that case, skip the residue.

The chosen policy must be documented in command help and metadata where practical.

Simpler acceptable initial policy:

- Skip any residue involved in phi or psi calculation if any required atom has multiple alternate conformations.

## Chain Break Policy

The implementation must avoid calculating phi or psi across chain breaks.

Recommended initial detection:

1. Never bridge different chains or models.
2. Treat missing previous or next parsed residue as a break.
3. Treat non-adjacent polymer residues as a break when the parser exposes chain discontinuity information.
4. Optionally apply a configurable peptide continuity check using the distance between previous `C` and current `N`, and between current `C` and next `N`.

If a distance threshold is used, it must be documented as a pragmatic continuity heuristic, not a full validation method.

## Error Handling

The command should:

- Fail clearly when no input files are provided.
- Fail clearly when an input path does not exist.
- Continue processing other files when one file cannot be parsed, unless `--strict` is supplied.
- Write an empty CSV with headers if no glycine residues pass the filters.
- Generate an empty plot with the correct axes if no residues pass the filters, unless the user disables plot generation.
- Return a non-zero exit code for invalid arguments or complete failure to read all inputs.

## Proposed Folder Structure

```text
GlyPhiPsi/
  README.md
  PROJECT_SPEC.md
  pyproject.toml
  src/
    glyphipsi/
      __init__.py
      cli.py
      io.py
      angles.py
      residues.py
      plotting.py
  tests/
    data/
      minimal_gly.pdb
      missing_atom_gly.pdb
      chain_break_gly.pdb
      altloc_gly.pdb
      nonstandard_residue.pdb
    test_angles.py
    test_parser_filters.py
    test_cli_outputs.py
```

## Suggested Module Responsibilities

### `cli.py`

- Parse command-line arguments.
- Expand input paths.
- Coordinate parsing, CSV writing, and plotting.
- Return clear exit codes.

### `io.py`

- Load PDB and mmCIF files through a structure parser.
- Expose normalized model, chain, residue, and atom records to the rest of the application.
- Handle parser exceptions and file-level errors.

### `residues.py`

- Identify standard residues.
- Select deterministic atom coordinates under the configured alternate conformation policy.
- Detect whether a glycine residue has valid previous and next residues.
- Apply skip rules.

### `angles.py`

- Calculate dihedral angles from four 3D points.
- Normalize angles to the documented output range.
- Keep geometry math independent from file parsing.

### `plotting.py`

- Read accepted angle rows or receive them in memory.
- Generate the Ramachandran scatter plot.
- Enforce axis limits and labels.

## Test Cases

### Angle Math

1. Calculates a known dihedral angle from four synthetic points.
2. Normalizes angles consistently to the documented -180 to +180 range.
3. Handles floating-point precision with explicit tolerances.

### Parser And Filtering

1. Includes a standard `GLY` residue with complete `C(i-1)`, `N`, `CA`, `C`, and `N(i+1)` atoms.
2. Skips a glycine residue missing current `N`.
3. Skips a glycine residue missing current `CA`.
4. Skips a glycine residue missing current `C`.
5. Skips a glycine residue missing previous residue `C`.
6. Skips a glycine residue missing next residue `N`.
7. Skips first residue in a chain.
8. Skips last residue in a chain.
9. Skips glycine adjacent to a chain break.
10. Skips non-standard glycine-like residue names such as `DGL`, `GLX`, or modified residues unless explicitly mapped later.
11. Skips glycine when alternate conformations cannot be resolved deterministically.
12. Selects the expected alternate conformation when the configured policy can resolve it deterministically.
13. Does not calculate across chains.
14. Handles insertion codes without merging distinct residues.
15. Handles multiple input files and preserves source file identity.

### CLI Outputs

1. Writes a CSV with the required header.
2. Writes one row per accepted glycine residue.
3. Writes no skipped residues to the primary CSV.
4. Produces a plot file with phi and psi axes fixed at -180 to +180.
5. Produces valid empty outputs when no residues pass filtering.
6. Returns a non-zero exit code for missing input files.

## Implementation Dependencies

Recommended dependencies:

- `biopython` for PDB and mmCIF parsing.
- `numpy` for vector math.
- `matplotlib` for plotting.

Recommended development dependencies:

- `pytest` for tests.
- `ruff` for linting and formatting.

The project should pin or bound dependencies only when necessary for compatibility. The specification does not require a specific parser library if equivalent behavior is implemented and tested.

## Acceptance Criteria

The initial implementation is complete when:

1. A user can run the CLI against local PDB or mmCIF files.
2. The command produces a CSV with one clean row per accepted glycine residue.
3. The command produces a scatter plot with phi and psi axes from -180 to +180 degrees.
4. Required skip rules are implemented and covered by tests.
5. The implementation makes no unsupported scientific claims about Ramachandran regions or structural quality.
6. Empty-result cases are handled gracefully.

