# MDR_RA TMA Builder

This repository designs Xenium-compatible tissue microarrays (TMAs) for the MDR_RA consortium.

The default study design is:

```text
consortium/project = MDR_RA
total patients = 80
number of TMAs = 2
target patients per TMA = 40
usable Xenium chamber per TMA = 22.45 x 10.45 mm
default core diameter = 0.6 mm
minimum edge-to-edge spacing = 1.0 mm
default layout = hexagonal
default orientation strategy = compact asymmetric empty holes
```

The repository is organized from the true repo root, not from `scripts/`:

```text
.
|-- README.md
|-- LICENSE
|-- Makefile
|-- docs/
|-- results/
|   |-- tables/
|   `-- figures/
`-- scripts/
    |-- TMA_design_geomx.R
    |-- tma_xenium_core_optimizer.py
    `-- tma_layout_generator.py
```

## Scientific Posture

This is a descriptive design-optimization and method-validation repository. It is not a causal, diagnostic-performance, prognostic, or treatment-effect analysis.

The clinically meaningful unit is the patient. Core-level and TMA-position-level rows are physical design units used to create patient-level representation. Do not interpret core-level counts as independent patient-level observations in downstream biological analyses.

## Chamber And Spacing Constraints

The optimizer uses a Xenium chamber of `22.45 x 10.45 mm`, or `234.60 mm2` per TMA. With two TMAs, the total chamber area is `469.21 mm2`.

The spacing rule is:

```text
minimum center-to-center distance = core diameter + minimum edge-to-edge spacing
```

With the default `0.6 mm` core diameter and `1.0 mm` edge-to-edge spacing:

```text
minimum center-to-center distance = 1.6 mm
hexagonal vertical pitch = 1.6 * sqrt(3) / 2 = 1.3856 mm
```

The TMA Grand Master-compatible diameter options currently encoded for comparison are:

```text
0.6 mm, 1.0 mm, 1.5 mm, 2.0 mm
```

These are project assumptions. Confirm punch availability, minimum local spacing, recipient block constraints, and final Xenium sample-area constraints with the TMA Grand Master operator before construction.

## Orientation Principle

Orientation markers follow the asymmetry principles described by Pilla et al. 2012 (PMCID: PMC3551499):

- use an intrinsically asymmetric design;
- keep the asymmetry source in a fixed corner;
- avoid full empty rows or columns in the Xenium chamber;
- prefer compact marker patterns that remain interpretable if a few spots are lost.

For orientation-only use, empty holes are preferred. Use fiducial/control tissue cores only when they also serve assay QC, positive/negative control, or histologic control purposes.

## Main Commands

Run the optimizer from the repository root:

```bash
python scripts/tma_xenium_core_optimizer.py
```

Or use:

```bash
make optimizer
```

Generate final per-TMA maps:

```bash
make layout
make layout-tma2
```

The first command generates `MDR_RA_TMA_001` with synthetic IDs `P001-P040`. The second generates `MDR_RA_TMA_002` with synthetic IDs `P041-P080`.

## Main Outputs

Primary optimizer tables:

```text
results/tables/mdr_ra_tma_xenium_core_optimization_area_efficiency.csv
results/tables/mdr_ra_tma_xenium_core_optimization_comparison.csv
results/tables/mdr_ra_tma_xenium_core_optimization_best_allowed.csv
results/tables/mdr_ra_tma_xenium_core_optimization_recommended_map.csv
results/tables/mdr_ra_tma_xenium_core_optimization_summary.txt
```

Primary figures:

```text
results/figures/mdr_ra_tma_xenium_core_optimization_recommended_map.png
results/figures/mdr_ra_tma_xenium_core_optimization_surface_vs_diameter.png
results/figures/mdr_ra_tma_xenium_core_optimization_cores_per_patient_vs_diameter.png
```

Final single-TMA maps:

```text
results/tables/mdr_ra_tma_01_map.csv
results/tables/mdr_ra_tma_01_map_metadata.txt
results/figures/mdr_ra_tma_01_map.png
results/tables/mdr_ra_tma_02_map.csv
results/tables/mdr_ra_tma_02_map_metadata.txt
results/figures/mdr_ra_tma_02_map.png
```

Compiled trace:

```text
docs/mdr_ra_tma_trace.md
```

## Current Optimizer Recommendation

The current default optimizer run recommends:

```text
2 TMAs
40 patients per TMA
0.6 mm cores
hexagonal layout
1.0 mm minimum edge-to-edge spacing
4 compact asymmetric empty diagonal holes per TMA
108 patient cores per TMA
216 patient cores across both TMAs
mean cores per patient = 2.70
total core surface / total chamber area = 0.1302
```

Each TMA has 40 unique patients, with balanced replication of 2 to 3 cores per patient.

## Using Real Patient IDs

Create two CSV files with no protected health information, one for each 40-patient TMA:

```csv
patient_id
P001
P002
P003
...
P040
```

Then run:

```bash
python scripts/tma_layout_generator.py ^
  --patients-csv patients_tma_01.csv ^
  --patient-id-column patient_id ^
  --output-prefix mdr_ra_tma_01_map
```

Use a second 40-patient CSV for TMA 2. Keep the patient-to-TMA split auditable and do not commit identifiable patient data.

## Validation

Run:

```bash
make validate
```

If `make` is unavailable:

```bash
python -m py_compile scripts/tma_xenium_core_optimizer.py scripts/tma_layout_generator.py
python scripts/tma_xenium_core_optimizer.py
```

## Documentation

Read these before final block construction:

- `docs/mdr_ra_tma_trace.md`
- `docs/research_trace.md`
- `docs/assumptions.md`
- `docs/resource_gaps.md`
- `results/tables/mdr_ra_tma_xenium_core_optimization_summary.txt`

## License

This repository is licensed under the Creative Commons Attribution 4.0 International License. See `LICENSE`.
