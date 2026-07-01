# Xenium Synovial TMA Builder

This repository contains reproducible design tools for planning Xenium-compatible tissue microarrays (TMAs) for a synovial biopsy cohort.

The current study design is:

```text
total patients = 80
number of TMAs = 2
target patients per TMA = 40
usable Xenium chamber per TMA = 22.45 x 10.45 mm
default core diameter = 0.6 mm
minimum edge-to-edge spacing = 1.0 mm
default layout = hexagonal
default orientation strategy = compact asymmetric empty holes
```

The repository is organized so that source scripts stay at the top level, research trace files live in `docs/`, tabular outputs live in `results/tables/`, and figure outputs live in `results/figures/`.

## Scientific Posture

This is a descriptive design-optimization and method-validation project. It is not a causal, diagnostic-performance, prognostic, or treatment-effect analysis.

The clinically meaningful unit is the patient. Core-level and TMA-position-level rows are physical design units used to create patient-level representation. Do not interpret core-level counts as independent patient-level observations in downstream analyses.

## Repository Layout

```text
.
|-- README.md
|-- LICENSE
|-- Makefile
|-- TMA_design_geomx.R
|-- tma_xenium_core_optimizer.py
|-- tma_layout_generator.py
|-- docs/
|   |-- research_trace.md
|   |-- research_brief.md
|   |-- analysis_plan.md
|   |-- assumptions.md
|   |-- resources.md
|   |-- data_dictionary.csv
|   |-- resource_gaps.md
|   `-- codex_intake.md
`-- results/
    |-- tables/
    |-- figures/
    `-- tma_xenium_core_optimization_outputs.zip
```

## Design Constraints

The optimizer uses a Xenium chamber of `22.45 x 10.45 mm`, or `234.60 mm2` per TMA. With 2 TMAs, the total chamber area is `469.21 mm2`.

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

Confirm final punch availability, minimum spacing, recipient block constraints, and Xenium sample-area constraints with the TMA Grand Master operator before construction.

## Orientation Principle

Orientation markers follow the asymmetry principles described by Pilla et al. 2012 (PMCID: PMC3551499):

- use an intrinsically asymmetric design;
- keep the asymmetry source in a fixed corner;
- avoid full empty rows or columns in the Xenium chamber;
- prefer compact marker patterns that remain interpretable if a few spots are lost.

For orientation-only use, empty holes are preferred. Use fiducial/control tissue cores only when they also serve assay QC, positive/negative control, or histologic control purposes.

## Main Scripts

### Strategy Optimizer

Run the design comparison:

```bash
python tma_xenium_core_optimizer.py
```

Primary outputs:

```text
results/tables/tma_xenium_core_optimization_area_efficiency.csv
results/tables/tma_xenium_core_optimization_comparison.csv
results/tables/tma_xenium_core_optimization_best_allowed.csv
results/tables/tma_xenium_core_optimization_recommended_map.csv
results/tables/tma_xenium_core_optimization_summary.txt
results/figures/tma_xenium_core_optimization_recommended_map.png
results/figures/tma_xenium_core_optimization_surface_vs_diameter.png
results/figures/tma_xenium_core_optimization_cores_per_patient_vs_diameter.png
```

The most useful first table is:

```text
results/tables/tma_xenium_core_optimization_area_efficiency.csv
```

It compares core diameter, rectangular versus hexagonal disposition, empty holes, fiducial/control cores, patient slots, replicated cores per patient, and total core surface divided by total chamber area.

### Final Map Generator

Generate a final single-TMA map for 40 patients:

```bash
python tma_layout_generator.py
```

Default outputs:

```text
results/tables/synovial_tma_01_map.csv
results/tables/synovial_tma_01_map_metadata.txt
results/figures/synovial_tma_01_map.png
```

Generate TMA 2 with synthetic IDs P041-P080:

```bash
python tma_layout_generator.py ^
  --tma-id SYNOVIAL_TMA_002 ^
  --patient-start-index 41 ^
  --output-prefix synovial_tma_02_map
```

On bash-like shells, replace `^` with `\`.

## Using Real Patient IDs

Create a CSV with no protected health information:

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
python tma_layout_generator.py ^
  --patients-csv patients_tma_01.csv ^
  --patient-id-column patient_id ^
  --output-prefix synovial_tma_01_map
```

Use a second 40-patient CSV for TMA 2. Keep the patient-to-TMA split auditable and avoid committing identifiable patient data.

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

## Validation

Run:

```bash
make validate
```

If `make` is unavailable:

```bash
python -m py_compile tma_xenium_core_optimizer.py tma_layout_generator.py
python tma_xenium_core_optimizer.py
```

## Documentation

Read these before final block construction:

- `docs/research_trace.md`
- `docs/assumptions.md`
- `docs/resource_gaps.md`
- `results/tables/tma_xenium_core_optimization_summary.txt`

## License

This repository is licensed under the Creative Commons Attribution 4.0 International License. See `LICENSE`.
