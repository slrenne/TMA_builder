# Xenium Synovial TMA Design Tools

This folder contains two Python scripts for designing and comparing tissue microarray (TMA) layouts for a Xenium experiment using synovial biopsy samples.

The historical `TMA_design.R` file is **not required**. It is treated only as a previous GeoMx-style reference design if you explicitly provide it with `--r-design`.

## Files to use

### 1. Final randomized TMA map generator

Use this script when you already know the design parameters and want to produce the final TMA map.

```bash
python tma_layout_generator.py
```

Main outputs:

- `<output-prefix>.csv`: randomized TMA coordinate map
- `<output-prefix>.png`: visual TMA map
- `<output-prefix>_metadata.txt`: design summary

This is the script to use for the final map to provide to the TMA Grand Master operator.

### 2. Xenium TMA strategy optimizer

Use this script to compare different design strategies before choosing the final design.

```bash
python tma_xenium_core_optimizer.py
```

Main outputs:

- `<output-prefix>_comparison.csv`: full strategy comparison table
- `<output-prefix>_best_allowed.csv`: best designs among allowed core sizes
- `<output-prefix>_recommended_map.csv`: randomized map for the recommended design
- `<output-prefix>_recommended_map.png`: PNG map for the recommended design
- `<output-prefix>_surface_vs_diameter.png`: tissue surface versus core diameter
- `<output-prefix>_cores_per_patient_vs_diameter.png`: cores per patient versus core diameter
- `<output-prefix>_allowed_sizes.png`: comparison of allowed core sizes
- `<output-prefix>_orientation_penalty.png`: effect of empty-orientation strategies
- `<output-prefix>_tma_count_scan.csv`: comparison across different numbers of TMAs
- `<output-prefix>_tma_count_scan.png`: plot comparing 1, 2, 3, and 4 TMAs
- `<output-prefix>_summary.txt`: text summary
- `<output-prefix>_outputs.zip`: compressed bundle of outputs

This is the script to use when deciding the best approach in terms of:

- total tissue surface
- total number of cores
- cores per patient
- feasibility in a single TMA or multiple TMAs
- cost of orientation holes
- rectangular versus hexagonal layout

## Requirements

Python 3 with the following packages:

```bash
pip install pandas numpy matplotlib
```

No internet connection is required to run the scripts.

## Recommended Xenium assumptions

The optimizer defaults are Xenium-oriented rather than GeoMx-oriented:

```text
usable width  = 22.45 mm
usable height = 10.45 mm
minimum core diameter = 0.6 mm
minimum edge-to-edge spacing = 1.0 mm
number of patients = 80
allowed core diameters = 0.6, 1.0, 1.5, 2.0 mm
```

The 0.5 mm core size is no longer used because the practical minimum available core size is assumed to be 0.6 mm.

## A. Produce a final randomized TMA map

Example for one Xenium TMA using 0.6 mm cores, 1 mm spacing, hexagonal packing, and asymmetric empty holes:

```bash
python tma_layout_generator.py \
  --n-patients 80 \
  --usable-width-mm 22.45 \
  --usable-height-mm 10.45 \
  --core-diameter-mm 0.6 \
  --min-edge-space-mm 1.0 \
  --layout hexagonal \
  --fiducials 4 \
  --fiducial-mode holes \
  --seed 20260630 \
  --output-prefix synovial_xenium_tma
```

This produces:

```text
synovial_xenium_tma.csv
synovial_xenium_tma.png
synovial_xenium_tma_metadata.txt
```

### Use real patient IDs

Create a CSV such as:

```csv
patient_id
P001
P002
P003
...
P080
```

Then run:

```bash
python tma_layout_generator.py \
  --patients-csv patients.csv \
  --patient-id-column patient_id \
  --usable-width-mm 22.45 \
  --usable-height-mm 10.45 \
  --core-diameter-mm 0.6 \
  --min-edge-space-mm 1.0 \
  --layout hexagonal \
  --fiducials 4 \
  --fiducial-mode holes \
  --seed 20260630 \
  --output-prefix synovial_xenium_tma
```

The `--seed` value makes the randomization reproducible. Use a different seed to generate a different randomization.

## B. Compare strategies for a single TMA

Use this command to compare strategies for one TMA:

```bash
python tma_xenium_core_optimizer.py \
  --number-of-tmas 1 \
  --n-patients 80 \
  --min-edge-space-mm 1.0 \
  --min-core-diameter-mm 0.6 \
  --allowed-core-diameters-mm 0.6,1.0,1.5,2.0 \
  --output-prefix single_tma_xenium
```

Review these files first:

```text
single_tma_xenium_best_allowed.csv
single_tma_xenium_comparison.csv
single_tma_xenium_surface_vs_diameter.png
single_tma_xenium_cores_per_patient_vs_diameter.png
single_tma_xenium_orientation_penalty.png
single_tma_xenium_summary.txt
```

For a single TMA, the optimizer will usually favor 0.6 mm cores if the goal is to maximize the number of patient cores. Larger cores provide more surface per core, but sharply reduce the number of patients/replicates that fit in the Xenium sample area.

## C. Compare strategies requiring a minimum number of cores per patient

For example, require at least 3 cores per patient across all TMAs:

```bash
python tma_xenium_core_optimizer.py \
  --number-of-tmas 3 \
  --n-patients 80 \
  --min-cores-per-patient 3 \
  --output-prefix my_xenium_tma
```

This produces a recommended map and statistics for the best feasible design under that constraint.

## D. Compare 1, 2, 3, and 4 TMAs

The optimizer automatically produces a TMA-count scan. You can control the scan with:

```bash
python tma_xenium_core_optimizer.py \
  --number-of-tmas 3 \
  --tma-counts-to-plot 1,2,3,4 \
  --min-cores-per-patient 3 \
  --output-prefix tma_count_scan
```

Review:

```text
tma_count_scan_tma_count_scan.csv
tma_count_scan_tma_count_scan.png
```

## E. Optional: compare against the old GeoMx R design

Only use this if `TMA_design.R` is in your current folder:

```bash
python tma_xenium_core_optimizer.py \
  --r-design TMA_design.R \
  --number-of-tmas 1 \
  --output-prefix xenium_with_geomx_reference
```

If the R file is not present, the optimizer continues normally.

## Orientation strategy recommendation

For Xenium, avoid a full empty row or full empty columns because the sample area is limited and every empty position costs tissue.

Preferred approach:

```text
4 to 6 compact asymmetric empty holes
```

A diagonal or staircase pattern is preferable to a straight empty row because it provides orientation asymmetry while using fewer positions and is less likely to be confused with sectioning loss, tissue dropout, or a systematic row artifact.

Use empty holes when the purpose is only orientation. Use control tissue cores only if you also need histology or assay quality-control tissue.

## Practical workflow

1. Run `tma_xenium_core_optimizer.py` to compare strategies.
2. Inspect the CSV and PNG summary outputs.
3. Choose the design: core size, layout, number of TMAs, and orientation-hole strategy.
4. Run `tma_layout_generator.py` with those final parameters and your real patient list.
5. Send the final PNG and CSV map to the TMA Grand Master operator.

## Current practical recommendation

For 80 synovial-biopsy patients in a Xenium experiment, start with:

```text
core size: 0.6 mm
layout: hexagonal
spacing: 1.0 mm edge-to-edge
orientation: 4 to 6 compact asymmetric empty holes
randomization: stratified/randomized with fixed seed
```

If at least 3 cores per patient are required, use multiple TMAs rather than increasing the core size.
