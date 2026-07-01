# Research Trace: Xenium TMA Chamber Design

## Design Question

For a Xenium experiment using synovial biopsy tissue, compare tissue microarray layouts that represent 80 patients by splitting them across 2 TMAs, with 40 patients assigned to each TMA.

The immediate design objective is descriptive and validation-oriented: maximize usable patient tissue surface within the Xenium chamber while preserving a robust orientation signature for de-arraying and avoiding avoidable sample-identification risk.

## Unit Of Analysis

- Patient: the clinically meaningful unit represented by replicated cores.
- TMA: one physical recipient block/section intended for one Xenium sample placement area.
- Core position: one physical spot in the array.
- Design scenario: one combination of chamber size, core diameter, spacing, layout, and orientation-marker strategy.

## Xenium Chamber Limitation

The optimizer uses a Xenium usable sample placement area of:

```text
22.45 mm x 10.45 mm = 234.60 mm2 per TMA
2 TMAs = 469.21 mm2 total chamber area
```

This is treated as a design constraint, not as a biological measurement. Core centers are generated only where the complete circular core fits inside this rectangle, using zero additional margin unless `--margin-mm` is specified.

## Spacing Constraint

The current design rule is:

```text
minimum edge-to-edge spacing = 1.0 mm
minimum center-to-center distance = core diameter + 1.0 mm
```

Examples:

```text
0.6 mm core -> minimum center distance 1.6 mm
1.0 mm core -> minimum center distance 2.0 mm
1.5 mm core -> minimum center distance 2.5 mm
2.0 mm core -> minimum center distance 3.0 mm
```

For hexagonal packing, the horizontal pitch is the minimum center-to-center distance and the vertical pitch is:

```text
(core diameter + spacing) * sqrt(3) / 2
```

For 0.6 mm cores with 1.0 mm edge spacing, that gives:

```text
horizontal pitch = 1.6000 mm
vertical pitch = 1.3856 mm
```

## TMA Grand Master Options Encoded

The scripts compare the TMA Grand Master-compatible core diameters currently assumed for this project:

```text
0.6 mm, 1.0 mm, 1.5 mm, 2.0 mm
```

The design scan includes:

- rectangular grid disposition
- hexagonal/staggered disposition
- compact fixed-corner empty orientation holes
- compact fixed-corner fiducial/control cores
- mixed empty-hole plus fiducial/control-core marker combinations

These are operational assumptions to confirm with the TMA Grand Master operator before construction. The code does not infer vendor availability; it only evaluates the allowed diameters supplied by `--allowed-core-diameters-mm`.

## Patient Split

Default study-level plan:

```text
total patients = 80
number of TMAs = 2
target patients per TMA = 40
TMA 1 generated IDs = P001-P040
TMA 2 generated IDs = P041-P080
```

For the current recommended optimizer run:

```text
layout = hexagonal
core diameter = 0.6 mm
orientation = 4 compact asymmetric empty diagonal holes
patient positions = 108 per TMA
total patient positions = 216 across 2 TMAs
mean cores per patient = 2.70
balanced patient replication per TMA = 12 patients with 2 cores and 28 patients with 3 cores
total core surface / total chamber area = 0.1302
```

## Orientation Logic

The orientation strategy follows the principles in Pilla et al. 2012:

- the array should have intrinsic asymmetry;
- the asymmetry source should be in a fixed corner;
- the marker pattern should remain interpretable if a few spots or peripheral positions are lost;
- full empty rows or columns should be avoided for Xenium because they consume scarce chamber positions and may be confused with tissue loss, section waviness, or de-arraying artifacts.

For orientation-only use, empty holes are preferred over control/fiducial cores. Control cores should be used only if they also serve histology, assay QC, or positive/negative control purposes.

## Measurement And Bias Considerations

The optimizer is geometric. It does not model:

- tissue exhaustion in donor blocks;
- core dropout during arraying, sectioning, staining, or Xenium processing;
- paraffin compression, section wrinkles, shear, or tearing;
- scanner/acquisition failures;
- segmentation or de-arraying failures;
- biological heterogeneity within synovial biopsies.

These are downstream technical and biological risks. The final design should be reviewed by the pathologist and the TMA operator before construction.
