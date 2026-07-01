# Assumptions

## Xenium Chamber

- Usable chamber width: 22.45 mm.
- Usable chamber height: 10.45 mm.
- Usable chamber area: 234.60 mm2 per TMA.
- Default number of TMAs: 2.
- Total chamber area for the default design: 469.21 mm2.

## Core Geometry

- Default core diameter: 0.6 mm.
- Allowed core diameters explored: 0.6, 1.0, 1.5, and 2.0 mm.
- Default minimum edge-to-edge spacing: 1.0 mm.
- Center-to-center spacing is core diameter plus edge-to-edge spacing.
- Default margin is 0.0 mm beyond the core radius required to keep cores inside the chamber.

## Patient Allocation

- Default cohort size: 80 patients.
- Default split: 40 patients per TMA.
- Generated IDs for TMA 1 are P001-P040.
- Generated IDs for TMA 2 are P041-P080 when `--patient-start-index 41` is used.
- Real patient IDs should be provided through CSV before construction.

## Orientation

- Orientation-only markers should usually be empty holes.
- Fiducial/control cores should be used only when they add assay or histologic QC value.
- Full empty rows or columns are discouraged for Xenium chamber-constrained layouts.

## Privacy

- Generated IDs are synthetic.
- Real patient identifiers or protected health information should not be committed.
