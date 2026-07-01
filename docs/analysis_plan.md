# Analysis Plan

## Goal

Generate reproducible TMA design comparisons and final candidate maps for Xenium-compatible chamber dimensions.

## Inputs

- Script defaults and command-line parameters.
- Optional patient CSV for final maps.
- Optional historical `scripts/TMA_design_geomx.R` as a GeoMx-style orientation reference.

## Workflow

1. Run `python scripts/tma_xenium_core_optimizer.py`.
2. Review `results/tables/mdr_ra_tma_xenium_core_optimization_area_efficiency.csv`.
3. Confirm core diameter, spacing, orientation-marker strategy, and number of TMAs with the TMA operator.
4. Generate final per-TMA maps with `python scripts/tma_layout_generator.py`, using real patient IDs or controlled generated ID ranges.
5. Send the final CSV and PNG maps to the TMA Grand Master operator.

## Validation

Run:

```bash
make validate
```

If `make` is unavailable on Windows, run:

```bash
python -m py_compile scripts/tma_xenium_core_optimizer.py scripts/tma_layout_generator.py
python scripts/tma_xenium_core_optimizer.py
```

## Human Review

Confirm the final usable Xenium chamber dimensions, punch diameter availability, minimum allowable spacing, recipient block constraints, and whether fiducial/control cores are needed.
