# Codex Intake: Xenium TMA Area-Efficiency Comparison

## Research Question

Compare Xenium TMA design options for 80 synovial-biopsy patients split across 2 TMAs, targeting 40 patients per TMA, with emphasis on total core surface area divided by total Xenium chamber area.

## Task Type

Descriptive design optimization and method-validation support. This is not a causal analysis.

## Unit Of Analysis

- Design-level rows in the comparison table.
- Physical TMA positions within each Xenium chamber.
- Patient-level representation summarized as cores per patient within each TMA.

## Variables And Roles

- Design variables: core diameter, core disposition/layout, empty orientation-hole count, fiducial/control-core count.
- Outcome metrics: patient core count, total core surface area, total chamber area, total core surface/chamber area ratio, orientation-marker position cost, asymmetry score.
- Technical constraints: Xenium usable chamber dimensions, edge-to-edge spacing, available core diameters, marker visibility, de-arraying/orientation robustness.

## Measurement And Bias Considerations

- The optimizer uses idealized geometric positions and does not model tissue loss, section distortion, scanner artifacts, or segmentation failures directly.
- Orientation markers are assumed to be visible at de-arraying. Fiducial/control cores are treated as distinguishable from patient cores.
- Empty full rows or columns are penalized conceptually because they consume chamber positions and can create orientation ambiguity if tissue sections are wavy or disrupted.
- Final choices require histotechnologist/pathologist review against available punch size, chamber dimensions, tissue availability, and Xenium assay constraints.

## Resources Read

- `README.md`
- `tma_xenium_core_optimizer.py`
- `tma_layout_generator.py`
- `TMA_design_geomx.R`
- Pilla D, Bosisio FM, Marotta R, et al. Tissue microarray design and construction for scientific, industrial and diagnostic use. J Pathol Inform. 2012;3:42. PMCID: PMC3551499. DOI: 10.4103/2153-3539.104904.
- Local CC-BY-4.0 license source: `C:\Users\srenne\GitHubRepo\Dispensa_PaDi\LICENSE`.

## Resources Missing Or To Confirm

- Final vendor/operator constraints for the TMA Grand Master punch, recipient block, and Xenium chamber usable area.
- Whether fiducial/control cores will be visually and computationally distinguishable from patient tissue in the final scanned image.
- Patient ID list and any stratification variables needed for the final randomized map.
