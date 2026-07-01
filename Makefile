.PHONY: validate optimizer layout

validate:
	python -m py_compile scripts/tma_xenium_core_optimizer.py scripts/tma_layout_generator.py
	python scripts/tma_xenium_core_optimizer.py

optimizer:
	python scripts/tma_xenium_core_optimizer.py

layout:
	python scripts/tma_layout_generator.py

layout-tma2:
	python scripts/tma_layout_generator.py --tma-id MDR_RA_TMA_002 --patient-start-index 41 --output-prefix mdr_ra_tma_02_map
