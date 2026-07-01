#!/usr/bin/env python3
"""
Xenium-specific TMA core-size, orientation-strategy, and multi-TMA optimizer.

The uploaded R design is treated as a historical GeoMx reference only. It is
parsed only to compare its old empty-row/empty-column orientation strategy; its
geometry is NOT used to define the Xenium chamber.

Defaults target the Xenium optimal sample placement area: 22.45 x 10.45 mm.
Default TMA Grand Master-compatible core diameters are: 0.6, 1.0, 1.5, 2.0 mm.
The minimum core diameter is 0.6 mm.

Main outputs:
  results/tables/<prefix>_comparison.csv
  results/tables/<prefix>_area_efficiency.csv
  results/tables/<prefix>_best_allowed.csv
  results/figures/<prefix>_surface_vs_diameter.png
  results/figures/<prefix>_cores_per_patient_vs_diameter.png
  results/figures/<prefix>_allowed_sizes.png
  results/figures/<prefix>_orientation_penalty.png
  results/tables/<prefix>_tma_count_scan.csv
  results/figures/<prefix>_tma_count_scan.png
  results/tables/<prefix>_recommended_map.csv
  results/figures/<prefix>_recommended_map.png
  results/tables/<prefix>_summary.txt
  results/<prefix>_outputs.zip
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import random
import re
import zipfile
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Patch, Rectangle


@dataclass(frozen=True)
class RDesignInfo:
    n_rows: int
    n_cols: int
    reference_core_diameter_mm: float
    empty_rows: Tuple[int, ...]
    empty_cols: Tuple[int, ...]


@dataclass(frozen=True)
class Position:
    row: int
    col: int
    x_mm: float
    y_mm: float


@dataclass(frozen=True)
class StrategyResult:
    strategy: str
    layout: str
    number_of_tmas: int
    n_patients: int
    patients_per_tma_target: float
    core_diameter_mm: float
    usable_width_mm: float
    usable_height_mm: float
    min_edge_space_mm: float
    positions_per_tma: int
    patient_positions_per_tma: int
    empty_positions_per_tma: int
    control_positions_per_tma: int
    total_positions_all_tmas: int
    patient_positions_all_tmas: int
    empty_positions_all_tmas: int
    control_positions_all_tmas: int
    mean_cores_per_patient: float
    min_balanced_cores_per_patient: int
    patients_with_extra_core: int
    patient_tissue_area_mm2: float
    all_tissue_area_mm2: float
    total_chamber_area_mm2: float
    patient_tissue_fraction_of_chamber_area_all_tmas: float
    total_core_surface_fraction_of_chamber_area_all_tmas: float
    tissue_fraction_of_sample_area_all_tmas: float
    orientation_penalty_percent: float
    asymmetry_score: int
    feasible: bool


STRATEGY_LABELS = {
    "none": "No orientation holes",
    "four_empty_row": "4 empty holes in one row",
    "four_empty_diagonal": "4 empty diagonal holes",
    "six_empty_staircase": "6 empty staircase holes",
    "four_control_asymmetric": "4 asymmetric control cores",
    "r_empty_row_cols": "Old R/GeoMx empty row + columns",
}

DEFAULT_STRATEGIES = [
    "none",
    "four_empty_row",
    "four_empty_diagonal",
    "six_empty_staircase",
    "four_control_asymmetric",
    "r_empty_row_cols",
]


def default_r_design_info(fallback_core_diameter_mm: float = 0.6) -> RDesignInfo:
    return RDesignInfo(
        n_rows=14,
        n_cols=31,
        reference_core_diameter_mm=fallback_core_diameter_mm,
        empty_rows=(11,),
        empty_cols=(16, 27),
    )


def parse_float_list(text: str) -> List[float]:
    return [float(x) for x in re.split(r"[,; ]+", text.strip()) if x]


def parse_int_list(text: str) -> List[int]:
    return [int(x) for x in re.split(r"[,; ]+", text.strip()) if x]


def patients_per_tma_target(n_patients: int, number_of_tmas: int, override: Optional[float] = None) -> float:
    if override is not None:
        return float(override)
    if number_of_tmas <= 0:
        return float("nan")
    return n_patients / number_of_tmas


def output_path(output_prefix: str, suffix: str, default_dir: str) -> str:
    prefix_dir = os.path.dirname(output_prefix)
    filename = f"{os.path.basename(output_prefix)}{suffix}"
    return os.path.join(prefix_dir or default_dir, filename)


def parse_strategy_list(text: str) -> List[str]:
    if text.strip().lower() == "all":
        return list(DEFAULT_STRATEGIES)
    out = [x.strip() for x in text.split(",") if x.strip()]
    unknown = sorted(set(out) - set(DEFAULT_STRATEGIES))
    if unknown:
        raise ValueError(f"Unknown strategies: {unknown}. Valid: {DEFAULT_STRATEGIES}")
    return out


def read_r_design(path: str, fallback_core_diameter_mm: float = 0.6) -> RDesignInfo:
    text = open(path, "r", encoding="utf-8", errors="replace").read()
    matrix_match = re.search(
        r"matrix\s*\(\s*nrow\s*=\s*(\d+)\s*,\s*ncol\s*=\s*(\d+)",
        text,
        flags=re.IGNORECASE | re.DOTALL,
    )
    if not matrix_match:
        raise ValueError("Could not find matrix(nrow = ..., ncol = ...) in the R file.")
    n_rows = int(matrix_match.group(1))
    n_cols = int(matrix_match.group(2))

    core_match = re.search(
        r"chosen\s+design\s+is\s*([0-9]+(?:\.[0-9]+)?)\s*cores",
        text,
        flags=re.IGNORECASE,
    )
    reference_core = float(core_match.group(1)) if core_match else fallback_core_diameter_mm

    empty_cols: List[int] = []
    for match in re.finditer(r"TMA_map\s*\[\s*,\s*c\(([^)]*)\)\s*\]\s*<-\s*['\"]E['\"]", text):
        empty_cols.extend(int(x) for x in re.findall(r"\d+", match.group(1)))
    for match in re.finditer(r"TMA_map\s*\[\s*,\s*(\d+)\s*\]\s*<-\s*['\"]E['\"]", text):
        empty_cols.append(int(match.group(1)))

    empty_rows: List[int] = []
    for match in re.finditer(r"TMA_map\s*\[\s*(\d+)\s*,\s*\]\s*<-\s*['\"]E['\"]", text):
        empty_rows.append(int(match.group(1)))

    return RDesignInfo(
        n_rows=n_rows,
        n_cols=n_cols,
        reference_core_diameter_mm=reference_core,
        empty_rows=tuple(sorted(set(empty_rows))),
        empty_cols=tuple(sorted(set(empty_cols))),
    )


def resolve_r_design(path: Optional[str], fallback_core_diameter_mm: float = 0.6) -> RDesignInfo:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    search_paths = [path] if path else [
        "TMA_design.R",
        "TMA_design_geomx.R",
        os.path.join(script_dir, "TMA_design.R"),
        os.path.join(script_dir, "TMA_design_geomx.R"),
    ]
    for candidate in search_paths:
        if candidate and os.path.exists(candidate):
            return read_r_design(candidate, fallback_core_diameter_mm)
    if path:
        print(f"Warning: R design file not found: {path}. Using built-in 14 x 31 GeoMx orientation reference.")
    else:
        print("Warning: no R design file found. Using built-in 14 x 31 GeoMx orientation reference.")
    return default_r_design_info(fallback_core_diameter_mm)


def generate_positions(
    usable_width_mm: float,
    usable_height_mm: float,
    core_diameter_mm: float,
    min_edge_space_mm: float,
    margin_mm: float,
    layout: str,
) -> List[Position]:
    pitch_x = core_diameter_mm + min_edge_space_mm
    pitch_y = pitch_x if layout == "rectangular" else pitch_x * math.sqrt(3.0) / 2.0
    x_min = margin_mm + core_diameter_mm / 2.0
    x_max = usable_width_mm - margin_mm - core_diameter_mm / 2.0
    y_min = margin_mm + core_diameter_mm / 2.0
    y_max = usable_height_mm - margin_mm - core_diameter_mm / 2.0
    positions: List[Position] = []
    if x_min > x_max or y_min > y_max:
        return positions
    row = 1
    y = y_min
    while y <= y_max + 1e-9:
        offset = pitch_x / 2.0 if layout == "hexagonal" and row % 2 == 0 else 0.0
        x = x_min + offset
        col = 1
        while x <= x_max + 1e-9:
            positions.append(Position(row=row, col=col, x_mm=x, y_mm=y))
            col += 1
            x += pitch_x
        row += 1
        y += pitch_y
    return positions


def nearest_indices(positions: Sequence[Position], targets: Sequence[Tuple[float, float]]) -> List[int]:
    used = set()
    chosen: List[int] = []
    for tx, ty in targets:
        best_idx = None
        best_dist = float("inf")
        for idx, pos in enumerate(positions):
            if idx in used:
                continue
            dist = (pos.x_mm - tx) ** 2 + (pos.y_mm - ty) ** 2
            if dist < best_dist:
                best_idx = idx
                best_dist = dist
        if best_idx is not None:
            chosen.append(best_idx)
            used.add(best_idx)
    return chosen


def target_points(positions: Sequence[Position], strategy: str, pitch: float) -> List[Tuple[float, float]]:
    xs = [p.x_mm for p in positions]
    ys = [p.y_mm for p in positions]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    width = max_x - min_x
    height = max_y - min_y
    if strategy == "four_empty_row":
        return [(max_x - i * pitch, min_y) for i in range(4)]
    if strategy == "four_empty_diagonal":
        return [(max_x - i * pitch, min_y + i * pitch) for i in range(4)]
    if strategy in {"six_empty_staircase", "four_control_asymmetric"}:
        pts = [
            (max_x, min_y),
            (max_x - pitch, min_y),
            (max_x - pitch, min_y + pitch),
            (max_x - 2.0 * pitch, min_y + pitch),
            (max_x - 2.0 * pitch, min_y + 2.0 * pitch),
            (min_x + 0.34 * width, min_y + 0.72 * height),
        ]
        return pts[:4] if strategy == "four_control_asymmetric" else pts
    return []


def compact_corner_targets(positions: Sequence[Position], marker_count: int, pitch: float) -> List[Tuple[float, float]]:
    if marker_count <= 0 or not positions:
        return []
    xs = [p.x_mm for p in positions]
    ys = [p.y_mm for p in positions]
    max_x = max(xs)
    min_y = min(ys)
    points: List[Tuple[float, float]] = []
    for idx in range(marker_count):
        step = idx // 2
        if idx % 2 == 0:
            points.append((max_x - step * pitch, min_y + step * pitch))
        else:
            points.append((max_x - (step + 1) * pitch, min_y + step * pitch))
    return points


def assign_compact_orientation_roles(
    positions: Sequence[Position],
    empty_count: int,
    fiducial_count: int,
    core_diameter_mm: float,
    min_edge_space_mm: float,
) -> List[str]:
    roles = ["patient" for _ in positions]
    marker_count = max(0, empty_count) + max(0, fiducial_count)
    if not positions or marker_count <= 0:
        return roles
    pitch = core_diameter_mm + min_edge_space_mm
    marker_indices = nearest_indices(positions, compact_corner_targets(positions, marker_count, pitch))
    for idx in marker_indices[:empty_count]:
        roles[idx] = "empty"
    for idx in marker_indices[empty_count:empty_count + fiducial_count]:
        roles[idx] = "control"
    return roles


def scaled_r_empty_indices(positions: Sequence[Position], r_design: RDesignInfo) -> Tuple[set, set]:
    n_rows = max(p.row for p in positions) if positions else 0
    n_cols = max(p.col for p in positions) if positions else 0
    rows = {max(1, min(n_rows, int(round(r / r_design.n_rows * n_rows)))) for r in r_design.empty_rows}
    cols = {max(1, min(n_cols, int(round(c / r_design.n_cols * n_cols)))) for c in r_design.empty_cols}
    return rows, cols


def assign_roles(
    positions: Sequence[Position],
    strategy: str,
    r_design: RDesignInfo,
    core_diameter_mm: float,
    min_edge_space_mm: float,
) -> List[str]:
    roles = ["patient" for _ in positions]
    if not positions or strategy == "none":
        return roles
    if strategy == "r_empty_row_cols":
        rows, cols = scaled_r_empty_indices(positions, r_design)
        for idx, pos in enumerate(positions):
            if pos.row in rows or pos.col in cols:
                roles[idx] = "empty"
        return roles
    pitch = core_diameter_mm + min_edge_space_mm
    marker_indices = nearest_indices(positions, target_points(positions, strategy, pitch))
    marker_role = "control" if strategy == "four_control_asymmetric" else "empty"
    for idx in marker_indices:
        roles[idx] = marker_role
    return roles


def asymmetry_score(positions: Sequence[Position], roles: Sequence[str]) -> int:
    if not positions:
        return 0
    n_rows = max(p.row for p in positions)
    n_cols = max(p.col for p in positions)
    role_values = {"empty": 0, "patient": 1, "control": 2}
    m = [[-1 for _ in range(n_cols)] for _ in range(n_rows)]
    for pos, role in zip(positions, roles):
        m[pos.row - 1][pos.col - 1] = role_values.get(role, 1)

    def rotate180(mat: List[List[int]]) -> List[List[int]]:
        return [list(reversed(row)) for row in reversed(mat)]

    def rotate90(mat: List[List[int]]) -> List[List[int]]:
        return [list(row) for row in zip(*mat[::-1])]

    def rotate270(mat: List[List[int]]) -> List[List[int]]:
        return [list(row) for row in zip(*mat)][::-1]

    transforms = [
        [list(reversed(row)) for row in m],
        list(reversed(m)),
        rotate180(m),
    ]
    if n_rows == n_cols:
        transforms.extend([
            rotate90(m),
            rotate270(m),
            [list(row) for row in zip(*m)],
            [list(row) for row in zip(*rotate180(m))],
        ])
    return min(sum(1 for r in range(n_rows) for c in range(n_cols) if m[r][c] != t[r][c]) for t in transforms)


def has_full_empty_line(positions: Sequence[Position], roles: Sequence[str]) -> bool:
    row_counts: Dict[int, int] = {}
    col_counts: Dict[int, int] = {}
    row_empty: Dict[int, int] = {}
    col_empty: Dict[int, int] = {}
    for pos, role in zip(positions, roles):
        row_counts[pos.row] = row_counts.get(pos.row, 0) + 1
        col_counts[pos.col] = col_counts.get(pos.col, 0) + 1
        if role == "empty":
            row_empty[pos.row] = row_empty.get(pos.row, 0) + 1
            col_empty[pos.col] = col_empty.get(pos.col, 0) + 1
    return any(row_empty.get(row, 0) == count for row, count in row_counts.items()) or any(
        col_empty.get(col, 0) == count for col, count in col_counts.items()
    )


def evaluate(
    r_design: RDesignInfo,
    strategies: Sequence[str],
    layouts: Sequence[str],
    diameters: Sequence[float],
    number_of_tmas: int,
    usable_width_mm: float,
    usable_height_mm: float,
    min_edge_space_mm: float,
    margin_mm: float,
    n_patients: int,
    min_cores_per_patient: int,
) -> List[StrategyResult]:
    results: List[StrategyResult] = []
    total_sample_area = usable_width_mm * usable_height_mm * number_of_tmas
    patients_per_tma = patients_per_tma_target(n_patients, number_of_tmas)
    for layout in layouts:
        for diameter in diameters:
            positions = generate_positions(usable_width_mm, usable_height_mm, diameter, min_edge_space_mm, margin_mm, layout)
            for strategy in strategies:
                roles = assign_roles(positions, strategy, r_design, diameter, min_edge_space_mm)
                patient_positions_per = roles.count("patient")
                empty_positions_per = roles.count("empty")
                control_positions_per = roles.count("control")
                positions_per = len(positions)
                patient_positions_all = patient_positions_per * number_of_tmas
                empty_positions_all = empty_positions_per * number_of_tmas
                control_positions_all = control_positions_per * number_of_tmas
                total_positions_all = positions_per * number_of_tmas
                core_area = math.pi * (diameter / 2.0) ** 2
                patient_tissue_area = patient_positions_all * core_area
                all_core_surface_area = (patient_positions_all + control_positions_all) * core_area
                patient_fraction = patient_tissue_area / total_sample_area if total_sample_area else 0.0
                all_core_fraction = all_core_surface_area / total_sample_area if total_sample_area else 0.0
                results.append(
                    StrategyResult(
                        strategy=strategy,
                        layout=layout,
                        number_of_tmas=number_of_tmas,
                        n_patients=n_patients,
                        patients_per_tma_target=patients_per_tma,
                        core_diameter_mm=diameter,
                        usable_width_mm=usable_width_mm,
                        usable_height_mm=usable_height_mm,
                        min_edge_space_mm=min_edge_space_mm,
                        positions_per_tma=positions_per,
                        patient_positions_per_tma=patient_positions_per,
                        empty_positions_per_tma=empty_positions_per,
                        control_positions_per_tma=control_positions_per,
                        total_positions_all_tmas=total_positions_all,
                        patient_positions_all_tmas=patient_positions_all,
                        empty_positions_all_tmas=empty_positions_all,
                        control_positions_all_tmas=control_positions_all,
                        mean_cores_per_patient=patient_positions_all / n_patients,
                        min_balanced_cores_per_patient=patient_positions_all // n_patients,
                        patients_with_extra_core=patient_positions_all % n_patients,
                        patient_tissue_area_mm2=patient_tissue_area,
                        all_tissue_area_mm2=all_core_surface_area,
                        total_chamber_area_mm2=total_sample_area,
                        patient_tissue_fraction_of_chamber_area_all_tmas=patient_fraction,
                        total_core_surface_fraction_of_chamber_area_all_tmas=all_core_fraction,
                        tissue_fraction_of_sample_area_all_tmas=patient_fraction,
                        orientation_penalty_percent=(empty_positions_per + control_positions_per) / positions_per * 100.0 if positions_per else 0.0,
                        asymmetry_score=asymmetry_score(positions, roles),
                        feasible=patient_positions_all >= n_patients * min_cores_per_patient,
                    )
                )
    return results


def frange(start: float, stop: float, step: float) -> Iterable[float]:
    n = int(math.floor((stop - start) / step + 1e-9))
    for i in range(n + 1):
        yield round(start + i * step, 6)


def result_to_dict(r: StrategyResult) -> Dict[str, object]:
    return {
        "strategy": r.strategy,
        "strategy_label": STRATEGY_LABELS[r.strategy],
        "layout": r.layout,
        "number_of_tmas": r.number_of_tmas,
        "n_patients": r.n_patients,
        "patients_per_tma_target": f"{r.patients_per_tma_target:.4f}",
        "core_diameter_mm": f"{r.core_diameter_mm:.4f}",
        "usable_width_mm": f"{r.usable_width_mm:.4f}",
        "usable_height_mm": f"{r.usable_height_mm:.4f}",
        "min_edge_space_mm": f"{r.min_edge_space_mm:.4f}",
        "positions_per_tma": r.positions_per_tma,
        "patient_positions_per_tma": r.patient_positions_per_tma,
        "empty_positions_per_tma": r.empty_positions_per_tma,
        "control_positions_per_tma": r.control_positions_per_tma,
        "total_positions_all_tmas": r.total_positions_all_tmas,
        "patient_positions_all_tmas": r.patient_positions_all_tmas,
        "empty_positions_all_tmas": r.empty_positions_all_tmas,
        "control_positions_all_tmas": r.control_positions_all_tmas,
        "mean_cores_per_patient": f"{r.mean_cores_per_patient:.4f}",
        "min_balanced_cores_per_patient": r.min_balanced_cores_per_patient,
        "patients_with_extra_core": r.patients_with_extra_core,
        "patient_tissue_area_mm2": f"{r.patient_tissue_area_mm2:.4f}",
        "all_tissue_area_mm2": f"{r.all_tissue_area_mm2:.4f}",
        "total_chamber_area_mm2": f"{r.total_chamber_area_mm2:.4f}",
        "patient_tissue_fraction_of_chamber_area_all_tmas": f"{r.patient_tissue_fraction_of_chamber_area_all_tmas:.6f}",
        "total_core_surface_fraction_of_chamber_area_all_tmas": f"{r.total_core_surface_fraction_of_chamber_area_all_tmas:.6f}",
        "total_core_surface_to_total_chamber_area": f"{r.total_core_surface_fraction_of_chamber_area_all_tmas:.6f}",
        "tissue_fraction_of_sample_area_all_tmas": f"{r.tissue_fraction_of_sample_area_all_tmas:.6f}",
        "orientation_penalty_percent": f"{r.orientation_penalty_percent:.4f}",
        "asymmetry_score": r.asymmetry_score,
        "feasible": int(r.feasible),
    }


def write_csv(results: Sequence[StrategyResult], path: str) -> None:
    fields = list(result_to_dict(results[0]).keys()) if results else []
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for r in results:
            writer.writerow(result_to_dict(r))


def write_area_efficiency_comparison(args: argparse.Namespace, layouts: Sequence[str], path: str) -> None:
    patients_per_tma = patients_per_tma_target(args.n_patients, args.number_of_tmas, args.patients_per_tma)
    total_chamber_area = args.usable_width_mm * args.usable_height_mm * args.number_of_tmas
    diameters = sorted(set(parse_float_list(args.allowed_core_diameters_mm)))
    empty_counts = sorted(set(parse_int_list(args.orientation_empty_counts)))
    fiducial_counts = sorted(set(parse_int_list(args.orientation_fiducial_counts)))
    rows: List[Dict[str, object]] = []
    for layout in layouts:
        core_disposition = "hexagonal_staggered" if layout == "hexagonal" else "rectangular_grid"
        for diameter in diameters:
            positions = generate_positions(
                args.usable_width_mm,
                args.usable_height_mm,
                diameter,
                args.min_edge_space_mm,
                args.margin_mm,
                layout,
            )
            core_area = math.pi * (diameter / 2.0) ** 2
            for empty_count in empty_counts:
                for fiducial_count in fiducial_counts:
                    roles = assign_compact_orientation_roles(
                        positions,
                        empty_count,
                        fiducial_count,
                        diameter,
                        args.min_edge_space_mm,
                    )
                    patient_positions_per = roles.count("patient")
                    empty_positions_per = roles.count("empty")
                    fiducial_positions_per = roles.count("control")
                    marker_count = empty_positions_per + fiducial_positions_per
                    patient_positions_all = patient_positions_per * args.number_of_tmas
                    empty_positions_all = empty_positions_per * args.number_of_tmas
                    fiducial_positions_all = fiducial_positions_per * args.number_of_tmas
                    patient_surface = patient_positions_all * core_area
                    fiducial_surface = fiducial_positions_all * core_area
                    total_core_surface = patient_surface + fiducial_surface
                    full_empty_line = has_full_empty_line(positions, roles)
                    score = asymmetry_score(positions, roles)
                    finite_patients_per_tma = math.isfinite(patients_per_tma) and patients_per_tma > 0
                    required_patient_positions = patients_per_tma * args.min_cores_per_patient if finite_patients_per_tma else float("inf")
                    feasible = patient_positions_per + 1e-9 >= required_patient_positions
                    passes_orientation = (
                        marker_count >= args.min_orientation_markers
                        and score >= args.min_orientation_asymmetry_score
                        and not full_empty_line
                    )
                    notes: List[str] = []
                    if marker_count < empty_count + fiducial_count:
                        notes.append("requested marker count exceeds available positions")
                    if marker_count == 0:
                        notes.append("no orientation markers")
                    if marker_count < args.min_orientation_markers:
                        notes.append(f"below minimum marker count {args.min_orientation_markers}")
                    if score < args.min_orientation_asymmetry_score:
                        notes.append("asymmetry score below threshold")
                    if full_empty_line:
                        notes.append("contains full empty row/column")
                    if not notes:
                        notes.append("compact asymmetric fixed-corner marker pattern")
                    mean_cores_per_patient = patient_positions_per / patients_per_tma if finite_patients_per_tma else float("nan")
                    min_balanced = math.floor(mean_cores_per_patient) if finite_patients_per_tma else 0
                    patients_per_tma_is_integer = finite_patients_per_tma and float(patients_per_tma).is_integer()
                    patients_with_extra = patient_positions_per % int(patients_per_tma) if patients_per_tma_is_integer else ""
                    rows.append({
                        "number_of_tmas": args.number_of_tmas,
                        "total_patients": args.n_patients,
                        "patients_per_tma": f"{patients_per_tma:.4f}",
                        "core_disposition": core_disposition,
                        "layout": layout,
                        "core_diameter_mm": f"{diameter:.4f}",
                        "orientation_empty_holes_per_tma": empty_positions_per,
                        "orientation_fiducial_cores_per_tma": fiducial_positions_per,
                        "orientation_markers_per_tma": marker_count,
                        "positions_per_tma": len(positions),
                        "patient_positions_per_tma": patient_positions_per,
                        "patient_positions_all_tmas": patient_positions_all,
                        "empty_positions_all_tmas": empty_positions_all,
                        "fiducial_positions_all_tmas": fiducial_positions_all,
                        "total_core_positions_all_tmas": patient_positions_all + fiducial_positions_all,
                        "mean_cores_per_patient": f"{mean_cores_per_patient:.4f}" if finite_patients_per_tma else "",
                        "min_balanced_cores_per_patient": min_balanced,
                        "patients_with_extra_core_per_tma": patients_with_extra,
                        "core_area_mm2": f"{core_area:.6f}",
                        "patient_core_surface_area_mm2": f"{patient_surface:.4f}",
                        "fiducial_core_surface_area_mm2": f"{fiducial_surface:.4f}",
                        "total_core_surface_area_mm2": f"{total_core_surface:.4f}",
                        "total_chamber_area_mm2": f"{total_chamber_area:.4f}",
                        "patient_core_surface_to_total_chamber_area": f"{patient_surface / total_chamber_area:.6f}" if total_chamber_area else "",
                        "total_core_surface_to_total_chamber_area": f"{total_core_surface / total_chamber_area:.6f}" if total_chamber_area else "",
                        "reserved_orientation_position_percent": f"{100.0 * marker_count / len(positions):.4f}" if positions else "",
                        "asymmetry_score": score,
                        "passes_orientation_principles": int(passes_orientation),
                        "feasible_for_patients_per_tma": int(feasible),
                        "minimum_valid_orientation_marker_count": "",
                        "is_minimum_valid_orientation_combo": 0,
                        "notes": "; ".join(notes),
                    })

    minimum_by_layout_diameter: Dict[Tuple[str, str], int] = {}
    for row in rows:
        if row["passes_orientation_principles"] and row["feasible_for_patients_per_tma"]:
            key = (str(row["layout"]), str(row["core_diameter_mm"]))
            marker_count = int(row["orientation_markers_per_tma"])
            minimum_by_layout_diameter[key] = min(marker_count, minimum_by_layout_diameter.get(key, marker_count))
    for row in rows:
        key = (str(row["layout"]), str(row["core_diameter_mm"]))
        minimum = minimum_by_layout_diameter.get(key)
        if minimum is not None:
            row["minimum_valid_orientation_marker_count"] = minimum
            row["is_minimum_valid_orientation_combo"] = int(
                bool(row["passes_orientation_principles"])
                and bool(row["feasible_for_patients_per_tma"])
                and int(row["orientation_markers_per_tma"]) == minimum
            )

    fieldnames = [
        "number_of_tmas",
        "total_patients",
        "patients_per_tma",
        "core_disposition",
        "layout",
        "core_diameter_mm",
        "orientation_empty_holes_per_tma",
        "orientation_fiducial_cores_per_tma",
        "orientation_markers_per_tma",
        "positions_per_tma",
        "patient_positions_per_tma",
        "patient_positions_all_tmas",
        "empty_positions_all_tmas",
        "fiducial_positions_all_tmas",
        "total_core_positions_all_tmas",
        "mean_cores_per_patient",
        "min_balanced_cores_per_patient",
        "patients_with_extra_core_per_tma",
        "core_area_mm2",
        "patient_core_surface_area_mm2",
        "fiducial_core_surface_area_mm2",
        "total_core_surface_area_mm2",
        "total_chamber_area_mm2",
        "patient_core_surface_to_total_chamber_area",
        "total_core_surface_to_total_chamber_area",
        "reserved_orientation_position_percent",
        "asymmetry_score",
        "passes_orientation_principles",
        "feasible_for_patients_per_tma",
        "minimum_valid_orientation_marker_count",
        "is_minimum_valid_orientation_combo",
        "notes",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def best_allowed(results: Sequence[StrategyResult], allowed: set) -> List[StrategyResult]:
    grouped: Dict[Tuple[str, str], List[StrategyResult]] = {}
    for r in results:
        if r.feasible and round(r.core_diameter_mm, 6) in allowed:
            grouped.setdefault((r.strategy, r.layout), []).append(r)
    out = []
    for vals in grouped.values():
        out.append(max(vals, key=lambda r: (r.min_balanced_cores_per_patient, r.patient_tissue_area_mm2, r.patient_positions_all_tmas, r.asymmetry_score, -r.empty_positions_all_tmas)))
    return sorted(out, key=lambda r: (r.layout, -r.patient_tissue_area_mm2, r.strategy))


def plot_lines(results: Sequence[StrategyResult], metric: str, ylabel: str, title: str, path: str, min_cores: int) -> None:
    fig, ax = plt.subplots(figsize=(11, 7))
    for layout, strategy in sorted(set((r.layout, r.strategy) for r in results)):
        sub = sorted([r for r in results if r.layout == layout and r.strategy == strategy], key=lambda r: r.core_diameter_mm)
        x = [r.core_diameter_mm for r in sub]
        y = [getattr(r, metric) if r.feasible else float("nan") for r in sub]
        ax.plot(x, y, linewidth=1.6, label=f"{layout}: {STRATEGY_LABELS[strategy]}")
    ax.set_xlabel("Core diameter, mm")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if metric == "mean_cores_per_patient":
        for h in [1, 2, 3, 4]:
            ax.axhline(h, linestyle="--", linewidth=0.8)
        ax.text(0.61, min_cores + 0.02, f"minimum requested = {min_cores}", fontsize=8)
    ax.grid(True, linewidth=0.3)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_allowed(best: Sequence[StrategyResult], path: str) -> None:
    labels = [f"{r.layout}\n{STRATEGY_LABELS[r.strategy]}\nd={r.core_diameter_mm:g}" for r in best]
    values = [r.patient_tissue_area_mm2 for r in best]
    fig, ax = plt.subplots(figsize=(max(10, 0.75 * len(labels)), 6.5))
    ax.bar(range(len(labels)), values)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Total patient tissue surface area across all TMAs, mm^2")
    ax.set_title("Best feasible standard core-size choice per strategy")
    for i, r in enumerate(best):
        ax.text(i, values[i], f"{r.patient_positions_all_tmas} slots\n{r.mean_cores_per_patient:.2f}/pt", ha="center", va="bottom", fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_orientation_penalty(best: Sequence[StrategyResult], path: str) -> None:
    labels = [f"{r.layout}\n{STRATEGY_LABELS[r.strategy]}" for r in best]
    values = [r.orientation_penalty_percent for r in best]
    fig, ax = plt.subplots(figsize=(max(10, 0.75 * len(labels)), 6.5))
    ax.bar(range(len(labels)), values)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Reserved orientation/control positions per TMA, %")
    ax.set_title("Position cost of orientation strategies")
    for i, r in enumerate(best):
        ax.text(i, values[i], f"{r.empty_positions_per_tma + r.control_positions_per_tma}/{r.positions_per_tma}", ha="center", va="bottom", fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def recreate(r: StrategyResult, r_design: RDesignInfo, margin_mm: float) -> Tuple[List[Position], List[str]]:
    pos = generate_positions(r.usable_width_mm, r.usable_height_mm, r.core_diameter_mm, r.min_edge_space_mm, margin_mm, r.layout)
    return pos, assign_roles(pos, r.strategy, r_design, r.core_diameter_mm, r.min_edge_space_mm)


def assign_patients_across_tmas(roles: Sequence[str], number_of_tmas: int, n_patients: int, seed: int) -> Dict[Tuple[int, int], str]:
    rng = random.Random(seed)
    slots: List[Tuple[int, int]] = []
    for tma_index in range(1, number_of_tmas + 1):
        slots.extend((tma_index, i) for i, role in enumerate(roles) if role == "patient")
    rng.shuffle(slots)
    reps = len(slots) // n_patients
    extra = len(slots) % n_patients
    patient_ids = [f"P{i:03d}" for i in range(1, n_patients + 1)]
    ids: List[str] = []
    for idx, pid in enumerate(patient_ids):
        ids.extend([pid] * (reps + (1 if idx < extra else 0)))
    rng.shuffle(ids)
    return {slot: pid for slot, pid in zip(slots, ids)}


def assign_patients_by_tma_groups(
    roles: Sequence[str],
    number_of_tmas: int,
    n_patients: int,
    patients_per_tma: float,
    seed: int,
) -> Dict[Tuple[int, int], str]:
    if not (math.isfinite(patients_per_tma) and patients_per_tma > 0 and float(patients_per_tma).is_integer()):
        print("Warning: patients per TMA is not an integer. Using global patient distribution across all TMAs.")
        return assign_patients_across_tmas(roles, number_of_tmas, n_patients, seed)
    patients_per_tma_int = int(patients_per_tma)
    if patients_per_tma_int * number_of_tmas != n_patients:
        print("Warning: patients per TMA does not multiply to total patients. Using global patient distribution across all TMAs.")
        return assign_patients_across_tmas(roles, number_of_tmas, n_patients, seed)
    rng = random.Random(seed)
    lookup: Dict[Tuple[int, int], str] = {}
    patient_slots = [idx for idx, role in enumerate(roles) if role == "patient"]
    for tma_index in range(1, number_of_tmas + 1):
        slots = [(tma_index, idx) for idx in patient_slots]
        rng.shuffle(slots)
        start = (tma_index - 1) * patients_per_tma_int + 1
        patient_ids = [f"P{i:03d}" for i in range(start, start + patients_per_tma_int)]
        reps = len(slots) // patients_per_tma_int
        extra = len(slots) % patients_per_tma_int
        ids: List[str] = []
        for idx, pid in enumerate(patient_ids):
            ids.extend([pid] * (reps + (1 if idx < extra else 0)))
        rng.shuffle(ids)
        lookup.update({slot: pid for slot, pid in zip(slots, ids)})
    return lookup


def write_map_csv(positions: Sequence[Position], roles: Sequence[str], patient_lookup: Dict[Tuple[int, int], str], r: StrategyResult, path: str) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["tma_index", "row", "col", "x_mm", "y_mm", "role", "patient_id", "core_diameter_mm", "min_edge_space_mm", "layout", "strategy", "number_of_tmas"])
        writer.writeheader()
        for tma_index in range(1, r.number_of_tmas + 1):
            for idx, (pos, role) in enumerate(zip(positions, roles)):
                writer.writerow({
                    "tma_index": tma_index,
                    "row": pos.row,
                    "col": pos.col,
                    "x_mm": f"{pos.x_mm:.4f}",
                    "y_mm": f"{pos.y_mm:.4f}",
                    "role": role,
                    "patient_id": patient_lookup.get((tma_index, idx), ""),
                    "core_diameter_mm": f"{r.core_diameter_mm:.4f}",
                    "min_edge_space_mm": f"{r.min_edge_space_mm:.4f}",
                    "layout": r.layout,
                    "strategy": r.strategy,
                    "number_of_tmas": r.number_of_tmas,
                })


def plot_layout(positions: Sequence[Position], roles: Sequence[str], patient_lookup: Dict[Tuple[int, int], str], r: StrategyResult, path: str, annotate: bool = False) -> None:
    fig_width = min(18, max(9, 5 * r.number_of_tmas))
    fig, axes = plt.subplots(1, r.number_of_tmas, figsize=(fig_width, 5.8), squeeze=False)
    radius = r.core_diameter_mm / 2.0
    for tma_index in range(1, r.number_of_tmas + 1):
        ax = axes[0][tma_index - 1]
        ax.add_patch(Rectangle((0, 0), r.usable_width_mm, r.usable_height_mm, fill=False, linewidth=1.2))
        for idx, (pos, role) in enumerate(zip(positions, roles)):
            if role == "patient":
                face, edge = "0.80", "0.25"
            elif role == "empty":
                face, edge = "white", "0.0"
            else:
                face, edge = "0.35", "0.0"
            ax.add_patch(Circle((pos.x_mm, pos.y_mm), radius=radius, facecolor=face, edgecolor=edge, linewidth=0.5))
            pid = patient_lookup.get((tma_index, idx), "")
            if annotate and role == "patient" and pid:
                ax.text(pos.x_mm, pos.y_mm, pid[1:], ha="center", va="center", fontsize=3.5)
        ax.set_xlim(-0.6, r.usable_width_mm + 0.6)
        ax.set_ylim(-0.6, r.usable_height_mm + 0.6)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x, mm")
        ax.set_ylabel("y, mm")
        ax.set_title(f"TMA {tma_index}")
    handles = [
        Patch(facecolor="0.80", edgecolor="0.25", label="Patient core"),
        Patch(facecolor="white", edgecolor="0.0", label="Empty orientation hole"),
        Patch(facecolor="0.35", edgecolor="0.0", label="Control/fiducial core"),
    ]
    fig.legend(handles=handles, fontsize=8, loc="upper center", ncol=3)
    fig.suptitle(f"Recommended Xenium TMA map | {r.layout}, {STRATEGY_LABELS[r.strategy]}, d={r.core_diameter_mm:g} mm", y=1.03)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def choose_recommended(best: Sequence[StrategyResult]) -> StrategyResult:
    oriented = [r for r in best if r.strategy != "none"]
    candidates = oriented if oriented else list(best)
    preferred_strategy_rank = {
        "four_empty_diagonal": 5,
        "four_control_asymmetric": 4,
        "six_empty_staircase": 3,
        "four_empty_row": 2,
        "r_empty_row_cols": 1,
        "none": 0,
    }
    return max(
        candidates,
        key=lambda r: (
            r.min_balanced_cores_per_patient,
            r.patient_tissue_area_mm2,
            r.patient_positions_all_tmas,
            -(r.empty_positions_all_tmas + r.control_positions_all_tmas),
            preferred_strategy_rank[r.strategy],
            r.asymmetry_score,
        ),
    )


def scan_tma_counts(args: argparse.Namespace, r_design: RDesignInfo, strategies: Sequence[str], layouts: Sequence[str], allowed: set, path_csv: str, path_png: str) -> None:
    rows: List[StrategyResult] = []
    allowed_diameters = sorted(allowed)
    for nt in parse_int_list(args.tma_counts_to_plot):
        results = evaluate(r_design, strategies, layouts, allowed_diameters, nt, args.usable_width_mm, args.usable_height_mm, args.min_edge_space_mm, args.margin_mm, args.n_patients, args.min_cores_per_patient)
        best = best_allowed(results, allowed)
        if best:
            rows.append(choose_recommended(best))
    write_csv(rows, path_csv)
    fig, ax = plt.subplots(figsize=(9, 6))
    x = [r.number_of_tmas for r in rows]
    y = [r.mean_cores_per_patient for r in rows]
    ax.plot(x, y, marker="o", linewidth=1.8)
    for r in rows:
        ax.text(r.number_of_tmas, r.mean_cores_per_patient, f"d={r.core_diameter_mm:g}\n{r.layout}\n{STRATEGY_LABELS[r.strategy].split()[0]}", ha="center", va="bottom", fontsize=8)
    for h in [1, 2, 3, 4]:
        ax.axhline(h, linestyle="--", linewidth=0.8)
    ax.set_xlabel("Number of TMAs / Xenium sample areas")
    ax.set_ylabel("Mean patient cores per patient")
    ax.set_title("Effect of number of TMAs on patient-level replication")
    ax.grid(True, linewidth=0.3)
    fig.tight_layout()
    fig.savefig(path_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_summary(path: str, r_design: RDesignInfo, best: Sequence[StrategyResult], recommended: StrategyResult, args: argparse.Namespace) -> None:
    r_total = r_design.n_rows * r_design.n_cols
    r_empty = r_design.n_rows * len(r_design.empty_cols) + r_design.n_cols * len(r_design.empty_rows) - len(r_design.empty_rows) * len(r_design.empty_cols)
    patients_per_tma = patients_per_tma_target(args.n_patients, args.number_of_tmas, args.patients_per_tma)
    with open(path, "w") as f:
        f.write("Xenium multi-TMA core-size optimizer summary\n")
        f.write("==========================================\n\n")
        f.write("This run is Xenium-specific. The uploaded R file is used only as a historical GeoMx reference for the old empty-row/empty-column orientation strategy.\n\n")
        f.write(f"Xenium sample placement area used per TMA: {args.usable_width_mm:g} x {args.usable_height_mm:g} mm = {args.usable_width_mm * args.usable_height_mm:.2f} mm2\n")
        f.write(f"Number of TMAs / Xenium sample areas: {args.number_of_tmas}\n")
        f.write(f"Patients: {args.n_patients}\n")
        f.write(f"Target patients per TMA: {patients_per_tma:.2f}\n")
        f.write(f"Minimum edge-to-edge spacing: {args.min_edge_space_mm:g} mm\n")
        f.write(f"Allowed core diameters: {args.allowed_core_diameters_mm} mm\n")
        f.write(f"Minimum core diameter explored: {args.min_core_diameter_mm:g} mm\n")
        f.write(f"Minimum cores/patient constraint across all TMAs: {args.min_cores_per_patient}\n\n")
        f.write("Orientation principles encoded from Pilla et al. 2012 (PMCID: PMC3551499): use a visible intrinsically asymmetric pattern, keep the asymmetry source in a fixed corner, avoid full empty rows/lines, and prefer compact markers that remain interpretable if a few spots or peripheral rows are lost.\n\n")
        f.write("Historical R design, not used as Xenium geometry:\n")
        f.write(f"- Matrix: {r_design.n_rows} x {r_design.n_cols} = {r_total} positions\n")
        f.write(f"- Empty rows: {list(r_design.empty_rows)}; empty columns: {list(r_design.empty_cols)}\n")
        f.write(f"- Empty-position penalty in that design: {r_empty}/{r_total} = {100*r_empty/r_total:.1f}%\n\n")
        f.write("Best feasible allowed-size choices:\n")
        for r in best:
            f.write(f"- {r.layout}, {STRATEGY_LABELS[r.strategy]}, d={r.core_diameter_mm:g} mm: {r.patient_positions_per_tma} patient slots/TMA, {r.patient_positions_all_tmas} total patient slots, mean {r.mean_cores_per_patient:.2f}/patient, total patient area {r.patient_tissue_area_mm2:.2f} mm2, orientation/control cost {r.orientation_penalty_percent:.1f}%\n")
        f.write("\nRecommended map:\n")
        f.write(f"- {recommended.number_of_tmas} TMA(s), {recommended.layout}, {STRATEGY_LABELS[recommended.strategy]}, d={recommended.core_diameter_mm:g} mm\n")
        f.write(f"- {recommended.patient_positions_per_tma} patient cores/TMA, {recommended.patient_positions_all_tmas} total patient cores, mean {recommended.mean_cores_per_patient:.2f}/patient, minimum balanced cores/patient {recommended.min_balanced_cores_per_patient}\n")
        f.write(f"- Total core surface / total chamber area: {recommended.total_core_surface_fraction_of_chamber_area_all_tmas:.4f}; patient core surface / total chamber area: {recommended.patient_tissue_fraction_of_chamber_area_all_tmas:.4f}\n")
        f.write("- For one Xenium sample area and 80 patients, 0.6 mm is usually required to keep all patients represented with 1 mm spacing. More TMAs increase patient-level replication while keeping the same per-slide constraints.\n")
        f.write("- Prefer compact asymmetric empty-hole markers over the historical full row/column empty pattern for Xenium, because the Xenium sample area is smaller and orientation holes are expensive.\n")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Xenium-specific multi-TMA core optimizer")
    p.add_argument("--r-design", default=None, help="Optional historical R TMA design. If omitted, local TMA_design.R or TMA_design_geomx.R is used when present.")
    p.add_argument("--output-prefix", default="mdr_ra_tma_xenium_core_optimization")
    p.add_argument("--results-dir", default="results")
    p.add_argument("--tables-dir", default=None)
    p.add_argument("--figures-dir", default=None)
    p.add_argument("--n-patients", type=int, default=80)
    p.add_argument("--number-of-tmas", type=int, default=2, help="Number of TMAs/Xenium sample areas. Default is 2 for 40 patients per TMA with 80 total patients.")
    p.add_argument("--patients-per-tma", type=float, default=None, help="Target patients assigned to each TMA. Defaults to n_patients / number_of_tmas.")
    p.add_argument("--tma-counts-to-plot", default="1,2,3,4", help="Comma-separated number_of_tmas values to compare in the scan plot.")
    p.add_argument("--min-cores-per-patient", type=int, default=1, help="Minimum cores per patient across all TMAs.")
    p.add_argument("--min-edge-space-mm", type=float, default=1.0)
    p.add_argument("--margin-mm", type=float, default=0.0)
    p.add_argument("--usable-width-mm", type=float, default=22.45)
    p.add_argument("--usable-height-mm", type=float, default=10.45)
    p.add_argument("--fallback-r-core-diameter-mm", type=float, default=0.6)
    p.add_argument("--min-core-diameter-mm", type=float, default=0.6)
    p.add_argument("--max-core-diameter-mm", type=float, default=2.0)
    p.add_argument("--diameter-step-mm", type=float, default=0.01)
    p.add_argument("--allowed-core-diameters-mm", default="0.6,1.0,1.5,2.0")
    p.add_argument("--orientation-empty-counts", default="0,1,2,3,4,5,6", help="Empty orientation-hole counts to scan in the area-efficiency comparison.")
    p.add_argument("--orientation-fiducial-counts", default="0,1,2,3,4", help="Fiducial/control-core counts to scan in the area-efficiency comparison.")
    p.add_argument("--min-orientation-markers", type=int, default=4, help="Minimum visible empty+fiducial markers required for the compact corner orientation comparison.")
    p.add_argument("--min-orientation-asymmetry-score", type=int, default=1, help="Minimum asymmetry score required for the compact corner orientation comparison.")
    p.add_argument("--strategies", default="all")
    p.add_argument("--layouts", default="rectangular,hexagonal")
    p.add_argument("--seed", type=int, default=20260630)
    p.add_argument("--annotate-map", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    tables_dir = args.tables_dir or os.path.join(args.results_dir, "tables")
    figures_dir = args.figures_dir or os.path.join(args.results_dir, "figures")
    os.makedirs(tables_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(args.results_dir, exist_ok=True)
    r_design = resolve_r_design(args.r_design, args.fallback_r_core_diameter_mm)
    strategies = parse_strategy_list(args.strategies)
    layouts = [x.strip() for x in args.layouts.split(",") if x.strip()]
    diameters = sorted(set(list(frange(args.min_core_diameter_mm, args.max_core_diameter_mm, args.diameter_step_mm)) + parse_float_list(args.allowed_core_diameters_mm)))
    allowed = {round(x, 6) for x in parse_float_list(args.allowed_core_diameters_mm)}

    results = evaluate(r_design, strategies, layouts, diameters, args.number_of_tmas, args.usable_width_mm, args.usable_height_mm, args.min_edge_space_mm, args.margin_mm, args.n_patients, args.min_cores_per_patient)
    best = best_allowed(results, allowed)
    if not best:
        raise RuntimeError("No feasible design found. Try increasing --number-of-tmas, lowering --min-cores-per-patient, or reducing spacing/margins.")
    recommended = choose_recommended(best)
    rec_pos, rec_roles = recreate(recommended, r_design, args.margin_mm)
    patients_per_tma = patients_per_tma_target(args.n_patients, args.number_of_tmas, args.patients_per_tma)
    patient_lookup = assign_patients_by_tma_groups(rec_roles, args.number_of_tmas, args.n_patients, patients_per_tma, args.seed)

    outputs = {
        "comparison": output_path(args.output_prefix, "_comparison.csv", tables_dir),
        "area_efficiency": output_path(args.output_prefix, "_area_efficiency.csv", tables_dir),
        "best": output_path(args.output_prefix, "_best_allowed.csv", tables_dir),
        "surface": output_path(args.output_prefix, "_surface_vs_diameter.png", figures_dir),
        "cores": output_path(args.output_prefix, "_cores_per_patient_vs_diameter.png", figures_dir),
        "allowed": output_path(args.output_prefix, "_allowed_sizes.png", figures_dir),
        "penalty": output_path(args.output_prefix, "_orientation_penalty.png", figures_dir),
        "tma_count_csv": output_path(args.output_prefix, "_tma_count_scan.csv", tables_dir),
        "tma_count_png": output_path(args.output_prefix, "_tma_count_scan.png", figures_dir),
        "map_csv": output_path(args.output_prefix, "_recommended_map.csv", tables_dir),
        "map_png": output_path(args.output_prefix, "_recommended_map.png", figures_dir),
        "summary": output_path(args.output_prefix, "_summary.txt", tables_dir),
        "zip": output_path(args.output_prefix, "_outputs.zip", args.results_dir),
    }
    for path in outputs.values():
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    write_csv(results, outputs["comparison"])
    write_area_efficiency_comparison(args, layouts, outputs["area_efficiency"])
    write_csv(best, outputs["best"])
    plot_lines(results, "patient_tissue_area_mm2", "Total patient tissue surface area across all TMAs, mm^2", f"{args.number_of_tmas} Xenium TMA(s): tissue surface by core size", outputs["surface"], args.min_cores_per_patient)
    plot_lines(results, "mean_cores_per_patient", "Mean patient cores per patient across all TMAs", f"{args.number_of_tmas} Xenium TMA(s): replication trade-off", outputs["cores"], args.min_cores_per_patient)
    plot_allowed(best, outputs["allowed"])
    plot_orientation_penalty(best, outputs["penalty"])
    scan_tma_counts(args, r_design, strategies, layouts, allowed, outputs["tma_count_csv"], outputs["tma_count_png"])
    write_map_csv(rec_pos, rec_roles, patient_lookup, recommended, outputs["map_csv"])
    plot_layout(rec_pos, rec_roles, patient_lookup, recommended, outputs["map_png"], annotate=args.annotate_map)
    write_summary(outputs["summary"], r_design, best, recommended, args)

    with zipfile.ZipFile(outputs["zip"], "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.write(__file__, arcname=os.path.basename(__file__))
        for path in outputs.values():
            if path != outputs["zip"] and os.path.exists(path):
                zf.write(path, arcname=os.path.basename(path))

    print(f"Recommended: {recommended.number_of_tmas} TMA(s), {recommended.layout}, {STRATEGY_LABELS[recommended.strategy]}, d={recommended.core_diameter_mm:g} mm")
    print(f"Patient slots: {recommended.patient_positions_all_tmas}; mean cores/patient: {recommended.mean_cores_per_patient:.2f}")
    print(f"Patients per TMA target: {patients_per_tma:.2f}")
    print(f"Total core surface / total chamber area: {recommended.total_core_surface_fraction_of_chamber_area_all_tmas:.4f}")
    for key, path in outputs.items():
        print(f"Wrote {key}: {path}")


if __name__ == "__main__":
    main()
