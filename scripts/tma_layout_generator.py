#!/usr/bin/env python3
"""
Generate a tissue microarray (TMA) layout map as CSV and PNG.

Default design target for this project:
- 80 patients
- 0.5 mm diameter cores
- 1.0 mm minimum edge-to-edge spacing
- asymmetrical fiducial/orientation positions
- randomized and spatially stratified patient assignment

The default usable chamber is 25 x 20 mm = 500 mm^2. Change this with
--usable-width-mm and --usable-height-mm after confirming the actual instrument
or slide imageable dimensions.

Example:
    python tma_layout_generator.py --output-prefix synovial_tma

Using actual patient IDs:
    python tma_layout_generator.py --patients-csv patients.csv \
        --patient-id-column patient_id --output-prefix synovial_tma
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import random
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, Patch


Position = Dict[str, object]
PatientRecord = Dict[str, str]


def row_label_from_index(index_1based: int) -> str:
    """Convert 1 -> A, 26 -> Z, 27 -> AA."""
    if index_1based < 1:
        raise ValueError("Row index must be 1-based and positive.")
    label = ""
    n = index_1based
    while n:
        n, remainder = divmod(n - 1, 26)
        label = chr(65 + remainder) + label
    return label


def float_range(start: float, stop: float, step: float, eps: float = 1e-9) -> Iterable[float]:
    x = start
    while x <= stop + eps:
        yield x
        x += step


def generate_positions(
    usable_width_mm: float,
    usable_height_mm: float,
    core_diameter_mm: float,
    min_edge_space_mm: float,
    margin_mm: float,
    layout: str,
) -> Tuple[List[Position], float, float]:
    """Generate core center positions in mm.

    min_edge_space_mm is the requested minimum space between core edges. The
    nearest-neighbor center distance must therefore be core_diameter + spacing.
    """
    if usable_width_mm <= 0 or usable_height_mm <= 0:
        raise ValueError("Usable width and height must be positive.")
    if core_diameter_mm <= 0:
        raise ValueError("Core diameter must be positive.")
    if min_edge_space_mm < 0:
        raise ValueError("Minimum edge spacing cannot be negative.")
    if margin_mm < 0:
        raise ValueError("Margin cannot be negative.")

    pitch_x = core_diameter_mm + min_edge_space_mm
    if layout == "rectangular":
        pitch_y = pitch_x
    elif layout == "hexagonal":
        pitch_y = pitch_x * math.sqrt(3.0) / 2.0
    else:
        raise ValueError("layout must be 'rectangular' or 'hexagonal'.")

    x_min = margin_mm + core_diameter_mm / 2.0
    x_max = usable_width_mm - margin_mm - core_diameter_mm / 2.0
    y_min = margin_mm + core_diameter_mm / 2.0
    y_max = usable_height_mm - margin_mm - core_diameter_mm / 2.0

    if x_max < x_min or y_max < y_min:
        raise ValueError("No core centers fit in the requested chamber with the requested margin.")

    positions: List[Position] = []
    for row_idx, y in enumerate(float_range(y_min, y_max, pitch_y), start=1):
        row_label = row_label_from_index(row_idx)
        row_offset = 0.0
        if layout == "hexagonal" and row_idx % 2 == 0:
            row_offset = pitch_x / 2.0
        x_start = x_min + row_offset
        col_idx = 1
        for x in float_range(x_start, x_max, pitch_x):
            positions.append(
                {
                    "row_number": row_idx,
                    "row_label": row_label,
                    "col_number": col_idx,
                    "position_id": f"{row_label}{col_idx:02d}",
                    "x_mm": round(x, 6),
                    "y_mm": round(y, 6),
                }
            )
            col_idx += 1

    return positions, pitch_x, pitch_y


def nearest_unassigned_position(
    positions: Sequence[Position],
    target_x: float,
    target_y: float,
    used_indices: set,
) -> int:
    best_idx = -1
    best_dist = float("inf")
    for idx, pos in enumerate(positions):
        if idx in used_indices:
            continue
        dx = float(pos["x_mm"]) - target_x
        dy = float(pos["y_mm"]) - target_y
        dist = dx * dx + dy * dy
        if dist < best_dist:
            best_idx = idx
            best_dist = dist
    if best_idx < 0:
        raise RuntimeError("Could not find an unassigned position.")
    return best_idx


def choose_asymmetric_fiducials(
    positions: Sequence[Position],
    n_fiducials: int,
    pitch_x: float,
    pitch_y: float,
) -> Dict[int, str]:
    """Choose asymmetrical fiducial sites.

    The first fiducial is anchored near the lower-right corner. Additional
    fiducials are deliberately non-symmetric so the array orientation is visible
    even after reflection or 180-degree rotation.
    """
    if n_fiducials <= 0:
        return {}

    xs = [float(p["x_mm"]) for p in positions]
    ys = [float(p["y_mm"]) for p in positions]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)

    targets = [
        ("FID_A_lower_right_anchor", max_x, min_y),
        ("FID_B_upper_left", min_x, max_y),
        ("FID_C_upper_right_inset", max_x - pitch_x, max_y),
        ("FID_D_lower_left_offset", min_x, min_y + pitch_y),
        ("FID_E_internal", min_x + 2.0 * pitch_x, min_y + 3.0 * pitch_y),
        ("FID_F_right_mid", max_x, min_y + 0.55 * (max_y - min_y)),
        ("FID_G_left_lower_mid", min_x, min_y + 0.30 * (max_y - min_y)),
        ("FID_H_top_internal", min_x + 0.40 * (max_x - min_x), max_y),
    ]

    fiducials: Dict[int, str] = {}
    used_indices: set = set()
    for label, tx, ty in targets[:n_fiducials]:
        idx = nearest_unassigned_position(positions, tx, ty, used_indices)
        fiducials[idx] = label
        used_indices.add(idx)
    return fiducials


def load_patient_records(
    patients_csv: Optional[str],
    patient_id_column: str,
    n_patients: int,
    patient_prefix: str,
) -> List[PatientRecord]:
    if patients_csv is None:
        return [{"patient_id": f"{patient_prefix}{i:03d}"} for i in range(1, n_patients + 1)]

    records: List[PatientRecord] = []
    with open(patients_csv, "r", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError("Patient CSV is empty or has no header row.")
        if patient_id_column not in reader.fieldnames:
            raise ValueError(
                f"Column '{patient_id_column}' was not found in {patients_csv}. "
                f"Available columns: {reader.fieldnames}"
            )
        for row in reader:
            pid = str(row.get(patient_id_column, "")).strip()
            if not pid:
                continue
            row = {str(k): str(v) for k, v in row.items()}
            row["patient_id"] = pid
            records.append(row)

    if not records:
        raise ValueError("No patient records were loaded.")
    if len(records) != n_patients:
        print(
            f"WARNING: loaded {len(records)} patients from {patients_csv}; "
            f"overriding --n-patients={n_patients}."
        )
    return records


def balanced_patient_core_list(
    patient_records: Sequence[PatientRecord],
    n_patient_slots: int,
    rng: random.Random,
    max_cores_per_patient: Optional[int],
) -> Tuple[List[Dict[str, object]], Dict[str, int]]:
    n_patients = len(patient_records)
    if n_patient_slots < n_patients:
        raise ValueError(
            f"Only {n_patient_slots} patient slots fit, but {n_patients} patients were requested."
        )

    if max_cores_per_patient is not None:
        if max_cores_per_patient < 1:
            raise ValueError("--max-cores-per-patient must be at least 1.")
        n_patient_slots = min(n_patient_slots, n_patients * max_cores_per_patient)

    base_reps, extra = divmod(n_patient_slots, n_patients)
    patient_indices = list(range(n_patients))
    rng.shuffle(patient_indices)
    extra_indices = set(patient_indices[:extra])

    replicate_counts: Dict[str, int] = {}
    for i, record in enumerate(patient_records):
        replicate_counts[record["patient_id"]] = base_reps + (1 if i in extra_indices else 0)

    core_entries: List[Dict[str, object]] = []
    max_reps = max(replicate_counts.values())
    patient_index_by_id = {record["patient_id"]: i for i, record in enumerate(patient_records)}

    # Round-wise patient ordering avoids placing all replicates for a patient in
    # one consecutive block before positional assignment.
    for rep_idx in range(1, max_reps + 1):
        eligible = [
            record for record in patient_records if replicate_counts[record["patient_id"]] >= rep_idx
        ]
        rng.shuffle(eligible)
        for record in eligible:
            pid = record["patient_id"]
            core_entries.append(
                {
                    "patient_id": pid,
                    "replicate_index": rep_idx,
                    "total_replicates_for_patient": replicate_counts[pid],
                    "patient_input_order": patient_index_by_id[pid] + 1,
                    "patient_metadata": record,
                }
            )

    if len(core_entries) != n_patient_slots:
        raise RuntimeError("Internal error: patient core count does not match target slot count.")

    return core_entries, replicate_counts


def add_strata(
    positions: List[Position],
    usable_width_mm: float,
    usable_height_mm: float,
    strata_x: int,
    strata_y: int,
) -> None:
    for pos in positions:
        x_bin = min(int(float(pos["x_mm"]) / usable_width_mm * strata_x), strata_x - 1) + 1
        y_bin = min(int(float(pos["y_mm"]) / usable_height_mm * strata_y), strata_y - 1) + 1
        pos["stratum_x"] = x_bin
        pos["stratum_y"] = y_bin
        pos["stratum_id"] = f"Y{y_bin}_X{x_bin}"


def stratified_position_order(
    positions: Sequence[Position],
    candidate_indices: Sequence[int],
    rng: random.Random,
) -> List[int]:
    groups: Dict[str, List[int]] = defaultdict(list)
    for idx in candidate_indices:
        groups[str(positions[idx]["stratum_id"])].append(idx)
    for group in groups.values():
        rng.shuffle(group)
    strata = list(groups.keys())
    rng.shuffle(strata)

    ordered: List[int] = []
    remaining = True
    while remaining:
        remaining = False
        for stratum in strata:
            if groups[stratum]:
                ordered.append(groups[stratum].pop())
                remaining = True
    return ordered


def validate_minimum_spacing(
    positions: Sequence[Position],
    core_diameter_mm: float,
    min_edge_space_mm: float,
) -> Tuple[float, float]:
    min_center_distance = float("inf")
    for i in range(len(positions)):
        xi, yi = float(positions[i]["x_mm"]), float(positions[i]["y_mm"])
        for j in range(i + 1, len(positions)):
            xj, yj = float(positions[j]["x_mm"]), float(positions[j]["y_mm"])
            dist = math.hypot(xi - xj, yi - yj)
            if dist < min_center_distance:
                min_center_distance = dist
    min_edge_distance = min_center_distance - core_diameter_mm
    if min_edge_distance + 1e-4 < min_edge_space_mm:
        raise RuntimeError(
            f"Spacing validation failed: minimum edge spacing is {min_edge_distance:.4f} mm, "
            f"but requested {min_edge_space_mm:.4f} mm."
        )
    return min_center_distance, min_edge_distance


def build_layout_rows(
    positions: List[Position],
    fiducial_indices: Dict[int, str],
    fiducial_mode: str,
    patient_core_entries: Sequence[Dict[str, object]],
    patient_position_order: Sequence[int],
    args: argparse.Namespace,
) -> List[Dict[str, object]]:
    patient_assignment: Dict[int, Dict[str, object]] = {}
    for idx, core_entry in zip(patient_position_order, patient_core_entries):
        patient_assignment[idx] = core_entry

    rows: List[Dict[str, object]] = []
    for idx, pos in enumerate(positions):
        fid_label = fiducial_indices.get(idx)
        role = "patient"
        fiducial_id = ""
        patient_id = ""
        replicate_index = ""
        total_reps = ""
        patient_input_order = ""
        note = ""
        metadata: PatientRecord = {}

        if fid_label is not None:
            fiducial_id = fid_label
            role = "fiducial_core" if fiducial_mode == "cores" else "fiducial_hole"
            patient_id = fid_label if fiducial_mode == "cores" else ""
            note = "Asymmetric orientation position"
        elif idx in patient_assignment:
            entry = patient_assignment[idx]
            patient_id = str(entry["patient_id"])
            replicate_index = int(entry["replicate_index"])
            total_reps = int(entry["total_replicates_for_patient"])
            patient_input_order = int(entry["patient_input_order"])
            metadata = dict(entry["patient_metadata"])  # type: ignore[arg-type]
        else:
            role = "unused"
            note = "Unused because --max-cores-per-patient limited patient assignment"

        row: Dict[str, object] = {
            "tma_id": args.tma_id,
            "position_index": idx + 1,
            "position_id": pos["position_id"],
            "row_number": pos["row_number"],
            "row_label": pos["row_label"],
            "col_number": pos["col_number"],
            "x_mm": pos["x_mm"],
            "y_mm": pos["y_mm"],
            "core_diameter_mm": args.core_diameter_mm,
            "min_edge_space_mm": args.min_edge_space_mm,
            "layout": args.layout,
            "stratum_x": pos["stratum_x"],
            "stratum_y": pos["stratum_y"],
            "stratum_id": pos["stratum_id"],
            "role": role,
            "patient_id": patient_id,
            "replicate_index": replicate_index,
            "total_replicates_for_patient": total_reps,
            "patient_input_order": patient_input_order,
            "fiducial_id": fiducial_id,
            "random_seed": args.seed,
            "note": note,
        }
        for key, value in metadata.items():
            if key == "patient_id":
                continue
            row[f"patient_{key}"] = value
        rows.append(row)
    return rows


def write_csv(rows: Sequence[Dict[str, object]], output_csv: str) -> None:
    base_fields = [
        "tma_id",
        "position_index",
        "position_id",
        "row_number",
        "row_label",
        "col_number",
        "x_mm",
        "y_mm",
        "core_diameter_mm",
        "min_edge_space_mm",
        "layout",
        "stratum_x",
        "stratum_y",
        "stratum_id",
        "role",
        "patient_id",
        "replicate_index",
        "total_replicates_for_patient",
        "patient_input_order",
        "fiducial_id",
        "random_seed",
        "note",
    ]
    extra_fields = sorted({key for row in rows for key in row.keys()} - set(base_fields))
    fieldnames = base_fields + extra_fields
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def label_for_row(row: Dict[str, object], label_mode: str) -> str:
    if label_mode == "none":
        return ""
    if label_mode == "position":
        return str(row["position_id"])
    if label_mode == "patient":
        return str(row["patient_id"])
    if label_mode == "patient_rep":
        pid = str(row["patient_id"])
        rep = str(row["replicate_index"])
        return f"{pid}\nR{rep}" if rep else pid
    if label_mode == "role":
        role = str(row["role"])
        if role.startswith("fiducial"):
            fid = str(row["fiducial_id"])
            return fid.split("_")[1] if "_" in fid else fid
        return str(row["position_id"])
    raise ValueError(f"Unknown label mode: {label_mode}")


def plot_layout(
    rows: Sequence[Dict[str, object]],
    output_png: str,
    usable_width_mm: float,
    usable_height_mm: float,
    core_diameter_mm: float,
    args: argparse.Namespace,
    n_patient_cores: int,
    n_fiducials: int,
    min_center_distance: float,
    min_edge_distance: float,
) -> None:
    fig_w = max(8.0, usable_width_mm / 2.2)
    fig_h = max(6.0, usable_height_mm / 2.2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    ax.add_patch(
        Rectangle(
            (0, 0),
            usable_width_mm,
            usable_height_mm,
            fill=False,
            linewidth=1.5,
            edgecolor="black",
        )
    )

    for row in rows:
        x = float(row["x_mm"])
        y = float(row["y_mm"])
        role = str(row["role"])
        if role == "patient":
            face = "#d9d9d9"
            edge = "#404040"
            alpha = 0.9
        elif role == "fiducial_core":
            face = "#e41a1c"
            edge = "#7f0000"
            alpha = 0.95
        elif role == "fiducial_hole":
            face = "white"
            edge = "#e41a1c"
            alpha = 1.0
        else:
            face = "#ffffff"
            edge = "#bdbdbd"
            alpha = 0.8
        ax.add_patch(
            Circle(
                (x, y),
                radius=core_diameter_mm / 2.0,
                facecolor=face,
                edgecolor=edge,
                linewidth=0.5,
                alpha=alpha,
            )
        )
        label = label_for_row(row, args.label_mode)
        if label:
            ax.text(
                x,
                y,
                label,
                ha="center",
                va="center",
                fontsize=args.label_fontsize,
                color="black",
            )

    ax.set_xlim(-0.6, usable_width_mm + 0.6)
    ax.set_ylim(-0.6, usable_height_mm + 0.6)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x coordinate, mm")
    ax.set_ylabel("y coordinate, mm")
    ax.set_title(
        f"{args.tma_id}: {n_patient_cores} patient cores + {n_fiducials} fiducials | "
        f"d={core_diameter_mm} mm, edge spacing>={args.min_edge_space_mm} mm"
    )

    summary = (
        f"seed={args.seed}; layout={args.layout}; chamber={usable_width_mm:g} x {usable_height_mm:g} mm; "
        f"min center distance={min_center_distance:.3f} mm; min edge distance={min_edge_distance:.3f} mm"
    )
    ax.text(
        0.5,
        -0.10,
        summary,
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=8,
    )

    legend_handles = [
        Patch(facecolor="#d9d9d9", edgecolor="#404040", label="Patient core"),
        Patch(facecolor="#e41a1c", edgecolor="#7f0000", label="Fiducial/core or anchor"),
        Patch(facecolor="white", edgecolor="#bdbdbd", label="Unused/hole"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=8, frameon=True)
    fig.tight_layout()
    fig.savefig(output_png, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)


def write_metadata_txt(
    output_txt: str,
    args: argparse.Namespace,
    rows: Sequence[Dict[str, object]],
    replicate_counts: Dict[str, int],
    min_center_distance: float,
    min_edge_distance: float,
) -> None:
    patient_counts = sorted(replicate_counts.values())
    n_patients = len(replicate_counts)
    n_patient_cores = sum(1 for r in rows if r["role"] == "patient")
    n_fiducials = sum(1 for r in rows if str(r["role"]).startswith("fiducial"))
    n_unused = sum(1 for r in rows if r["role"] == "unused")
    with open(output_txt, "w") as f:
        f.write(f"TMA ID: {args.tma_id}\n")
        f.write(f"Output prefix: {args.output_prefix}\n")
        f.write(f"Random seed: {args.seed}\n")
        f.write(f"Layout: {args.layout}\n")
        f.write(f"Usable chamber: {args.usable_width_mm} x {args.usable_height_mm} mm\n")
        f.write(f"Core diameter: {args.core_diameter_mm} mm\n")
        f.write(f"Minimum edge spacing requested: {args.min_edge_space_mm} mm\n")
        f.write(f"Minimum center distance observed: {min_center_distance:.4f} mm\n")
        f.write(f"Minimum edge distance observed: {min_edge_distance:.4f} mm\n")
        f.write(f"Total physical positions: {len(rows)}\n")
        f.write(f"Patients: {n_patients}\n")
        f.write(f"Patient cores assigned: {n_patient_cores}\n")
        f.write(f"Fiducials: {n_fiducials}\n")
        f.write(f"Unused positions: {n_unused}\n")
        if patient_counts:
            f.write(
                f"Replicates per patient: min={min(patient_counts)}, "
                f"median={patient_counts[len(patient_counts)//2]}, max={max(patient_counts)}\n"
            )
        f.write("\nNotes:\n")
        f.write("- Patient assignment is randomized and interleaved across spatial strata.\n")
        f.write("- Fiducial positions create a non-symmetric orientation signature.\n")
        f.write("- Verify the actual usable chamber dimensions and available punch diameter before construction.\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a randomized, asymmetrical TMA map as PNG and CSV."
    )
    parser.add_argument("--tma-id", default="SYNOVIAL_TMA_001")
    parser.add_argument("--n-patients", type=int, default=80)
    parser.add_argument("--patient-prefix", default="P")
    parser.add_argument("--patients-csv", default=None, help="Optional CSV with patient/sample metadata.")
    parser.add_argument("--patient-id-column", default="patient_id")
    parser.add_argument("--usable-width-mm", type=float, default=25.0)
    parser.add_argument("--usable-height-mm", type=float, default=20.0)
    parser.add_argument("--margin-mm", type=float, default=0.0)
    parser.add_argument("--core-diameter-mm", type=float, default=0.6)
    parser.add_argument("--min-edge-space-mm", type=float, default=1.0)
    parser.add_argument(
        "--layout",
        choices=["hexagonal", "rectangular"],
        default="hexagonal",
        help="Hexagonal maximizes positions while preserving minimum edge spacing.",
    )
    parser.add_argument(
        "--fiducials",
        type=int,
        default=4,
        help="Number of asymmetric fiducial/orientation positions.",
    )
    parser.add_argument(
        "--fiducial-mode",
        choices=["cores", "holes"],
        default="cores",
        help="Use fiducial control cores, or leave fiducial positions empty as holes.",
    )
    parser.add_argument(
        "--max-cores-per-patient",
        type=int,
        default=None,
        help="Optional cap. Without a cap, all non-fiducial positions are assigned to patients.",
    )
    parser.add_argument("--strata-x", type=int, default=4)
    parser.add_argument("--strata-y", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20260630)
    parser.add_argument("--output-prefix", default="synovial_tma_map")
    parser.add_argument(
        "--label-mode",
        choices=["none", "position", "patient", "patient_rep", "role"],
        default="position",
    )
    parser.add_argument("--label-fontsize", type=float, default=3.8)
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.n_patients < 1:
        raise ValueError("--n-patients must be positive.")
    if args.strata_x < 1 or args.strata_y < 1:
        raise ValueError("--strata-x and --strata-y must be positive.")

    rng = random.Random(args.seed)

    positions, pitch_x, pitch_y = generate_positions(
        usable_width_mm=args.usable_width_mm,
        usable_height_mm=args.usable_height_mm,
        core_diameter_mm=args.core_diameter_mm,
        min_edge_space_mm=args.min_edge_space_mm,
        margin_mm=args.margin_mm,
        layout=args.layout,
    )
    add_strata(positions, args.usable_width_mm, args.usable_height_mm, args.strata_x, args.strata_y)
    min_center_distance, min_edge_distance = validate_minimum_spacing(
        positions, args.core_diameter_mm, args.min_edge_space_mm
    )

    fiducial_indices = choose_asymmetric_fiducials(positions, args.fiducials, pitch_x, pitch_y)
    fiducial_slot_indices = set(fiducial_indices.keys())
    patient_candidate_indices = [idx for idx in range(len(positions)) if idx not in fiducial_slot_indices]

    patient_records = load_patient_records(
        args.patients_csv, args.patient_id_column, args.n_patients, args.patient_prefix
    )
    patient_core_entries, replicate_counts = balanced_patient_core_list(
        patient_records,
        n_patient_slots=len(patient_candidate_indices),
        rng=rng,
        max_cores_per_patient=args.max_cores_per_patient,
    )
    patient_position_order = stratified_position_order(positions, patient_candidate_indices, rng)
    patient_position_order = patient_position_order[: len(patient_core_entries)]

    rows = build_layout_rows(
        positions=positions,
        fiducial_indices=fiducial_indices,
        fiducial_mode=args.fiducial_mode,
        patient_core_entries=patient_core_entries,
        patient_position_order=patient_position_order,
        args=args,
    )

    output_csv = f"{args.output_prefix}.csv"
    output_png = f"{args.output_prefix}.png"
    output_txt = f"{args.output_prefix}_metadata.txt"
    os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)
    write_csv(rows, output_csv)
    plot_layout(
        rows=rows,
        output_png=output_png,
        usable_width_mm=args.usable_width_mm,
        usable_height_mm=args.usable_height_mm,
        core_diameter_mm=args.core_diameter_mm,
        args=args,
        n_patient_cores=len(patient_core_entries),
        n_fiducials=len(fiducial_indices),
        min_center_distance=min_center_distance,
        min_edge_distance=min_edge_distance,
    )
    write_metadata_txt(output_txt, args, rows, replicate_counts, min_center_distance, min_edge_distance)

    patient_count_values = sorted(replicate_counts.values())
    print(f"Wrote CSV: {output_csv}")
    print(f"Wrote PNG: {output_png}")
    print(f"Wrote metadata: {output_txt}")
    print(f"Total positions: {len(rows)}")
    print(f"Patient cores: {len(patient_core_entries)}")
    print(f"Fiducials: {len(fiducial_indices)} ({args.fiducial_mode})")
    print(
        f"Replicates per patient: min={min(patient_count_values)}, "
        f"max={max(patient_count_values)}"
    )
    print(
        f"Minimum edge spacing observed: {min_edge_distance:.4f} mm "
        f"(requested {args.min_edge_space_mm:.4f} mm)"
    )


if __name__ == "__main__":
    main()
