from __future__ import annotations

import csv
from pathlib import Path
from time import perf_counter

from .curve import EllipticCurve
from .models import Point
from .solver import RoAlgorithmECDLPSolver, RoAlgorithmSettings

start_bound = 10


def run_test_file(
    path: str = "data/ro_alg_ecdlp_data.txt", csv_path: str = "data/test_output.csv"
) -> None:
    test_path = Path(path)
    if not test_path.exists():
        raise FileNotFoundError(f"Test data file not found: {test_path}")
    result_csv_path = Path(csv_path)
    result_csv_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "bits",
        "p",
        "a",
        "b",
        "px",
        "py",
        "qx",
        "qy",
        "expected",
        "provided",
        "success",
        "elapsed_seconds",
        "one_cycle_iteration_count",
    ]

    with (
        test_path.open("r", encoding="utf-8") as f,
        result_csv_path.open("w", newline="", encoding="utf-8") as csv_file,
    ):
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for idx, line in enumerate(f, start=1):
            parts = [int(x) for x in line.strip().split()]
            if len(parts) != 8:
                continue

            p, a, b, px, py, qx, qy, expected = parts
            curve = EllipticCurve(a, b, p)
            p_point = Point(px, py)
            q_point = Point(qx, qy)

            settings = RoAlgorithmSettings(
                p_point=p_point,
                q_point=q_point,
                curve=curve,
                one_cycle_iteration_count=2 * int(p**0.5),
            )

            start_time = perf_counter()
            provided = RoAlgorithmECDLPSolver(settings).solve_ecdlp()
            end_time = perf_counter()

            elapsed_time = end_time - start_time
            provided_number = None if provided is None else provided.number
            success = provided is not None and provided.number == expected
            bits = p.bit_length()

            if provided is None or provided.number != expected:
                print(
                    f"[FAIL] #{idx + start_bound - 1}: expected={expected}, got={None if provided is None else provided.number}"
                )
            else:
                print(f"[OK] #{idx + start_bound - 1}")

            print(f"Execution time: {elapsed_time:.6f} seconds")
            writer.writerow(
                {
                    "bits": bits,
                    "p": p,
                    "a": a,
                    "b": b,
                    "px": px,
                    "py": py,
                    "qx": qx,
                    "qy": qy,
                    "expected": expected,
                    "provided": provided_number,
                    "success": success,
                    "elapsed_seconds": f"{elapsed_time:.6f}",
                    "one_cycle_iteration_count": settings.one_cycle_iteration_count,
                }
            )
