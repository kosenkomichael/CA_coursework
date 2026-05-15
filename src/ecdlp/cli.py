from __future__ import annotations

import csv
from pathlib import Path
from math import isqrt
from time import perf_counter

from .curve import EllipticCurve
from .models import Point, StepResult
from .solver import RoAlgorithmECDLPSolver, RoAlgorithmSettings


def ask_int(prompt: str) -> int:
    return int(input(prompt).strip())


def ask_yes_no(prompt: str) -> bool:
    return input(prompt).strip().lower() == "y"


def main() -> None:
    print("ECDLP solver (Pollard rho)\n")

    csv_path = Path("data/global_data.csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "bits",
        "p",
        "a",
        "b",
        "px",
        "py",
        "qx",
        "qy",
        "provided",
        "elapsed_seconds",
        "one_cycle_iteration_count",
    ]

    p = ask_int("Введите p (модуль): ")
    a = ask_int("Введите a (коэффициент кривой): ")
    b = ask_int("Введите b (коэффициент кривой): ")

    curve = EllipticCurve(a, b, p)

    x_p = ask_int("\nВведите x_P: ")
    y_p = ask_int("Введите y_P: ")
    p_point = Point(x_p, y_p)

    if not curve.is_point_on_curve(p_point):
        print("Ошибка: точка P не лежит на кривой")
        return

    x_q = ask_int("\nВведите x_Q: ")
    y_q = ask_int("Введите y_Q: ")
    q_point = Point(x_q, y_q)

    if not curve.is_point_on_curve(q_point):
        print("Ошибка: точка Q не лежит на кривой")
        return

    should_print = ask_yes_no("\nРаспечатывать подробную информацию? (y/N): ")

    start_x = None
    start_y = None

    if ask_yes_no("\nЗадать начальные точки X0/Y0? (y/N): "):
        print("\nВвод X0:")
        a1 = ask_int("Введите α: ")
        b1 = ask_int("Введите β: ")

        print("\nВвод Y0:")
        a2 = ask_int("Введите γ: ")
        b2 = ask_int("Введите δ: ")

        x0 = curve.add_points(
            curve.multiply_point(p_point, a1), curve.multiply_point(q_point, b1)
        )
        y0 = curve.add_points(
            curve.multiply_point(p_point, a2), curve.multiply_point(q_point, b2)
        )

        start_x = StepResult(x0, a1, b1)
        start_y = StepResult(y0, a2, b2)

    retry_count = None
    if start_x is None and ask_yes_no("\nЗадать число попыток? (y/N): "):
        retry_count = ask_int("Введите количество попыток: ")

    if ask_yes_no("\nЗадать лимит итераций в одной попытке? (y/N): "):
        iteration_count = ask_int("Введите число итераций: ")
    else:
        iteration_count = 2 * isqrt(p)

    settings = RoAlgorithmSettings(
        p_point=p_point,
        q_point=q_point,
        curve=curve,
        start_x=start_x,
        start_y=start_y,
        should_print=should_print,
        one_cycle_iteration_count=iteration_count,
        retry_count=retry_count,
    )

    solver = RoAlgorithmECDLPSolver(settings)
    start_time = perf_counter()
    result = solver.solve_ecdlp()
    end_time = perf_counter()

    elapsed_time = end_time - start_time

    bits = p.bit_length()
    provided_number = None if result is None else result.number

    print("\n")
    if result is None:
        print("Не удалось решить задачу ECDLP")
    else:
        print(f"Ответ: k = {result.number} (mod {result.modulus})")
    print(f"Execution time: {elapsed_time:.6f} seconds")

    file_exists_and_not_empty = csv_path.exists() and csv_path.stat().st_size > 0

    with csv_path.open("a", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

        if not file_exists_and_not_empty:
            writer.writeheader()

        writer.writerow(
            {
                "bits": bits,
                "p": p,
                "a": a,
                "b": b,
                "px": x_p,
                "py": y_p,
                "qx": x_q,
                "qy": y_q,
                "provided": provided_number,
                "elapsed_seconds": f"{elapsed_time:.6f}",
                "one_cycle_iteration_count": settings.one_cycle_iteration_count,
            }
        )
