from __future__ import annotations

import subprocess

from .curve import EllipticCurve
from .models import Point
from .solver import RoAlgorithmSettings


class SageProvider:
    @staticmethod
    def _execute_sage_code(code: str) -> str:
        result = subprocess.run(
            ["sage", "-c", code],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip()

    @staticmethod
    def find_elliptic_curve_order(curve: EllipticCurve) -> int:
        code = (
            f"p={curve.p}; a={curve.a}; b={curve.b}; "
            f"E=EllipticCurve(GF(p), [a, b]); print(E.order())"
        )
        return int(SageProvider._execute_sage_code(code))

    @staticmethod
    def find_point_order(curve: EllipticCurve, point: Point) -> int:
        if point.is_infinity:
            return 1
        code = (
            f"p={curve.p}; a={curve.a}; b={curve.b}; "
            f"E=EllipticCurve(GF(p), [a, b]); "
            f"P=E({point.x}, {point.y}); print(P.order())"
        )
        return int(SageProvider._execute_sage_code(code))

    @staticmethod
    def solve_ecdlp(settings: RoAlgorithmSettings) -> int:
        curve = settings.curve
        p_point = settings.p_point
        q_point = settings.q_point
        code = (
            f"p={curve.p}; a={curve.a}; b={curve.b}; "
            f"E=EllipticCurve(GF(p), [a, b]); "
            f"P=E({p_point.x}, {p_point.y}); "
            f"Q=E({q_point.x}, {q_point.y}); "
            f"print(Q.log(P))"
        )
        return int(SageProvider._execute_sage_code(code))
