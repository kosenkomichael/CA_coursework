from __future__ import annotations

import random

from .curve import EllipticCurve
from .models import Point, StepResult


class StepCalculator:

    def __init__(
        self, p_point: Point, q_point: Point, n: int, curve: EllipticCurve
    ) -> None:

        self.p_point = p_point
        self.q_point = q_point
        self.n = n
        self.curve = curve

        third = curve.p // 3
        self.first_right = third
        self.second_right = 2 * third

    def generate_random_zero_step(self) -> StepResult:

        a = random.randrange(0, self.n)
        b = random.randrange(0, self.n)
        x = self.curve.add_points(
            self.curve.multiply_point(self.p_point, a),
            self.curve.multiply_point(self.q_point, b),
        )
        return StepResult(x, a, b)

    def _range_of(self, point: Point) -> int:

        if point.is_infinity:
            return 0
        assert point.x is not None
        x = point.x
        if 0 <= x < self.first_right:
            return 0
        if self.first_right <= x < self.second_right:
            return 1
        return 2

    def calculate(self, previous_step: StepResult) -> StepResult:

        point = previous_step.point
        bucket = self._range_of(point)

        if bucket == 0:
            return StepResult(
                self.curve.add_points(point, self.p_point),
                (previous_step.a + 1) % self.n,
                previous_step.b % self.n,
            )
        if bucket == 1:
            return StepResult(
                self.curve.add_points(point, point),
                (2 * previous_step.a) % self.n,
                (2 * previous_step.b) % self.n,
            )
        return StepResult(
            self.curve.add_points(point, self.q_point),
            previous_step.a % self.n,
            (previous_step.b + 1) % self.n,
        )
