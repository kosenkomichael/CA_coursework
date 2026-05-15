from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Point:

    x: int | None = None
    y: int | None = None
    is_infinity: bool = False

    @staticmethod
    def infinity() -> "Point":

        return Point(is_infinity=True)

    def __str__(self) -> str:

        if self.is_infinity:
            return "(INF)"
        return f"({self.x}, {self.y})"


@dataclass(frozen=True)
class NumberModulo:

    number: int
    modulus: int


@dataclass(frozen=True)
class StepResult:

    point: Point
    a: int
    b: int


@dataclass(frozen=True)
class CycleResult:

    step_x: StepResult
    step_y: StepResult
