from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Point:
    """Точка на эллиптической кривой"""

    x: int | None = None
    y: int | None = None
    is_infinity: bool = False

    @staticmethod
    def infinity() -> "Point":
        """Генерирует точку на бесконечности = нейтральный элемент

        Returns:
            Point: точку
        """

        return Point(is_infinity=True)

    def __str__(self) -> str:
        """Функция для print точки

        Returns:
            str: строка вида (x,y)
        """

        if self.is_infinity:
            return "(INF)"
        return f"({self.x}, {self.y})"


@dataclass(frozen=True)
class NumberModulo:
    """Класс для записи "число по модулю" """

    number: int
    modulus: int


@dataclass(frozen=True)
class StepResult:
    """Класс для промежуточного результата одного шага rho-алгоритма"""

    point: Point
    a: int
    b: int


@dataclass(frozen=True)
class CycleResult:
    """Результат при обнаружении коллизииы"""

    step_x: StepResult
    step_y: StepResult
