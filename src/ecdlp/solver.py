from __future__ import annotations

from dataclasses import dataclass
from math import isqrt
from typing import Optional

from .curve import EllipticCurve
from .models import CycleResult, NumberModulo, Point, StepResult
from .step_calculator import StepCalculator
from .utils import solve_linear_congruence


@dataclass
class RoAlgorithmSettings:
    """Параметры для запуска алгоритма"""

    p_point: Point
    q_point: Point
    curve: EllipticCurve
    start_x: StepResult | None = None
    start_y: StepResult | None = None
    should_print: bool = False
    one_cycle_iteration_count: int | None = None
    retry_count: int | None = None


class RoAlgorithmECDLPSolver:
    """Класс для решения ECDLP rho-алгоритмом"""

    def __init__(self, settings: RoAlgorithmSettings) -> None:
        """Инициализируем алгоритм настройками

        Args:
            settings (RoAlgorithmSettings): Класс настроек для алгоритма

        Raises:
            ValueError: Если хоть одна точка не на кривой
        """

        self.settings = settings
        self.p_point = settings.p_point
        self.q_point = settings.q_point
        self.curve = settings.curve

        if not self.curve.is_point_on_curve(self.p_point):
            raise ValueError("Point P is not on curve")
        if not self.curve.is_point_on_curve(self.q_point):
            raise ValueError("Point Q is not on curve")

        self.n = self.curve.find_point_order(self.p_point)
        self.step_calculator = StepCalculator(
            self.p_point, self.q_point, self.n, self.curve
        )

    def solve_ecdlp(self) -> Optional[NumberModulo]:
        """функция запуска циклов алгоритма

        Returns:
            Optional[NumberModulo]: Класс вычетов при успехе, None, если не хватило попыток
        """

        start_x = self.settings.start_x
        start_y = self.settings.start_y

        if start_x is not None and start_y is not None:
            return self._try_solve(start_x, start_y)

        retries_left = self.settings.retry_count
        while retries_left is None or retries_left >= 0:
            result = self._try_solve(
                self.step_calculator.generate_random_zero_step(),
                self.step_calculator.generate_random_zero_step(),
            )
            if result is not None:
                return result
            if retries_left is not None:
                retries_left -= 1

        return None

    def _try_solve(
        self, step_x: StepResult, step_y: StepResult
    ) -> Optional[NumberModulo]:
        """Одна попытка решения алгоритма

        Args:
            step_x (StepResult): Начальная точка x (черепаха)
            step_y (StepResult): Начальная точка y (заяц)

        Raises:
            RuntimeError: Если ни одно коллизионное сравнение не дало результата

        Returns:
            Optional[NumberModulo]: найденное решение в виде числа по модулю
        """

        cycle_result = self._cycle(step_x, step_y)
        if cycle_result is None:
            if self.settings.should_print:
                print()
            return None

        step_x = cycle_result.step_x
        step_y = cycle_result.step_y

        congruence_result = solve_linear_congruence(
            step_x.a - step_y.a,
            step_y.b - step_x.b,
            self.n,
        )
        if congruence_result is None:
            return None

        answer = congruence_result.number
        modulus = congruence_result.modulus

        if self.settings.should_print:
            print()
            print(f"{self.p_point} * k = {self.q_point}")
            print(f"({step_x.a} - {step_y.a}) * k = ({step_y.b} - {step_x.b})")
            print(f"k = {answer} (mod {modulus})")

        while answer < self.n:
            computed_q = self.curve.multiply_point(self.p_point, answer)
            if computed_q == self.q_point:
                if self.settings.should_print:
                    print(f"Answer: k = {answer} (mod {self.n})")
                return NumberModulo(answer, self.n)
            answer += modulus

        raise RuntimeError("ECDLP has no solution")

    def _cycle(self, step_x: StepResult, step_y: StepResult) -> Optional[CycleResult]:
        """Непосредственно процесс пересчета следующего щага

        Args:
            step_x (StepResult): Текцщее значение точки x
            step_y (StepResult): Текущее значение точки y

        Returns:
            Optional[CycleResult]: Пара совпавших x,y, либо none, если не встретилось за отведенное число попыток
        """

        self._print_step(0, step_x, step_y)

        i = 1
        limit = self.settings.one_cycle_iteration_count
        while True:
            if limit is not None and i > limit:
                return None

            step_x = self.step_calculator.calculate(step_x)
            step_y = self.step_calculator.calculate(
                self.step_calculator.calculate(step_y)
            )

            self._print_step(i, step_x, step_y)

            if step_x.point == step_y.point:
                if step_x.a == step_y.a or step_x.b == step_y.b:
                    return None
                return CycleResult(step_x, step_y)

            i += 1

    def _print_step(self, i: int, step_x: StepResult, step_y: StepResult) -> None:
        """Визуализация решения алгоритма

        Args:
            i (int): номер шага
            step_x (StepResult): текущий х
            step_y (StepResult): текцщий у
        """

        if self.settings.should_print:
            print(
                f"i = {i} | X = {step_x.point} | Y = {step_y.point} | "
                f"α = {step_x.a} | β = {step_x.b} | γ = {step_y.a} | δ = {step_y.b}"
            )
