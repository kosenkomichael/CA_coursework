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
    """
    Параметры запуска алгоритма Pollard rho для решения ECDLP.
    Задаёт исходные данные задачи
        Q = k P
    в подгруппе, порождённой точкой P на эллиптической кривой E(F_p),
    а также служебные настройки алгоритма: начальные состояния,
    лимит итераций, число перезапусков и режим подробного вывода.
    """

    p_point: Point
    q_point: Point
    curve: EllipticCurve
    start_x: StepResult | None = None
    start_y: StepResult | None = None
    should_print: bool = False
    one_cycle_iteration_count: int | None = None
    retry_count: int | None = None


class RoAlgorithmECDLPSolver:
    """
    Решатель задачи ECDLP методом Pollard rho.
    Ищет скаляр k из соотношения
        Q = k P,
    где P и Q — точки на эллиптической кривой, а вычисления ведутся
    в циклической подгруппе <P> порядка n = ord(P).
    Алгоритм строит псевдослучайную последовательность состояний
        X_i = a_i P + b_i Q,
    ищет коллизию X_i = X_j методом Флойда и затем восстанавливает k
    из линейного сравнения
        a_i - a_j ≡ (b_j - b_i) k (mod n).
    """

    def __init__(self, settings: RoAlgorithmSettings) -> None:
        """
        Инициализирует решатель ECDLP и подготавливает вспомогательные структуры.
        Проверяет, что точки P и Q лежат на заданной кривой, вычисляет
        порядок точки P:
            n = ord(P),
        и создаёт объект StepCalculator, который поддерживает инвариант
            X_i = a_i P + b_i Q
        на каждом шаге rho-процесса.
        Args:
            settings: Настройки алгоритма и исходные данные задачи.
        Raises:
            ValueError: Если P или Q не лежат на эллиптической кривой.
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
        """
        Запускает решение задачи дискретного логарифма на эллиптической кривой.
        Если начальные состояния X0 и Y0 заданы явно, выполняется один запуск
        Pollard rho. Иначе используются случайные стартовые состояния вида
            X_0 = a_0 P + b_0 Q,
            Y_0 = c_0 P + d_0 Q,
        и при необходимости выполняются повторные попытки.
        Теоретически алгоритм Pollard rho ищет коллизию в конечной группе
        размера n за ожидаемое время порядка O(√n), поэтому перезапуски
        повышают вероятность нахождения полезной коллизии.
        Returns:
            Класс вычетов для найденного решения k, либо None, если за
            отведённое число попыток решить ECDLP не удалось.
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
        """
        Выполняет одну попытку решения ECDLP из заданных стартовых состояний.
        Сначала ищет коллизию состояний
            X = a_x P + b_x Q,
            Y = a_y P + b_y Q
        такую, что X = Y. Тогда:
            a_x P + b_x Q = a_y P + b_y Q,
        откуда после переноса и подстановки Q = k P получается сравнение
            (a_x - a_y) ≡ (b_y - b_x) k (mod n).
        Это линейное сравнение решается относительно k. Если решение задано
        как класс вычетов
            k ≡ k_0 (mod m),
        то далее перебираются все кандидаты k_0 + t m в диапазоне [0, n),
        пока не будет найдено значение, удовлетворяющее равенству Q = k P.
        Args:
            step_x: Начальное состояние "черепахи".
            step_y: Начальное состояние "зайца".
        Returns:
            Найденное решение k в виде NumberModulo, либо None, если текущий
            запуск не дал полезной коллизии или сравнение неразрешимо.
        Raises:
            RuntimeError: Если после получения сравнения ни один кандидат не
                удовлетворяет уравнению Q = k P.
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
        """
        Ищет коллизию в rho-последовательности методом Флойда.
        Используется схема "черепаха-заяц":
            X_{i+1} = f(X_i),
            Y_{i+1} = f(f(Y_i)),
        где функция перехода f сохраняет представление
            X_i = a_i P + b_i Q.
        Если на некоторой итерации выполняется
            X_i.point = Y_i.point,
        то получена коллизия в группе точек. Такая коллизия даёт основу
        для восстановления дискретного логарифма. Вырожденные случаи,
        когда совпадение не позволяет построить полезное сравнение,
        отбрасываются.
        Args:
            step_x: Начальное состояние медленной последовательности.
            step_y: Начальное состояние быстрой последовательности.
        Returns:
            Пара столкнувшихся состояний CycleResult, либо None, если лимит
            итераций исчерпан или найдена бесполезная коллизия.
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
        """
        Печатает текущее состояние итераций Pollard rho в отладочном режиме.
        Выводит номер шага и коэффициенты представления
            X_i = α P + β Q,
            Y_i = γ P + δ Q,
        что позволяет визуально отслеживать сохранение инварианта и момент
        возникновения коллизии X_i = Y_i.
        Args:
            i: Номер итерации.
            step_x: Текущее состояние медленной последовательности.
            step_y: Текущее состояние быстрой последовательности.
        """
        if self.settings.should_print:
            print(
                f"i = {i} | X = {step_x.point} | Y = {step_y.point} | "
                f"α = {step_x.a} | β = {step_x.b} | γ = {step_y.a} | δ = {step_y.b}"
            )
