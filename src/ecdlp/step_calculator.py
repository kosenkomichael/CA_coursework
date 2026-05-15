from __future__ import annotations

import random

from .curve import EllipticCurve
from .models import Point, StepResult


class StepCalculator:
    """
    Вычисляет переходы rho-последовательности в алгоритме Pollard rho для ECDLP.
    Поддерживает инвариант представления состояния:
        X_i = a_i P + b_i Q,
    где P — базовая точка, Q = kP — целевая точка, а все коэффициенты
    берутся по модулю n = ord(P).
    Для построения псевдослучайной итерационной функции множество точек
    разбивается на три класса по x-координате, и для каждого класса
    применяется своё правило обновления. Такая схема является стандартной
    для метода Pollard rho и позволяет получать коллизии вида X_i = X_j,
    из которых затем восстанавливается дискретный логарифм k.
    """

    def __init__(
        self, p_point: Point, q_point: Point, n: int, curve: EllipticCurve
    ) -> None:
        """
        Инициализирует функцию перехода для rho-алгоритма.
        Хранит точки P и Q, порядок подгруппы
            n = ord(P),
        и параметры разбиения множества точек на три подмножества по
        x-координате:
            S_0 = {X : 0 <= x(X) < p / 3},
            S_1 = {X : p / 3 <= x(X) < 2p / 3},
            S_2 = {X : 2p / 3 <= x(X) < p}.
        Такое разбиение используется для задания итерационной функции f,
        которая должна быть достаточно "перемешивающей", но при этом
        сохранять инвариант X_i = a_i P + b_i Q.
        Args:
            p_point: Базовая точка P.
            q_point: Целевая точка Q = kP.
            n: Порядок точки P, то есть размер циклической подгруппы <P>.
            curve: Эллиптическая кривая, на которой выполняются вычисления.
        """
        self.p_point = p_point
        self.q_point = q_point
        self.n = n
        self.curve = curve

        third = curve.p // 3
        self.first_right = third
        self.second_right = 2 * third

    def generate_random_zero_step(self) -> StepResult:
        """
        Генерирует случайное начальное состояние rho-последовательности.
        Выбирает случайные коэффициенты
            a, b ∈ Z_n
        и строит точку
            X_0 = a P + b Q.
        Такое представление является основным инвариантом алгоритма Pollard rho:
        каждый шаг хранится не только как точка X_i, но и как коэффициенты
        a_i, b_i, чтобы после обнаружения коллизии X_i = X_j можно было
        вывести линейное сравнение на искомый скаляр k.
        Returns:
            Случайное начальное состояние StepResult с инвариантом
            point = a P + b Q.
        """
        a = random.randrange(0, self.n)
        b = random.randrange(0, self.n)
        x = self.curve.add_points(
            self.curve.multiply_point(self.p_point, a),
            self.curve.multiply_point(self.q_point, b),
        )
        return StepResult(x, a, b)

    def _range_of(self, point: Point) -> int:
        """
        Определяет класс точки для итерационной функции Pollard rho.
        Разбиение выполняется по x-координате точки:
            - bucket 0, если 0 <= x < p / 3;
            - bucket 1, если p / 3 <= x < 2p / 3;
            - bucket 2, иначе.
        Точка на бесконечности по соглашению относится к первому классу.
        Такое разбиение задаёт кусочно-определённую функцию перехода f(X),
        которая используется для построения псевдослучайной rho-траектории.
        Args:
            point: Точка на эллиптической кривой.

        Returns:
            Номер класса: 0, 1 или 2.
        """
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
        """
        Вычисляет следующий шаг rho-последовательности.
        Если текущее состояние имеет вид
            X = a P + b Q,
        то в зависимости от класса точки X применяется одно из правил:
        1) Если X ∈ S_0:
               X' = X + P,
               a' = a + 1,
               b' = b.
        2) Если X ∈ S_1:
               X' = 2X,
               a' = 2a,
               b' = 2b.
        3) Если X ∈ S_2:
               X' = X + Q,
               a' = a,
               b' = b + 1.
        Все коэффициенты берутся по модулю n. Эти формулы гарантируют
        сохранение инварианта:
            X' = a' P + b' Q.
        Именно это свойство делает возможным восстановление дискретного
        логарифма после нахождения коллизии в rho-последовательности.
        Args:
            previous_step: Предыдущее состояние алгоритма вида
                point = a P + b Q.
        Returns:
            Следующее состояние StepResult после применения функции перехода.
        """
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
