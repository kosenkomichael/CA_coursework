from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Point:
    """
    Точка на эллиптической кривой над полем F_p.
    Может представлять либо конечную точку P = (x, y),
    либо особый элемент is_infinity=True, соответствующий
    точке на бесконечности O, которая является нейтральным
    элементом группы E(F_p).
    """

    x: int | None = None
    y: int | None = None
    is_infinity: bool = False

    @staticmethod
    def infinity() -> "Point":
        """
        Создаёт точку на бесконечности O.
        В групповой структуре E(F_p) точка O играет роль
        нейтрального элемента:
            P + O = O + P = P
        для любой точки P на кривой.
        """
        return Point(is_infinity=True)

    def __str__(self) -> str:
        """
        Возвращает строковое представление точки.
        Для точки на бесконечности возвращает строку "(INF)".
        Для обычной точки возвращает координаты в виде "(x, y)".
        """
        if self.is_infinity:
            return "(INF)"
        return f"({self.x}, {self.y})"


@dataclass(frozen=True)
class NumberModulo:
    """
    Представление класса вычетов x (mod m).
    Используется для хранения решения линейного сравнения вида
        x ≡ number (mod modulus),
    возникающего, например, при восстановлении скаляра k в задаче
    дискретного логарифма из уравнения
        A ≡ B k (mod n).
    """

    number: int
    modulus: int


@dataclass(frozen=True)
class StepResult:
    """
    Состояние шага в алгоритме Pollard rho для ECDLP.
    Хранит:
        point = a P + b Q ∈ E(F_p),
    где P — базовая точка, Q — целевая точка, а a, b ∈ Z_n.
    Инвариант алгоритма:
        на каждом шаге существует представление
            point_i = a_i P + b_i Q,
        и переходы (из StepCalculator) обновляют a_i, b_i
        согласованно с групповым законом.
    """

    point: Point
    a: int
    b: int


@dataclass(frozen=True)
class CycleResult:
    """
    Результат обнаружения цикла (коллизии) в Pollard rho.
    Содержит два шага X и Y такие, что:
        X.point = Y.point,
        X.point = a_x P + b_x Q,
        Y.point = a_y P + b_y Q.
    Из равенства a_x P + b_x Q = a_y P + b_y Q следует:
        (a_x - a_y) P = (b_y - b_x) Q = (b_y - b_x) k P,
    а значит возникает линейное сравнение на k:
        a_x - a_y ≡ (b_y - b_x) k (mod n),
    по которому восстанавливается дискретный логарифм k.
    """

    step_x: StepResult
    step_y: StepResult
