from __future__ import annotations

import random
from math import isqrt

from .models import Point


class EllipticCurve:
    """
    Эллиптическая кривая вида y^2 ≡ x^3 + a x + b (mod p) над полем F_p.
    Кривая задаётся уравнением:
        y^2 ≡ x^3 + a x + b (mod p),
    где p — нечётное простое, а, b ∈ F_p. Для корректного группового закона
    требуется невырожденность, то есть дискриминант
        Δ = -16(4 a^3 + 27 b^2)
    не равен 0 по модулю p. Это гарантирует, что множество точек E(F_p)
    образует абелеву группу по сложению.
    """

    def __init__(self, a: int, b: int, p: int) -> None:
        """
        Инициализирует эллиптическую кривую y^2 ≡ x^3 + a x + b (mod p).
        Проверяет невырожденность кривой по условию:
            4 a^3 + 27 b^2 ≠ 0 (mod p),
        что эквивалентно Δ ≠ 0 для дискриминанта эллиптической кривой.
        При нарушении этого условия кривая имеет особые точки, и
        групповой закон на E(F_p) перестаёт быть корректным.
        Args:
            a: Коэффициент a в уравнении кривой.
            b: Коэффициент b в уравнении кривой.
            p: Модуль простого поля F_p.
        Raises:
            ValueError: Если 4 a^3 + 27 b^2 ≡ 0 (mod p), то есть кривая вырождена.
        """
        self.a = a
        self.b = b
        self.p = p

        discriminant = (4 * pow(a, 3, p) + 27 * pow(b, 2, p)) % p
        if discriminant == 0:
            raise ValueError("Elliptic curve is singular")

    def is_point_on_curve(self, point: Point) -> bool:
        """
        Проверяет, лежит ли точка на эллиптической кривой E(F_p).
        Для конечной точки (x, y) проверяется уравнение:
            y^2 ≡ x^3 + a x + b (mod p).
        Точка на бесконечности (нейтральный элемент группы) по определению
        считается принадлежащей кривой.
        Args:
            point: Точка, заданная координатами (x, y) или как точка на бесконечности.
        Returns:
            True, если point удовлетворяет уравнению кривой или является точкой
            на бесконечности, иначе False.
        """
        if point.is_infinity:
            return True

        assert point.x is not None and point.y is not None
        left = pow(point.y, 2, self.p)
        right = (pow(point.x, 3, self.p) + self.a * point.x + self.b) % self.p
        return left == right

    def add_points(self, p1: Point, p2: Point) -> Point:
        """
        Складывает две точки на эллиптической кривой по групповому закону.
        Реализует стандартный групповой закон для кривой y^2 = x^3 + a x + b (mod p):
        1) Нейтральный элемент:
           - если p1 = O, то p1 + p2 = p2;
           - если p2 = O, то p1 + p2 = p1.
        2) Противоположные точки:
           если x1 = x2 и y1 ≡ -y2 (mod p), то p1 + p2 = O.
        3) Удвоение точки P = (x1, y1) (случай p1 = p2):
           λ = (3 x1^2 + a) / (2 y1) (mod p),
           x3 = λ^2 - x1 - x2 (mod p),
           y3 = λ (x1 - x3) - y1 (mod p).
        4) Сложение разных точек P = (x1, y1), Q = (x2, y2), P ≠ ±Q:
           λ = (y2 - y1) / (x2 - x1) (mod p),
           x3 = λ^2 - x1 - x2 (mod p),
           y3 = λ (x1 - x3) - y1 (mod p).
        Во всех случаях деление по модулю p реализуется умножением на
        мультипликативный обратный элемент в F_p.
        Args:
            p1: Первая точка на кривой (может быть точкой на бесконечности).
            p2: Вторая точка на кривой (может быть точкой на бесконечности).
        Returns:
            Точка p1 + p2 в группе E(F_p).
        """
        if p1.is_infinity:
            return p2
        if p2.is_infinity:
            return p1

        assert p1.x is not None and p1.y is not None
        assert p2.x is not None and p2.y is not None

        x1, y1 = p1.x, p1.y
        x2, y2 = p2.x, p2.y

        if x1 == x2 and (y1 + y2) % self.p == 0:
            return Point.infinity()

        if p1 == p2:
            numerator = (3 * x1 * x1 + self.a) % self.p
            denominator = (2 * y1) % self.p
        else:
            numerator = (y2 - y1) % self.p
            denominator = (x2 - x1) % self.p

        lam = (numerator * pow(denominator, -1, self.p)) % self.p
        x3 = (lam * lam - x1 - x2) % self.p
        y3 = (lam * (x1 - x3) - y1) % self.p
        return Point(x3, y3)

    def multiply_point(self, point: Point, k: int) -> Point:
        """
        Вычисляет kP для точки P методом двоичного "double-and-add".
        Реализует скалярное умножение в циклической подгруппе:
            k P = P + P + ... + P  (k раз),
        но вместо наивного O(k) алгоритма используется двоичное разложение:
            k = ∑_{i} k_i 2^i,  k_i ∈ {0, 1},
        и схема:
            - если k_i = 1, добавить текущую точку tmp к результату;
            - на каждом шаге удвоить tmp (tmp ← 2 tmp).
        В результате сложность по времени — O(log k) операций сложения точек
        на кривой.
        Args:
            point: Точка P на кривой (элемент группы E(F_p)).
            k: Целое скалярное множитель (может быть большим).
        Returns:
            Точка kP. Для k = 0 или P = O возвращается точка на бесконечности.
        """
        if point.is_infinity or k == 0:
            return Point.infinity()

        result = Point.infinity()
        tmp = point

        while k > 0:
            if k & 1:
                result = self.add_points(result, tmp)
            tmp = self.add_points(tmp, tmp)
            k >>= 1

        return result

    def tonelli_shanks(self, n: int) -> int | None:
        """
        Находит квадратный корень r числа n по модулю p методом Tonelli–Shanks.
        Ищется r ∈ F_p такое, что:
            r^2 ≡ n (mod p).
        Алгоритм:
        1) Проверка существования корня по критерию Эйлера:
               n^{(p-1)/2} ≡ 1 (mod p)
           является необходимым и достаточным условием того, что n —
           квадратичный вычет (существует квадратный корень по модулю p).
        2) Разложение p - 1 = q · 2^s, q нечётно.
        3) Специальный случай s = 1 (p ≡ 3 mod 4):
               r ≡ n^{(p+1)/4} (mod p).
        4) Общий случай: используется элемент z, являющийся квадратичным
           невычетом, и итерационная схема Tonelli–Shanks, уменьшающая
           порядок элемента t, пока t не станет 1; тогда r — искомый корень.
        Args:
            n: Остаток по модулю p, для которого ищется квадратный корень.
        Returns:
            Любой квадратный корень r по модулю p (r^2 ≡ n mod p),
            либо None, если n не является квадратичным вычетом.
        """
        p = self.p
        n %= p
        if n == 0:
            return 0
        if pow(n, (p - 1) // 2, p) != 1:
            return None

        q = p - 1
        s = 0
        while q % 2 == 0:
            q //= 2
            s += 1

        if s == 1:
            return pow(n, (p + 1) // 4, p)

        z = 2
        while pow(z, (p - 1) // 2, p) == 1:
            z += 1

        c = pow(z, q, p)
        r = pow(n, (q + 1) // 2, p)
        t = pow(n, q, p)
        m = s

        while t != 1:
            t2i = t
            i = 0
            while i < m and t2i != 1:
                t2i = (t2i * t2i) % p
                i += 1

            b = pow(c, 1 << (m - i - 1), p)
            r = (r * b) % p
            c = (b * b) % p
            t = (t * c) % p
            m = i

        return r

    def get_random_point(self) -> Point:
        """
        Генерирует случайную точку на эллиптической кривой E(F_p).
        Алгоритм:
        1) Случайно выбирается x ∈ {0, 1, ..., p-1}.
        2) Вычисляется y^2 = x^3 + a x + b (mod p).
        3) Методом Tonelli–Shanks ищется квадратный корень y:
               y^2 ≡ x^3 + a x + b (mod p).
           Если корня нет (y = None), выбирается новый x.
        4) Проверяется принадлежность точки (x, y) кривой (формально избыточно,
           но даёт дополнительную защиту).
        В итоге возвращается точка P = (x, y) ∈ E(F_p), распределённая "равномерно"
        по выбранным x (при условии равномерности функции randrange).
        Returns:
            Случайная точка на кривой E(F_p).
        """
        while True:
            x = random.randrange(0, self.p)
            y_sq = (pow(x, 3, self.p) + self.a * x + self.b) % self.p
            y = self.tonelli_shanks(y_sq)
            if y is None:
                continue
            point = Point(x, y)
            if self.is_point_on_curve(point):
                return point

    def find_point_order(self, point: Point) -> int:
        """
        Вычисляет порядок точки P в группе E(F_p) методом baby-step/giant-step.
        Порядок точки P — это минимальное n ≥ 1, такое что:
            n P = O
        (O — точка на бесконечности). Порядок точки делит порядок группы E(F_p),
        а число точек на кривой удовлетворяет оценке Хассе:
            |#E(F_p) - (p + 1)| ≤ 2 √p.
        Алгоритм:
        1) Строится верхняя оценка order(P) через модифицированную границу Хассе:
               upper_bound ≈ p + 1 + 2 (⌊√p⌋ + 1).
        2) Выбирается m ≈ ⌊√upper_bound⌋, и используется baby-step/giant-step:
           - baby steps: записываются точки j P для j = 1..m в словарь;
           - giant steps: итеративно добавляется -m P к текущей точке и
             проверяется совпадение с одним из baby steps.
        3) Если найдены i, j такие что:
               j P = i m P,
           то (i m + j) P = O, и n = i m + j — порядок точки.
        Args:
            point: Точка P на кривой.
        Returns:
            Порядок точки P в группе E(F_p).
        Raises:
            ValueError: Если точка не лежит на кривой.
            RuntimeError: Если верхняя оценка оказалась недостаточной и
                порядок не найден в указанном диапазоне.
        """
        if point.is_infinity:
            return 1
        if not self.is_point_on_curve(point):
            raise ValueError("Point is not on curve")

        upper_bound = self.p + 1 + 2 * (isqrt(self.p) + 1)
        m = isqrt(upper_bound) + 1

        baby_steps: dict[Point, int] = {}
        current = point
        baby_steps[current] = 1

        for i in range(2, m + 1):
            current = self.add_points(current, point)
            baby_steps[current] = i

        assert point.x is not None and point.y is not None
        neg_point = Point(point.x, (-point.y) % self.p)
        m_times_neg_point = self.multiply_point(neg_point, m)
        current = Point.infinity()

        for i in range(0, m + 1):
            if current in baby_steps:
                j = baby_steps[current]
                return i * m + j
            current = self.add_points(current, m_times_neg_point)

        raise RuntimeError("Point order upper bound was too small")
