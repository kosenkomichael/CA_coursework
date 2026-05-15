from __future__ import annotations

import random
from math import isqrt

from .models import Point


class EllipticCurve:

    def __init__(self, a: int, b: int, p: int) -> None:

        self.a = a
        self.b = b
        self.p = p

        discriminant = (4 * pow(a, 3, p) + 27 * pow(b, 2, p)) % p
        if discriminant == 0:
            raise ValueError("Elliptic curve is singular")

    def is_point_on_curve(self, point: Point) -> bool:

        if point.is_infinity:
            return True

        assert point.x is not None and point.y is not None
        left = pow(point.y, 2, self.p)
        right = (pow(point.x, 3, self.p) + self.a * point.x + self.b) % self.p
        return left == right

    def add_points(self, p1: Point, p2: Point) -> Point:

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
