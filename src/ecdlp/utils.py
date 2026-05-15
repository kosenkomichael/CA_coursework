from __future__ import annotations

from math import gcd
from typing import Optional

from .models import NumberModulo


def normalize_mod(value: int, modulus: int) -> int:

    return value % modulus


def mod_inverse(number: int, modulus: int) -> int:

    number %= modulus
    if gcd(number, modulus) != 1:
        raise ArithmeticError("Modular inverse does not exist")
    return pow(number, -1, modulus)


def solve_linear_congruence(
    number: int, x_multiplier: int, modulus: int
) -> Optional[NumberModulo]:

    x_multiplier %= modulus
    number %= modulus

    d = gcd(x_multiplier, modulus)
    if number % d != 0:
        return None

    number //= d
    x_multiplier //= d
    modulus //= d

    inverse = mod_inverse(x_multiplier, modulus)
    x = (number * inverse) % modulus
    return NumberModulo(x, modulus)
