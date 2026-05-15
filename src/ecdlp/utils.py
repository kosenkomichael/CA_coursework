from __future__ import annotations

from math import gcd
from typing import Optional

from .models import NumberModulo


def normalize_mod(value: int, modulus: int) -> int:
    """Взятие числа по модулю

    Args:
        value (int): число
        modulus (int): модуль

    Returns:
        int: число по модулю
    """

    return value % modulus


def mod_inverse(number: int, modulus: int) -> int:
    """Обратный элемент по модулю p

    Args:
        number (int): число
        modulus (int): модуль

    Raises:
        ArithmeticError: обратного может не быть

    Returns:
        int: обратное к переданному число
    """

    number %= modulus
    if gcd(number, modulus) != 1:
        raise ArithmeticError("Modular inverse does not exist")
    return pow(number, -1, modulus)


def solve_linear_congruence(
    number: int, x_multiplier: int, modulus: int
) -> Optional[NumberModulo]:
    """Решение линейного сравнения

    Args:
        number (int):
        x_multiplier (int): x_multiplier * x = number mod p <- ищем x
        modulus (int):

    Returns:
        Optional[NumberModulo]: найденное решение
    """

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
