from __future__ import annotations

from math import gcd
from typing import Optional

from .models import NumberModulo


def normalize_mod(value: int, modulus: int) -> int:
    """
    Приводит целое число к представителю класса вычетов по модулю m.
    Реализует операцию:
        value (mod modulus),
    то есть возвращает остаток от деления value на modulus. В контексте
    арифметики на эллиптической кривой все вычисления ведутся в поле F_p
    или в кольце Z_n, поэтому такая нормализация используется для явного
    перехода к классу вычетов.
    Args:
        value: Целое число для приведения.
        modulus: Модуль, по которому берётся остаток.
    Returns:
        Остаток value % modulus в диапазоне [0, modulus).
    """
    return value % modulus


def mod_inverse(number: int, modulus: int) -> int:
    """
    Вычисляет мультипликативный обратный элемент по модулю m.
    Ищет такое x, что:
        number * x ≡ 1 (mod modulus),
    то есть x = number^{-1} в кольце Z_modulus. Обратимый элемент существует
    тогда и только тогда, когда gcd(number, modulus) = 1.
    В алгоритме Pollard rho это используется для решения линейных сравнений
    вида:
        a x ≡ b (mod n),
    а также при реализации группового закона на эллиптической кривой
    (деление заменяется умножением на обратный элемент в F_p).
    Args:
        number: Число, для которого ищется обратный элемент.
        modulus: Модуль, по которому считается обратимость.
    Returns:
        Мультипликативный обратный элемент number^{-1} (mod modulus).
    Raises:
        ArithmeticError: Если gcd(number, modulus) ≠ 1, и обратный элемент
            по модулю modulus не существует.
    """
    number %= modulus
    if gcd(number, modulus) != 1:
        raise ArithmeticError("Modular inverse does not exist")
    return pow(number, -1, modulus)


def solve_linear_congruence(
    number: int, x_multiplier: int, modulus: int
) -> Optional[NumberModulo]:
    """
    Решает линейное сравнение вида a x ≡ b (mod m).
    Формально решается уравнение:
        x_multiplier * x ≡ number (mod modulus).
    Теория:
        Пусть d = gcd(a, m). Тогда сравнение
            a x ≡ b (mod m)
        имеет решение тогда и только тогда, когда d | b. В этом случае
        все решения образуют класс вычетов:
            x ≡ x_0 (mod m / d),
        где x_0 — одно из решений сокращённого сравнения:
            (a / d) x ≡ (b / d) (mod m / d).
    Алгоритм функции:
    1) Нормализует коэффициенты a = x_multiplier и b = number по модулю m.
    2) Вычисляет d = gcd(a, m) и проверяет условие существования решения
       b ≡ 0 (mod d).
    3) Делит a, b и m на d, получая взаимно простые a' и m'.
    4) Находит обратный элемент a'^{-1} (mod m') и вычисляет
           x_0 ≡ b' * a'^{-1} (mod m'),
       где b' = b / d.
    5) Возвращает решение в виде NumberModulo(x_0, m').
    Args:
        number: b в сравнении a x ≡ b (mod m).
        x_multiplier: a в сравнении a x ≡ b (mod m).
        modulus: m — модуль сравнения.
    Returns:
        NumberModulo с одним из решений x ≡ x_0 (mod m'), где m' = modulus / d,
        либо None, если сравнение не имеет решений (d ∤ b).
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
