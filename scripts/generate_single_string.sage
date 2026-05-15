from sage.all import GF, EllipticCurve, randint, random_prime
import sys


def generate_case(bits: int) -> str:
    if bits < 2:
        raise ValueError("bits must be >= 2")

    while True:
        try:
            p = random_prime(2**bits - 1, lbound=2**(bits - 1))
            F = GF(p)

            while True:
                a = F(randint(0, p - 1))
                b = F(randint(0, p - 1))
                if 4 * a**3 + 27 * b**2 != 0:
                    break

            E = EllipticCurve(F, [a, b])

            while True:
                P = E.random_point()
                if P != E(0):
                    break
            Px, Py = P.xy()

            while True:
                k = randint(1, p - 1)
                Q = k * P
                if Q != E(0):
                    break
            Qx, Qy = Q.xy()

            answer = Q.log(P)

            return (
                f"{int(p)} {int(a)} {int(b)} "
                f"{int(Px)} {int(Py)} "
                f"{int(Qx)} {int(Qy)} "
                f"{int(answer)}"
            )
        except Exception:
            pass


print(generate_case(30))