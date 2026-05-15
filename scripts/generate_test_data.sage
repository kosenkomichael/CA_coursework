if __name__ == "__main__":
    with open("data/ro_alg_ecdlp_data.txt", "w", encoding="utf-8") as f:
        for bits in range(10, 31):
            for _ in range(100):
                while True:
                    try:
                        p = random_prime(2**bits, lbound=2**(bits - 1))
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

                        f.write(
                            f"{int(p)} {int(a)} {int(b)} "
                            f"{int(Px)} {int(Py)} "
                            f"{int(Qx)} {int(Qy)} "
                            f"{int(answer)}\n"
                        )
                        break
                    except Exception:
                        pass