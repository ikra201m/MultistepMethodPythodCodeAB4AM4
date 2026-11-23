def f(t, y):
    return y - t**2 + 1


def adams_moulton_rk4():
    t = [0.0] * 1001
    y = [0.0] * 1001

    a = 0.0
    b = 2.0
    t[0] = a
    y[0] = 0.5
    h = 0.25
    n = round((b - a) / h)

    # first 3 by using RK4
    for i in range(1, 4):
        t[i] = a + i * h
        k1 = h * f(t[i - 1], y[i - 1])
        k2 = h * f(t[i - 1] + h / 2.0, y[i - 1] + k1 / 2.0)
        k3 = h * f(t[i - 1] + h / 2.0, y[i - 1] + k2 / 2.0)
        k4 = h * f(t[i - 1] + h, y[i - 1] + k3)

        y[i] = y[i - 1] + (1.0 / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        print(f"{i:2d} {t[i]:12.6f} {y[i]:12.6f}")

    #  3-step (4th-order) Adams-Moulton
    for i in range(3, n):
        t[i + 1] = a + (i + 1) * h
        y[i + 1] = y[i]  # (initial guess)

        #  (corrector step)
        for _ in range(3):
            y[i + 1] = y[i] + (h / 24.0) * (
                9.0 * f(t[i + 1], y[i + 1])
                + 19.0 * f(t[i], y[i])
                - 5.0 * f(t[i - 1], y[i - 1])
                + f(t[i - 2], y[i - 2])
            )

        print(f"{i + 1:2d} {t[i + 1]:12.6f} {y[i + 1]:12.6f}")


if __name__ == "__main__":
    adams_moulton_rk4()