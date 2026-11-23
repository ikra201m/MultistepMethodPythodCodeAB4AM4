def f(t, y):
    return y - t**2 + 1


def adams_bashforth_rk4():
    t = [0.0] * 1001
    y = [0.0] * 1001

    t[0] = 0.0
    y[0] = 0.5
    h = 0.25
    n = round(2.0 / h)

    print(f"{0:2d} {t[0]:12.6f} {y[0]:12.6f}")
# for the for three using rk
    
    for i in range(1, 4):
        t[i] = i * h
        k1 = h * f(t[i - 1], y[i - 1])
        k2 = h * f(t[i - 1] + h / 2.0, y[i - 1] + k1 / 2.0)
        k3 = h * f(t[i - 1] + h / 2.0, y[i - 1] + k2 / 3.0)
        k4 = h * f(t[i - 1] + h, y[i - 1] + k3)

        y[i] = y[i - 1] + (1.0 / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        print(f"{i:2d} {t[i]:12.6f} {y[i]:12.6f}")

    # for the next 3
    for i in range(3, n):
        t[i + 1] = (i + 1) * h
        y[i + 1] = y[i] + (h / 24.0) * (
            55.0 * f(t[i], y[i])
            - 59.0 * f(t[i - 1], y[i - 1])
            + 37.0 * f(t[i - 2], y[i - 2])
            - 9.0 * f(t[i - 3], y[i - 3])
        )
        print(f"{i + 1:2d} {t[i + 1]:12.6f} {y[i + 1]:12.6f}")


if __name__ == "__main__":
    adams_bashforth_rk4()