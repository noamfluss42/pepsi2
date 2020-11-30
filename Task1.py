# imports
import numpy as np
import matplotlib.pyplot as plt

v_start = 200
C = 0.75
M = 32
A = 0.018
RHO = 1.2


def shooting_kassam(dt, x0, y0, v0, theta0, friction_coefficient=None):
    x = np.zeros(1)
    v_x = np.zeros(1)
    y = np.zeros(1)
    v_y = np.zeros(1)
    x[0] = x0
    v_x[0] = v0 * np.cos(theta0 * np.pi / 180)
    y[0] = y0
    v_y[0] = v0 * np.sin(theta0 * np.pi / 180)
    i = 0
    while i < 100 or y[-1] > 0.01:
        x = np.append(x, [dt * v_x[i] + x[i]])
        v_x = np.append(v_x, [v_x[i]])
        if friction_coefficient is not None:
            v_x[-1] -= friction_coefficient * ((v_x[i] ** 2 + v_y[i] ** 2) ** 0.5) * v_x[i] * dt
        y = np.append(y, [dt * v_y[i] + y[i]])
        v_y = np.append(v_y, v_y[i] - 10 * dt)

        if friction_coefficient is not None:
            v_y[-1] -= friction_coefficient * ((v_x[i] ** 2 + v_y[i] ** 2) ** 0.5) * v_y[i] * dt
        i += 1
    print(f"i is {i}")
    return x, y


def kassam_in_vaccum(dt, x0, y0, v0, theta0):
    return shooting_kassam(dt, x0, y0, v0, theta0)


def kassam_in_air(dt, x0, y0, v0, theta0, friction_coefficient):
    return shooting_kassam(dt, x0, y0, v0, theta0, friction_coefficient=friction_coefficient)


def q2_section6():
    x, y = kassam_in_vaccum(0.001, 0, 0, v_start, 10)
    if x[0] > 0:
        x_from_start = np.concatenate([[0], x])
    else:
        x_from_start = x
    plt.plot(x_from_start, y, label="theta 10")
    x, y = kassam_in_vaccum(0.001, 0, 0, v_start, 30)
    if x[0] > 0:
        x_from_start = np.concatenate([[0], x])
    else:
        x_from_start = x
    plt.plot(x_from_start, y, label="theta 30")
    x, y = kassam_in_vaccum(0.001, 0, 0, v_start / 2, 80)
    if x[0] > 0:
        x_from_start = np.concatenate([[0], x])
    else:
        x_from_start = x
    plt.plot(x_from_start, y, label="theta 89")
    plt.legend()
    plt.show()


def q3_section4():
    x, y = kassam_in_air(0.001, 0, 0, v_start, 50, friction_coefficient=C * A * RHO / M)
    if x[0] > 0:
        x_from_start = np.concatenate([[0], x])
    else:
        x_from_start = x
    plt.plot(x_from_start, y)

    x, y = kassam_in_air(0.001, 0, 0, v_start, 50, friction_coefficient=0)
    if x[0] > 0:
        x_from_start = np.concatenate([[0], x])
    else:
        x_from_start = x
    plt.plot(x_from_start, y, "")

    plt.show()


def x_hit(theta):
    x, y = kassam_in_air(0.001, 0, 0, v_start, theta, friction_coefficient=C * A * RHO / M)
    return x[-1]


def find_theta(x_dest, initial_theta1, initial_theta2):
    guesses = [initial_theta1, initial_theta2]
    consecutive_successes = 0

    def f(theta):
        return x_dest - x_hit(theta)
    f_prev = f(initial_theta1)
    while consecutive_successes < 3:
        curr = guesses[-1]
        prev = guesses[-2]
        f_curr = f(curr)
        nxt = curr - f_curr * (curr - prev) / (f_curr - f_prev)
        guesses.append(nxt)
        if abs(nxt-curr) < 0.01:
            consecutive_successes += 1
        else:
            consecutive_successes = 0

        f_prev = f_curr
    return guesses[-1]


def draw_x_hit():
    theta = np.linspace(5, 85, 100)
    x_hits = [x_hit(t) for t in theta]
    plt.plot(theta, x_hits)
    plt.xlabel("theta")
    plt.ylabel("x")
    plt.show()




#q2_section6()
#q3_section4()
#draw_x_hit()
print(f"theta is {find_theta(1750, 10, 30)}")
