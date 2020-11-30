# imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

v_start = 200
C = 0.75
M = 32
A = 0.018
RHO = 1.2
add_block_size = 50000


def shooting_kassam(dt, x0, y0, v0, theta0, friction_coefficient=None):
    x = np.zeros(add_block_size)
    v_x = np.zeros(add_block_size)
    y = np.zeros(add_block_size)
    v_y = np.zeros(add_block_size)
    x[0] = x0
    v_x[0] = v0 * np.cos(theta0 * np.pi / 180)
    y[0] = y0
    v_y[0] = v0 * np.sin(theta0 * np.pi / 180)
    i = 0
    while i < 100 or y[i] > 0.01:
        if i + 1 >= x.shape[0]:
            x = np.append(x, np.zeros(add_block_size))
            v_x = np.append(v_x, np.zeros(add_block_size))
            y = np.append(y, np.zeros(add_block_size))
            v_y = np.append(v_y, np.zeros(add_block_size))
        x[i + 1] = dt * v_x[i] + x[i]
        v_x[i + 1] = v_x[i]
        if friction_coefficient is not None:
            v_x[i + 1] -= friction_coefficient * ((v_x[i] ** 2 + v_y[i] ** 2) ** 0.5) * v_x[i] * dt
        y[i + 1] = dt * v_y[i] + y[i]
        v_y[i + 1] = v_y[i] - 10 * dt

        if friction_coefficient is not None:
            v_y[i + 1] -= friction_coefficient * ((v_x[i] ** 2 + v_y[i] ** 2) ** 0.5) * v_y[i] * dt
        i += 1
    print(f"i is {i}")
    return x[:i + 1], y[:i + 1]


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


def q3_section6():
    log_dt_values = np.linspace(start=-4, stop=-1, num=500)
    dt_values = np.power(10, log_dt_values)
    result = [kassam_in_air(dt, 0, 0, v_start, 50, friction_coefficient=C * A * RHO / M)[0][-1] for dt in dt_values]
    plt.loglog(dt_values, result)
    plt.xlabel('$dt$', size=15)
    plt.ylabel(r'x landing', size=15)
    plt.grid()
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
        if abs(nxt - curr) < 0.01:
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


def find_minimal_distance ( x0_first , theta0_first , x0_second , theta0_second ):
    x_first, y_first = kassam_in_air(0.001, x0_first, 0, v_start, theta0_first, friction_coefficient=C * A * RHO / M)
    x_second, y_second = kassam_in_air(0.001, x0_second, 0, v_start, theta0_second, friction_coefficient=0.7*C * A * RHO / M)
    x_first = np.append(x_first, abs(len(x_first) - len(x_second))*[x_first[-1]])
    y_first = np.append(y_first, abs(len(y_first) - len(y_second))*[y_first[-1]])
    x_second = np.append(x_second, abs(len(x_first) - len(x_second))*[x_second[-1]])
    y_second = np.append(y_second, abs(len(y_first) - len(y_second))*[y_second[-1]])
    dx_power_2 = np.power(x_first - x_second,2)
    dy_power_2 = np.power(y_first - y_second,2)
    return np.amin(dx_power_2 + dy_power_2)


# q2_section6()
# q3_section4()
# draw_x_hit()
print(f"theta is {find_theta(1750, 10, 30)}")
