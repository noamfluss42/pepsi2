# imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

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
    while i < 20 or y[-1] > 0.01:
        x = np.append(x, [dt * v_x[i] + x[i]])
        v_x = np.append(v_x, [v_x[i]])
        if friction_coefficient is not None:
            v_x[-1] -= friction_coefficient * ((v_x[i] ** 2 + v_y[i] ** 2) ** 0.5) * v_x[i] * dt
        y = np.append(y, [dt * v_y[i] + y[i]])
        v_y = np.append(v_y, v_y[i] - 10 * dt)

        if friction_coefficient is not None:
            v_y[-1] -= friction_coefficient * ((v_x[i] ** 2 + v_y[i] ** 2) ** 0.5) * v_y[i] * dt
        i += 1
    return x, y


def kassam_in_vaccum(dt, x0, y0, v0, theta0):
    return shooting_kassam(dt, x0, y0, v0, theta0)


def kassam_in_air(dt, x0, y0, v0, theta0, friction_coefficient):
    return shooting_kassam(dt, x0, y0, v0, theta0, friction_coefficient=friction_coefficient)


def q2_section6():
    x, y = kassam_in_vaccum(0.001, 0, 0, v_start, 50)
    if x[0] > 0:
        x_from_start = np.concatenate([[0], x])
    else:
        x_from_start = x
    plt.plot(x_from_start, y)
    plt.show()


def q3_section4():
    x, y = kassam_in_air(0.001, 0, 0, v_start, 50, friction_coefficient=C * A * RHO / M)
    if x[0] > 0:
        x_from_start = np.concatenate([[0], x])
    else:
        x_from_start = x
    plt.plot(x_from_start, y)
    plt.show()
q2_section6()
q3_section4()