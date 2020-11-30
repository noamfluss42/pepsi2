import numpy as np

v_start = 200


def kassam_in_vaccum(dt, x0, y0, v0, theta0):
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
        v_x = np.append(v_x, [v_x[0]])
        y = np.append(y, [dt * v_y[i] + y[i]])
        v_y = np.append(v_y, v_y[i] - 10 * dt)
        i += 1
    return x, y