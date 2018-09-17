import math

import numpy as np
from scipy.integrate import solve_ivp

import matplotlib.pyplot as plt


class Func:
    def __init__(self, lx, ly, nc, eps=1e-4, s=1e-2):
        self.lx = lx
        self.ly = ly
        self.nc = nc
        self.eps = eps
        self.s = s

    def f_wall(self, xc, yc, rc):
        return 0, 0

    def f_collision(self, xc_p, yc_p, rc_p, xc_q, yc_q, rc_q):
        return 0, 0

    def __call__(self, t, y):
        xc = y[:self.nc // 6]
        yc = y[self.nc // 6:self.nc // 3]
        vxc = y[self.nc // 3:self.nc // 2]
        vyc = y[self.nc // 2:2 * self.nc // 3]
        rc = y[2 * self.nc // 3:5 * self.nc // 6]
        dr = y[5 * self.nc // 6:]

        fxc = np.zeros_like(xc)
        fyc = np.zeros_like(yc)

        for i, xc_p, yc_p, rc_p in zip(range(self.nc), xc, yc, rc):
            fw = self.f_wall(xc, yc, rc)

            fxc[i] += fw[0]
            fyc[i] += fw[1]

            for j, xc_q, yc_q, rc_q in zip(range(self.nc, xc, yc, rc)):
                if i == j:
                    continue

                fc = self.f_collision(xc_p, yc_p, rc_p, xc_q, yc_q, rc_q)

                fxc[i] += fc[0]
                fyc[i] += fc[1]

        return vxc, vyc, fxc, fyc, dr, np.zeros_like(rc)


if __name__ == '__main__':
    lx = 2
    ly = 1
    nc = 36

    f = Func(lx, ly, nc)

    x = np.random.rand(nc) * lx
    y = np.random.rand(nc) * ly

    solve_ivp(f.__call__, (0, 5), np.zeros(6 * nc))
