import numpy as np
import scipy.integrate as ode
from shapes import Point2D

class Solver:
    def __init__(self, cylinders, domain, damping=0.05, dt=1e-3, eps=1e-3, s=0.01):
        self.cylinders = cylinders
        self.domain = domain
        self.damping = damping
        self.dt = dt
        self.eps = eps
        self.s = s

        self.m = np.array([c.area() for c in self.cylinders])

    def compute_forces(self, vx, vy):
        fx = np.zeros(len(self.cylinders))
        fy = np.zeros(len(self.cylinders))

        for i, cp in enumerate(self.cylinders):
            f = Point2D(0., 0.)

            # cylinder-cylinder collisions
            for cq in self.cylinders:
                if cq is cp:
                    continue

                zeta = (cq.xc - cp.xc).mag() - (cp.r + cq.r)
                f += (cp.xc - cq.xc) / self.eps * max(0., -(zeta - self.s)) ** 2

            # cylinder-wall collisions
            lc = self.domain.lc
            uc = self.domain.uc

            for pt in Point2D(lc.x, cp.xc.y), Point2D(uc.x, cp.xc.y), Point2D(cp.xc.x, lc.y), Point2D(cp.xc.x, uc.y):
                zeta = (cp.xc - pt).mag() - cp.r
                f += (cp.xc - pt) / self.eps * max(0., -(zeta - self.s / 2.)) ** 2

            fx[i] = f.x
            fy[i] = f.y

        return fx - self.damping * vx, fy - self.damping * vy

    def solve(self):
        def f(y, t):
            n = y.shape[0] // 4

            vx = y[:n]
            vy = y[n:2 * n]
            x = y[2 * n:3 * n]
            y = y[3 * n:]

            for xc, yc, c in zip(x, y, self.cylinders):
                c.xc.x = xc
                c.xc.y = yc

            fx, fy = self.compute_forces(vx, vy)

            print('t =', t)

            return np.concatenate((fx / self.m, fy / self.m, vx, vy))

        vx = np.random.normal(0, 0.15, len(self.cylinders))
        vy = np.random.normal(0, 0.15, len(self.cylinders))
        x = np.array([c.xc.x for c in self.cylinders])
        y = np.array([c.xc.y for c in self.cylinders])

        return ode.odeint(f, np.concatenate((vx, vy, x, y)), np.arange(0, 15, self.dt))
