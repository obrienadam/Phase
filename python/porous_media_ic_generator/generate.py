import math

import matplotlib.pyplot as plt

from shapes import *
from solver import *

if __name__ == '__main__':
    lx = 2
    ly = 1
    n = 32
    vf_target = 0.5

    domain = Box(np.array([0, 0]), np.array([lx, ly]))
    r = np.random.normal(0.1, 0.02, n)

    nx = math.ceil(math.sqrt(lx / ly * n))
    ny = math.ceil(math.sqrt(ly / lx * n))

    x = np.linspace(lx / nx / 2, lx - lx / nx / 2, nx)
    y = np.linspace(ly / ny / 2, ly - ly / ny / 2, ny)
    x, y = np.meshgrid(x, y, indexing='ij')
    x = x.flatten()[:n]
    y = y.flatten()[:n]

    cylinders = [Cylinder(r, np.array([x, y])) for r, x, y in zip(r, x, y)]

    vf = np.sum([c.area() for c in cylinders]) / (lx * ly)
    ratio = (vf_target / vf) ** 0.5

    for c in cylinders:
        c.r *= ratio

    solver = Solver(cylinders, domain)

    fig, ax = plt.subplots()

    for c in cylinders:
        c = plt.Circle((c.x[0], c.x[1]), c.r, color='red')
        ax.add_artist(c)

    plt.xlim(0, lx)
    plt.ylim(0, ly)
    plt.show()

    r = solver.solve()

    vx = r[:, :n]
    vy = r[:, n:2 * n]
    x = r[:, 2 * n:3 * n]
    y = r[:, 3 * n:]

    fig, ax = plt.subplots()

    for c in cylinders:
        c = plt.Circle((c.x[0], c.x[1]), c.r, color='blue')
        ax.add_artist(c)

    plt.xlim(0, lx)
    plt.ylim(0, ly)
    plt.show()

    print(r.shape)
