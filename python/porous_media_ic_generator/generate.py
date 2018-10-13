import math

import matplotlib.pyplot as plt

from shapes import *
from solver import *

if __name__ == '__main__':
    lx = 2
    ly = 1
    n = 72
    vf_target = 0.4
    offsets = (0.5, 0)

    domain = Box(Point2D(0, 0), Point2D(lx, ly))
    r = np.random.normal(0.1, 0.02, n)

    nx = math.ceil(math.sqrt(lx / ly * n))
    ny = math.ceil(math.sqrt(ly / lx * n))

    print(nx, ny)

    x = np.linspace(lx / nx / 2, lx - lx / nx / 2, nx)
    y = np.linspace(ly / ny / 2, ly - ly / ny / 2, ny)
    x, y = np.meshgrid(x, y, indexing='ij')
    x = x.flatten()[:n]
    y = y.flatten()[:n]

    cylinders = [Cylinder(r, Point2D(x, y)) for r, x, y in zip(r, x, y)]

    vf = np.sum([c.area() for c in cylinders]) / (lx * ly)
    ratio = (vf_target / vf) ** 0.5

    for c in cylinders:
        c.r *= ratio

    solver = Solver(cylinders, domain, eps=2.5e-3, damping=0.01, s=0.05)

    fig, ax = plt.subplots()

    for c in cylinders:
        c = plt.Circle((c.xc.x, c.xc.y), c.r, color='red')
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
        c = plt.Circle((c.xc.x, c.xc.y), c.r, color='blue')
        ax.add_artist(c)

    plt.xlim(0, lx)
    plt.ylim(0, ly)
    plt.show()

    sum = 0.
    for c in cylinders:
        sum += 2. * c.r

    print('mean diameter =', sum / len(cylinders))

    with open('cylinders.info', 'w') as f:
        for i, c in enumerate(cylinders):
            if c.xc.x > 0 and c.xc.x < lx and c.xc.y > 0 and c.xc.y < ly:
                f.write(
                    'Cylinder{}\n'.format(i) +
                    '{\n' +
                    '  geometry\n' +
                    '  {\n' +
                    '    center ({},{})\n'.format(c.xc.x + offsets[0], c.xc.y + offsets[1]) +
                    '    radius {}\n'.format(c.r) +
                    '  }\n'
                    '}\n\n'
                )
            else:
                print('Cylinder out of bounds, not writing')