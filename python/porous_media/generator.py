import scipy.stats as stats
import numpy as np

from point import Vector2D

def generate(void_fraction, num_cylinders, lx=2, ly=1):
    x = np.random.rand(num_cylinders) * lx
    y = np.random.rand(num_cylinders) * ly

    eps = 1e-4

    x = [Vector2D(x, y) for x, y in zip(x, y)]
    v = [Vector2D(0., 0.) for _ in num_cylinders]
    f = [Vector2D(0., 0.) for _ in num_cylinders]


print(generate(1, 10, 2, 1)[1])