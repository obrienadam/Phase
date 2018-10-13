import math
import numpy as np

class Point2D:
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y

    def mag(self):
        return math.sqrt(self.x ** 2 + self.y ** 2)

    def unit(self):
        return self / self.mag()

    def __str__(self):
        return '({},{})'.format(self.x, self.y)

    def __add__(self, u):
        return Point2D(self.x + u.x, self.y + u.y)

    def __sub__(self, u):
        return Point2D(self.x - u.x, self.y - u.y)

    def __mul__(self, a):
        return Point2D(a * self.x, a * self.y)

    def __rmul__(self, a):
        return self.__mul__(a)

    def __truediv__(self, a):
        return Point2D(self.x / a, self.y / a)

class Cylinder:
    def __init__(self, r=0.5, xc=Point2D(0, 0)):
        self.r = r
        self.xc = xc

    def area(self):
        return math.pi * self.r ** 2

class Box:
    def __init__(self, lc=Point2D(0, 0), uc=Point2D(0, 0)):
        self.lc = lc
        self.uc = uc