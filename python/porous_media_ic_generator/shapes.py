import math
import numpy as np

class Cylinder:
    def __init__(self, r=0.5, x=np.zeros(2)):
        self.r = r
        self.x = x

    def area(self):
        return math.pi * self.r ** 2

class Box:
    def __init__(self, lc=np.zeros(2), uc=np.zeros(2)):
        self.lc = lc
        self.uc = uc