import math

class Point2D:
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y

    def __add__(self, other):
        return Vector2D(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return Vector2D(self.x - other.x, self.y - other.y)

class Vector2D(Point2D):
    def mag(self):
        return math.sqrt(self.x ** 2 + self.y ** 2)
