"""Module that contains the most basic geometric shapes:
## 2D:
- Circle
- Rectangle
## 3D:
- Cylinder
- Sphere
- Cube
- Cone"""

import math

# 2-D Shapes #


class Circle:
    def __init__(self, diameter):
        self.Diameter = diameter
    @property
    def Radius(self):
        return self.Diameter/2
    @property
    def Area(self):
        return math.pi * self.Radius ** 2


class Rectangle:
    def __init__(self, height, length):
        self.Height = height
        self.Length = length
    @property
    def Diagonal(self):
        return math.sqrt(self.Height ** 2 + self.Length ** 2)
    @property
    def Area(self):
        return self.Height * self.Length


#  3-D Shapes #


class Cylinder:
    def __init__(self, diameter, height):
        self.Diameter = diameter
        self.Height = height
    @property
    def Radius(self):
        return self.Diameter/2
    @property
    def BaseArea(self):
        return math.pi * self.Radius ** 2
    @property
    def SideArea(self):
        return math.pi*self.Diameter*self.Height
    @property
    def Area(self):
        return 2 * self.BaseArea+self.SideArea
    @property
    def Volume(self):
        return self.BaseArea*self.Height


class Sphere:
    def __init__(self, diameter):
        self.Diameter = diameter
    @property
    def Radius(self):
        return self.Diameter/2
    @property
    def Area(self):
        return 4*math.pi*self.Radius**2
    @property
    def Volume(self):
        return math.pi*4/3*self.Radius**3


class Cube:
    def __init__(self, height, length, width):
        self.Height = height
        self.Length = length
        self.Width = width
    @property
    def BaseArea(self):
        return Rectangle(self.Height, self.Length).Area
    @property
    def Area(self):
        return 2*(self.Height * self.Length + self.Length * self.Width + self.Height * self.Width)
    @property
    def Volume(self):
        return self.Height * self.Length * self.Width


class Cone:
    def __init__(self, diameter=None, height=None, angle=None):
        self._D = diameter
        self._H = height
        self._angle = angle

    @property
    def Diameter(self):
        if self._D is None:
            return 2*self.Height/math.tan(math.radians(self.Angle))
        else:
            return self._D

    @property
    def Height(self):
        if self._H is None:
            return self.Radius * math.tan(math.radians(self.Angle))
        else:
            return self._H

    @property
    def Angle(self):
        if self._angle is None:
            return math.degrees(math.atan(self.Height / self.Radius))
        else:
            return self._angle

    @property
    def Radius(self):
        return self.Diameter / 2

    @property
    def Diagonal(self):
        return math.sqrt(self.Radius**2 + self.Height**2)

    @property
    def BaseArea(self):
        return math.pi * (self.Radius**2)

    @property
    def SideArea(self):
        return math.pi * self.Radius * self.Diagonal

    @property
    def Area(self):
        return self.BaseArea + self.SideArea

    @property
    def Volume(self):
        return (1.0/3.0) * self.BaseArea * self.Height
