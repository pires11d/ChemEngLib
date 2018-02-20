"""Module that contains classes that represent the chemical components of a system."""


from collections import OrderedDict
from Geometry import Circle


class Substance:
    def __init__(self, name):
        self.Name = name
        self.Temperature = 298.15
        self.Pressure = 101325.0
        self.MolarMass = None


class Mixture:
    def __init__(self, components, wi=None, zi=None):
        self.Components = components
        self.MassFractions = wi
        self.Temperature = 298.15
        self.Pressure = 101325.0
        self.Density0 = None
        self.Viscosity0 = None
        self.Diffusivity0 = None
    
    @property
    def Density(temperature):
        


        