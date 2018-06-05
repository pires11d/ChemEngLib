
class Substance:
    def __init__(self, name):
        self.Name = str(name)
        self.Temperature = 298.15
        self.Pressure = 101325.0
        self.Phase = None
        self.MolarMass = None
        self.Density0 = None
        self.Viscosity0 = None
        self.SpecificHeat0 = None
        self.Diffusivity0 = None

    @property
    def Density(self):
        return self.Density0
    @property
    def Viscosity(self):
        return self.Viscosity0
    @property
    def SpecificHeat(self):
        return self.SpecificHeat0
    @property
    def Diffusivity(self):
        return self.Diffusivity0

