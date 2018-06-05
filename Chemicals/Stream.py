from .Mixture import Mixture

class Stream(Mixture):
    def __init__(self, components):
        # Mixture attributes #
        self.Components = components
        self.ComponentNames = [c.Name for c in self.Components]
        self.wi = None
        self.zi = None
        self.Temperature = 298.15
        self.Pressure = 101325.0
        self.ParticleSize = None
        self.ParticleCharge = None
        # Stream attributes #
        self.Tag = None
        self.From = None
        self.To = None
        self.Mf = None
        self.Vf = None
        self.Vfo = None
        self.Nf = None
        # Equilibrium #
        self.Ki = None
        # Other attributes #
        self._phase = None
        self._To = 273.15
        self._Po = 101325.0
        self.isElectrolyte = False

    @property
    def MassFlow(self):
        if self.Mf is None:
            if self.Nf is None:
                return self.VolumeFlow * self.Density
            else:
                return self.Nf * self.MolarMass
        else:
            return self.Mf
    @property
    def VolumeFlow(self):
        if self.Vf is None:
            if self.Vfo is None:
                return self.MassFlow / self.Density
            else:
                return self.NormalVolumeFlow * (self._Po / self._To) * (self.Temperature / self.Pressure)
        else:
            return self.Vf
    @property
    def NormalVolumeFlow(self):
        if self.Vfo is None:
            return self.VolumeFlow * (self._To / self._Po) * (self.Pressure / self.Temperature)
        else:
            return self.Vfo
    @property
    def MolarFlow(self):
        if self.Nf is None:
            return self.MassFlow / self.MolarMass
        else:
            return self.Nf

    @property
    def MassFlows(self):
        return [wi * self.MassFlow for wi in self.MassFractions]
    @property
    def VolumeFlows(self):
        return [zi * self.VolumeFlow for zi in self.MolarFractions]
    @property
    def MolarFlows(self):
        return [zi * self.MolarFlow for zi in self.MolarFractions] 
