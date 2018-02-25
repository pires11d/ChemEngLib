"""Module that contains classes that represent the chemical components of a system."""


from collections import OrderedDict
from Geometry import Circle


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


class Mixture:
    def __init__(self, components):
        # Mixture attributes #
        self.Components = components
        self.ComponentNames = [c.Name for c in self.Components]
        self.wi = None
        self.zi = None
        self.V = None
        self.M = None
        self.N = None
        self.Temperature = 298.15
        self.Pressure = 101325.0
        self.ParticleSize = None
        self.ParticleCharge = None
        # Other attributes #
        self._To = 273.15
        self._Po = 101325.0
        self.isElectrolyte = False

    @property
    def Color(self):
        Tmin = 273.15
        Tmax = 373.15
        T = self.Temperature
        if T < Tmin:
            T = Tmin
        if T > Tmax:
            T = Tmax
        if len(self.MassFractions) > 1:
            w0 = self.MassFractions[0]
            color = (w0*1.0, w0*1.0, w0*0.5, w0)
        else:
            r = (T - Tmin)/(Tmax-Tmin)
            color = [r*1.0, 0.0, (1-r)*1.0, 0.5]
        return color
    
    @property
    def Phase(self):
        if self.GasContent > 0:
            if self.LiquidContent == 0:
                if self.SolidContent == 0:
                    return 'g'
                else:
                    return None
            else:
                return None
        if self.LiquidContent > 0:
            if self.GasContent == 0:
                if self.SolidContent == 0:
                    return 'l'
                else:
                    return None
            else:
                return None    
        if self.SolidContent > 0:
            if self.GasContent == 0:
                if self.LiquidContent == 0:
                    return 's'
                else:
                    return None
            else:
                return None    

    @property
    def MolarMass(self):
        MM = 0
        # TODO: Rever
        for i,MMi in enumerate(self.MolarMasses.values()):
            MM += self.MolarFractions[i] * MMi
        return MM
    @property
    def MolarMasses(self):
        MMi = []
        for c in self.Components:
            MMi.append(c.MolarMass)
        return OrderedDict(zip(self.ComponentNames, MMi))
    
    @property
    def MassFractions(self):
        wi = []
        if self.wi is None:
            for i, MMi in enumerate(self.MolarMasses.values()):
                den = 0
                num = self.MolarFractions[i] * MMi
                for j, MMj in enumerate(self.MolarMasses.values()):
                    den += self.MolarFractions[j] * MMj
                wi.append(num / den)
        else:
            wi = self.wi
        return wi
    @property
    def MolarFractions(self):
        zi = []
        if self.zi is None:
                for i, MMi in enumerate(self.MolarMasses.values()):
                    den = 0
                    num = self.MassFractions[i] / MMi
                    for j, MMj in enumerate(self.MolarMasses.values()):
                        den += self.MassFractions[j] / MMj
                    zi.append(num / den)
        else:
            zi = self.zi
        return zi

    @property
    def Mass(self):
        if self.M is None:
            if self.N is None:
                M = self.V * self.Density
            else:
                M = self.N * self.MolarMass
        else:
            M = self.M
        return M
    @property
    def Volume(self):
        if self.V is None:
            if self.N is None:
                V = self.M / self.Density
            else:
                # Ideal gas #
                V = self.N * 8.314 * self.Temperature / self.Pressure 
        else:
            V = self.V
        return V
    @property
    def Moles(self):
        if self.N is None:
            if self.V is None:
                N = self.M / self.MolarMass
            else:
                # Ideal gas #
                N = (self.Pressure * self.V) / (8.314 * self.Temperature)
        else:
            N = self.N
        return N

    @property
    def Density(self):
        rho = 0
        for i,c in enumerate(self.Components):
            rho += self.MassFractions[i] / c.Density
        return 1 / rho
    @property
    def Viscosity(self):
        mu = 1.0
        for i,c in enumerate(self.Components):
            mu = mu*(c.Viscosity)**self.MassFractions[i]
        return mu
    @property
    def SpecificHeat(self):
        Cp = 0
        for i,c in enumerate(self.Components):
            Cp += self.MassFractions[i] * c.SpecificHeat
        return Cp

    # Phase properties #

    @property
    def GasComponents(self):
        gas_list = []
        for c in self.Components:
            if c.Phase == 'g':
                gas_list.append(c)
        return gas_list
    @property
    def LiquidComponents(self):
        liquid_list = []
        for c in self.Components:
            if c.Phase == 'l':
                liquid_list.append(c)
        return liquid_list
    @property
    def SolidComponents(self):
        solid_list = []
        for c in self.Components:
            if c.Phase == 's':
                solid_list.append(c)
        return solid_list

    @property
    def GasFractions(self):
        if self.GasComponents == []:
            return [0]
        gf = []
        for i,c in enumerate(self.Components):
            if c.Phase == 'g':
                gf.append(self.MassFractions[i])
        gf = [gfi / sum(gf) for gfi in gf]
        return gf
    @property
    def LiquidFractions(self):
        if self.LiquidComponents == []:
            return [0]
        lf = []
        for i,c in enumerate(self.Components):
            if c.Phase == 'l':
                lf.append(self.MassFractions[i])
        lf = [lfi / sum(lf) for lfi in lf]
        return lf
    @property
    def SolidFractions(self):
        if self.SolidComponents == []:
            return [0]
        sf = []
        for i,c in enumerate(self.Components):
            if c.Phase == 's':
                sf.append(self.MassFractions[i])
        sf = [sfi / sum(sf) for sfi in sf]
        return sf

    @property
    def GasContent(self):
        gc = 0
        for i,c in enumerate(self.Components):
            if c.Phase == 'g':
                gc += self.MassFractions[i]
        return gc
    @property
    def LiquidContent(self):
        lc = 0
        for i,c in enumerate(self.Components):
            if c.Phase == 'l':
                lc += self.MassFractions[i]
        return lc
    @property
    def SolidContent(self):
        sc = 0
        for i,c in enumerate(self.Components):
            if c.Phase == 's':
                sc += self.MassFractions[i]
        return sc

    @property
    def GasDensity(self):
        if self.GasContent > 0:
            rho = 0
            for i,c in enumerate(self.GasComponents):
                rho += self.GasFractions[i] / c.Density
            return 1 / rho
        else:
            return None
    @property
    def LiquidDensity(self):
        if self.LiquidContent > 0:
            rho = 0
            for i,c in enumerate(self.LiquidComponents):
                rho += self.LiquidFractions[i] / c.Density
            return 1 / rho    
        else:
            return None
    @property
    def SolidDensity(self):
        if self.SolidContent > 0:
            rho = 0
            for i,c in enumerate(self.SolidComponents):
                rho += self.SolidFractions[i] / c.Density
            return 1 / rho    
        else:
            return None        

    @property
    def GasViscosity(self):
        mu = 1.0
        for i,c in enumerate(self.GasComponents):
            mu = mu*(c.Viscosity)**self.GasFractions[i]
        return mu
    @property
    def LiquidViscosity(self):
        mu = 1.0
        for i,c in enumerate(self.LiquidComponents):
            mu = mu*(c.Viscosity)**self.LiquidFractions[i]
        return mu

    @property
    def GasSpecificHeat(self):
        Cp = 0
        for i,c in enumerate(self.GasComponents):
            Cp += self.GasFractions[i] * c.SpecificHeat
        return Cp
    @property
    def LiquidSpecificHeat(self):
        Cp = 0
        for i,c in enumerate(self.LiquidComponents):
            Cp += self.LiquidFractions[i] * c.SpecificHeat
        return Cp
    @property
    def SolidSpecificHeat(self):
        Cp = 0
        for i,c in enumerate(self.SolidComponents):
            Cp += self.SolidFractions[i] * c.SpecificHeat
        return Cp


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
        # Other attributes #
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
                return self.Mf / self.Density
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
