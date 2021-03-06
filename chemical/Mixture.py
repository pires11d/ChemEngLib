from collections import OrderedDict

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
        # Equilibrium #
        self.Ki = None
        # Other attributes #
        self._phase = None
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
            color = (w0, w0, (1-w0)/2, 0.7)
        else:
            r = (T - Tmin)/(Tmax-Tmin)
            color = [r, 1-r, 1-r, 0.5]
        return color
    
    @property
    def Hatch(self):
        hatch = None
        if self.Phase == 'g':
            hatch = '..'
            return hatch
        if self.SolidContent > 0:
            hatch = '.'
            if self.SolidContent > 0.2:
                hatch = '..'
                if self.SolidContent > 0.4:
                    hatch = '...'
                    if self.SolidContent > 0.6:
                        hatch = '....'
                        if self.SolidContent > 0.8:
                            hatch = '.....'
        return hatch

    @property
    def Phase(self):
        if self._phase != None:
            return self._phase
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
                if self.Phase == 'g':
                    V = self.N * 8.314 * self.Temperature / self.Pressure 
                else:
                    V = self.N * self.MolarMass / self.Density
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
                if self.Phase == 'g':
                    N = (self.Pressure * self.V) / (8.314 * self.Temperature)
                else:
                    N = self.V * self.Density / self.MolarMass
        else:
            N = self.N
        return N

    @property
    def ComponentMass(self):
        return [wi*self.Mass for wi in self.MassFractions]
    @property
    def ComponentMoles(self):
        return [zi*self.Moles for zi in self.MolarFractions]

    @property
    def MassConcentrations(self):
        return [mi/self.Volume for mi in self.ComponentMass]
    @property
    def MolarConcentrations(self):
        return [ni/self.Volume for ni in self.ComponentMoles]

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