from collections import OrderedDict
from thermo import Chemical as che
from thermo import Mixture as mix
from thermo import electrochem as el
from thermo import UNIFAC
from Geometry import Circle
from Numerical import *


class Substance:
    def __init__(self, name, T=298.15, P=101325.0):
        self.Name = name
        self.Temperature = T
        self.Pressure = P

    @property
    def chemical(self):
        return che(self.Name, T=self.Temperature, P=self.Pressure)
    @property
    def CAS(self):
        return self.chemical.CAS
    @property
    def MolarMass(self):
        return self.chemical.MW
    @property
    def Phase(self):
        return self.chemical.phase
    @property
    def Density(self):
        return self.chemical.rho
    @property
    def HeatCapacity(self):
        return self.chemical.Cp
    @property
    def Viscosity(self):
        return self.chemical.mu


class Mixture:
    def __init__(self, components, wi, T=298.15, P=101325.0, electrolyte=False):
        self.Components = components
        self.MassFractions = wi
        self.Temperature = T
        self.Pressure = P
        self.isElectrolyte = electrolyte

    @property
    def mixture(self):
        return mix(self.Components, ws=self.MassFractions, T=self.Temperature, P=self.Pressure)
    @property
    def CAS(self):
        return self.mixture.CASs
    @property
    def MolarMass(self):
        return self.mixture.MW
    @property
    def Phase(self):
        return self.mixture.phase
    @property
    def Density(self):
        if self.isElectrolyte is False:
            return self.mixture.rho
        else:
            return el.Laliberte_density(self.Temperature, self._MassFractions, self._CAS)
    @property
    def SpecificHeat(self):
        if self.isElectrolyte is False:
            return self.mixture.Cp
        else:
            return el.Laliberte_heat_capacity(self.Temperature, self._MassFractions, self._CAS)
    @property
    def Viscosity(self):
        if self.isElectrolyte is False:
            return self.mixture.mu
        else:
            return el.Laliberte_viscosity(self.Temperature, self._MassFractions, self._CAS)
    @property
    def Enthalpy(self):
        return self.mixture.H

# HIDDEN PROPERTIES

    @property
    def _w(self):
        for i, c in enumerate(self.Components):
            if c == 'water':
                return i
    @property
    def _CAS(self):
        cas = list(self.CAS)
        if len(cas) == len(self.CAS):
            cas.remove(self.CAS[self._w])
        return cas
    @property
    def _MassFractions(self):
        mf = list(self.MassFractions)
        if len(mf) == len(self.MassFractions):
            mf.remove(self.MassFractions[self._w])
        return mf


class FlowStream:
    def __init__(self, components,
                 T=298.15, P=101325.0,
                 wi=None, vi=None, zi=None,
                 Mf=None, Vf=None, Vfo=None, Nf=None,
                 electrolyte=False):

        self.Components = components
        self.Temperature = T
        self.Pressure = P
        self.isElectrolyte = electrolyte

        self._wi = wi
        self._vi = vi
        self._zi = zi

        self._Mf = Mf
        self._Vf = Vf
        self._Vfo = Vfo
        self._Nf = Nf

        self._Po = 101325.0
        self._To = 273.15
        self._Fe2O3 = {'CASRN': '1309-37-1', 'name': 'hematite', 'MW': 159.688, 'formula': 'Fe2O3'}

        self.Tag = None
        self.From = None
        self.To = None
        self.ParticleSize = None
        self.ParticleCharge = None

    @property
    def Phase(self):
        if self._IronOxide is True:
            if self.SolidContent > 0:
                if self.SolidContent < 1:
                    return 'two-phase'
                else:
                    return 's'
            else:
                if self.GasContent == 1.0:
                    return 'g'
                elif self.LiquidContent == 1.0:
                    return 'l'
                elif self.SolidContent == 1.0:
                    return 's'
                else:
                    return 'two-phase'
        else:
            return mix(self._Components, ws=self.MassFractions,
                       T=self.Temperature, P=self.Pressure).phase
    @property
    def MolarMass(self):
        return mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).MW
    @property
    def Pbubble(self):
        return mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).Pbubble
    @property
    def Pdew(self):
        return mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).Pdew

    # TOTAL FLOW #

    @property
    def MassFlow(self):
        if self._Mf is None:
            if self._Nf is None:
                return self.VolumeFlow * self.Density
            else:
                return self._Nf * self.MolarMass
        else:
            return self._Mf
    @property
    def VolumeFlow(self):
        if self._Vf is None:
            if self._Vfo is None:
                return self._Mf / self.Density
            else:
                return self.NormalVolumeFlow * (self._Po / self._To) * (self.Temperature / self.Pressure)
        else:
            return self._Vf
    @property
    def NormalVolumeFlow(self):
        if self._Vfo is None:
            return self.VolumeFlow * (self._To / self._Po) * (self.Pressure / self.Temperature)
        else:
            return self._Vfo
    @property
    def MolarFlow(self):
        if self._Nf is None:
            return self.MassFlow / self.MolarMass
        else:
            return self._Nf

    # COMPONENT LISTS #

    @property
    def MassFractions(self):
        if self._wi is None:
            wi = []
            for i, MW in enumerate(self.MolarMasses.values()):
                den = 0
                num = self.MolarFractions[i] * MW
                for j, mw in enumerate(self.MolarMasses.values()):
                    den += self.MolarFractions[j] * mw
                wi.append(num / den)
            return wi
        else:
            return self._wi
    @property
    def MolarFractions(self):
        if self._zi is None:
            if self._vi is None:
                zi = []
                for i, MW in enumerate(self.MolarMasses.values()):
                    den = 0
                    num = self.MassFractions[i] / MW
                    for j, mw in enumerate(self.MolarMasses.values()):
                        den += self.MassFractions[j] / mw
                    zi.append(num / den)
                return zi
            else:
                return self._vi
        else:
            return self._zi
    @property
    def VolumeFractions(self):
        if self.Phase == 'g':
            if self._vi is None:
                return self.MolarFractions
            else:
                return self._vi
        else:
            return 'not a gas!'
    @property
    def ElementMassFractions(self):
        return mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).mass_fractions
    @property
    def ElementMolarFractions(self):
        return mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).atom_fractions

    @property
    def MassFlows(self):
        return [self.MassFlow * wi for wi in self.MassFractions]
    @property
    def MolarFlows(self):
        return [self.MolarFlow * zi for zi in self.MolarFractions]
    @property
    def VolumeFlows(self):
        if self.Phase == 'g':
            return [self.VolumeFlow * vi for vi in self.VolumeFractions]
        else:
            return 'not a gas!'
    @property
    def ElementMassFlows(self):
        values = [self.MassFlow * ei for ei in self.ElementMassFractions.values()]
        return OrderedDict(zip(self.Elements, values))
    @property
    def ElementMolarFlows(self):
        values = [self.MolarFlow * ei for ei in self.ElementMolarFractions.values()]
        return OrderedDict(zip(self.Elements, values))

    # PHASE LISTS #

    @property
    def GasComponents(self):
        gas_list = []
        for component in self._Components:
            if che(component, T=self.Temperature, P=self.Pressure).phase == 'g':
                gas_list.append(component)
        return gas_list
    @property
    def LiquidComponents(self):
        liquid_list = []
        for component in self._Components:
            if che(component, T=self.Temperature, P=self.Pressure).phase == 'l':
                liquid_list.append(component)
        return liquid_list
    @property
    def SolidComponents(self):
        solid_list = []
        for component in self._Components:
            if che(component, T=self.Temperature, P=self.Pressure).phase == 's':
                solid_list.append(component)
        return solid_list

    @property
    def GasFractions(self):
        gf = []
        for i, component in enumerate(self._Components):
            if che(component, T=self.Temperature, P=self.Pressure).phase == 'g':
                gf.append(self.MassFractions[i])
        gf = [gfi / sum(gf) for gfi in gf]
        return gf
    @property
    def LiquidFractions(self):
        lf = []
        for i, component in enumerate(self._Components):
            if che(component, T=self.Temperature, P=self.Pressure).phase == 'l':
                lf.append(self.MassFractions[i])
        lf = [lfi / sum(lf) for lfi in lf]
        return lf
    @property
    def SolidFractions(self):
        sf = []
        for i, component in enumerate(self._Components):
            if che(component, T=self.Temperature, P=self.Pressure).phase == 's':
                sf.append(self.MassFractions[i])
        sf = [sfi / sum(sf) for sfi in sf]
        return sf

    @property
    def GasContent(self):
        gc = 0
        for i, component in enumerate(self._Components):
            if che(component, T=self.Temperature, P=self.Pressure).phase == 'g':
                gc += self.MassFractions[i]
        return gc
    @property
    def LiquidContent(self):
        lc = 0
        for i, component in enumerate(self._Components):
            if che(component, T=self.Temperature, P=self.Pressure).phase == 'l':
                lc += self.MassFractions[i]
        return lc
    @property
    def SolidContent(self):
        sc = 0
        for i, component in enumerate(self._Components):
            if che(component, T=self.Temperature, P=self.Pressure).phase == 's' or che(component, T=self.Temperature, P=self.Pressure).phase is None:
                sc += self.MassFractions[i]
        return sc

    # PHYSICAL PROPERTIES #

    @property
    def GasDensity(self):
        if self.GasContent > 0:
            return mix(self.GasComponents, ws=self.GasFractions,
                       T=self.Temperature, P=self.Pressure).rhog
        else:
            return 'not a gas!'
    @property
    def LiquidDensity(self):
        if self.LiquidContent > 0:
            return mix(self.LiquidComponents, ws=self.LiquidFractions,
                       T=self.Temperature, P=self.Pressure).rhol
        else:
            return 'not a liquid!'
    @property
    def SolidDensity(self):
        if self.SolidContent > 0:
            rho = mix(self.SolidComponents, ws=self.SolidFractions,
                           T=self.Temperature, P=self.Pressure).rho
            if rho is None:
                rho = 0
                for i, c in enumerate(self.SolidComponents):
                    rho += self.SolidFractions[i] / che(c).rhos
                return 1/rho
            else:
                return mix(self.SolidComponents, ws=self.SolidFractions,
                           T=self.Temperature, P=self.Pressure).rho
        else:
            return 'not a solid!'
    @property
    def Density(self):
        if self.isElectrolyte is False:
            if self.Phase == 'g':
                return self.GasDensity
            elif self.Phase == 'l':
                return self.LiquidDensity
            elif self.Phase == 's':
                return self.SolidDensity
            else:
                if self.SolidContent > 0:
                    if self.LiquidContent > 0:
                        return 1 / (self.LiquidContent / self.LiquidDensity +
                                    self.SolidContent / self.SolidDensity)
                    elif self.GasContent > 0:
                        return 1/(self.GasContent / self.GasDensity +
                                  self.SolidContent / self.SolidDensity)
                else:
                    return 1 / (self.GasContent / self.GasDensity +
                                self.LiquidContent / self.LiquidDensity)
        else:
            return el.Laliberte_density(self.Temperature, self._MassFractions, self._CAS)

    @property
    def GasHeatCapacity(self):
        if self.GasContent > 0:
            return mix(self.GasComponents, ws=self.GasFractions,
                       T=self.Temperature, P=self.Pressure).Cpg
        else:
            return 0.0
    @property
    def LiquidHeatCapacity(self):
        if self.LiquidContent > 0:
            return mix(self.LiquidComponents, ws=self.LiquidFractions,
                       T=self.Temperature, P=self.Pressure).Cpl
        else:
            return 0.0
    @property
    def SolidHeatCapacity(self):
        if self.SolidContent > 0:
            if self._IronOxide is True:
                return mix(self.SolidComponents, ws=self.SolidFractions, T=self.Temperature, P=self.Pressure).\
                           HeatCapacitySolidMixture(T=self.Temperature,P=self.Pressure,
                                                    zs=self.SolidFractions,ws=self.SolidFractions)
            else:
                return mix(self.SolidComponents, ws=self.SolidFractions,
                           T=self.Temperature, P=self.Pressure).Cp
        else:
            return 0.0
    @property
    def HeatCapacity(self):
        if self.isElectrolyte is False:
            if self.GasContent > 0:
                if self.SolidContent > 0:
                    return (self.GasContent * self.GasHeatCapacity +
                            self.SolidContent * self.SolidHeatCapacity)
                else:
                    return self.GasHeatCapacity
            if self.LiquidContent > 0:
                if self.SolidContent > 0:
                    return (self.LiquidContent * self.LiquidHeatCapacity +
                            self.SolidContent * self.SolidHeatCapacity)
                else:
                    return self.LiquidHeatCapacity
            if self.SolidContent == 1:
                return self.SolidHeatCapacity
        else:
            return el.Laliberte_heat_capacity(self.Temperature, self._MassFractions, self._CAS)

    @property
    def GasViscosity(self):
        if self.GasContent > 0:
            if self.SolidContent != 1:
                mu = mix(self.GasComponents, ws=self.GasFractions,
                           T=self.Temperature, P=self.Pressure).mug
                return mu
            else:
                return 'not a fluid!'
    @property
    def LiquidViscosity(self):
        if self.LiquidContent > 0:
            if self.SolidContent != 1:
                if self.isElectrolyte is False:
                    mu = mix(self.LiquidComponents, ws=self.LiquidFractions,
                               T=self.Temperature, P=self.Pressure).mul
                    return mu
                else:
                    return el.Laliberte_viscosity(self.Temperature, self._MassFractions, self._CAS)
            return 'not a fluid!'
    @property
    def Viscosity(self):
        if self.Phase == 'g':
            return self.GasViscosity
        elif self.Phase == 'l':
            return self.LiquidViscosity
        elif self.Phase == 'two-phase':
            return [self.GasViscosity, self.LiquidViscosity]
        else:
            if self.GasViscosity is None:
                return self.LiquidViscosity
            else:
                return self.GasViscosity

    @property
    def GasConductivity(self):
        if self.SolidContent > 0:
            if self.GasContent > 0:
                return mix(self.GasComponents, ws=self.GasFractions,
                           T=self.Temperature, P=self.Pressure).kg
        else:
            return mix(self._Components, ws=self.MassFractions,
                       T=self.Temperature, P=self.Pressure).kg
    @property
    def LiquidConductivity(self):
        if self.SolidContent > 0:
            if self.LiquidContent > 0:
                return mix(self.LiquidComponents, ws=self.LiquidFractions,
                           T=self.Temperature, P=self.Pressure).kl
        else:
            return mix(self._Components, ws=self.MassFractions,
                       T=self.Temperature, P=self.Pressure).kl
    @property
    def Conductivity(self):
        if self.Phase == 'g':
            return self.GasConductivity
        elif self.Phase == 'l':
            return self.LiquidConductivity
        elif self.Phase == 'two-phase':
            return [self.GasConductivity, self.LiquidConductivity]

    @property
    def GasAlpha(self):
        if self.SolidContent > 0:
            if self.GasContent > 0:
                return mix(self.GasComponents, ws=self.GasFractions,
                           T=self.Temperature, P=self.Pressure).alphal
        else:
            return mix(self._Components, ws=self.MassFractions,
                       T=self.Temperature, P=self.Pressure).alphal
    @property
    def LiquidAlpha(self):
        if self.SolidContent > 0:
            if self.LiquidContent > 0:
                return mix(self.LiquidComponents, ws=self.LiquidFractions,
                           T=self.Temperature, P=self.Pressure).alphal
        else:
            return mix(self._Components, ws=self.MassFractions,
                       T=self.Temperature, P=self.Pressure).alphal
    @property
    def Alpha(self):
        if self.Phase == 'g':
            return self.GasAlpha
        elif self.Phase == 'l':
            return self.LiquidAlpha
        elif self.Phase == 'two-phase':
            return [self.GasAlpha, self.LiquidAlpha]

    @property
    def Zg(self):
        if self.SolidContent != 1:
            return mix(self._Components, ws=self.MassFractions,
                       T=self.Temperature, P=self.Pressure).Zg
        else:
            return 'not a fluid!'
    @property
    def Zl(self):
        if self.SolidContent != 1:
            return mix(self._Components, ws=self.MassFractions,
                       T=self.Temperature, P=self.Pressure).Zl
        else:
            return 'not a fluid!'
    @property
    def Z(self):
        if self.Phase == 'g':
            return self.Zg
        elif self.Phase == 'l':
            return self.Zl
        elif self.Phase == 'two-phase':
            return [self.Zg, self.Zl]

    @property
    def SurfaceTension(self):
        return mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).sigma

    @property
    def Enthalpy(self):
        return mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).H

    # PURE COMPONENT LISTS #

    @property
    def CAS(self):
        return mix(self._Components, self.MassFractions).CASs
    @property
    def Formulas(self):
        return mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).formulas
    @property
    def Elements(self):
        return self.ElementMolarFractions.keys()
    @property
    def MolarMasses(self):
        if self._wi is None:
            values = mix(self._Components, zs=self.MolarFractions,
                       T=self.Temperature, P=self.Pressure).MWs
            return OrderedDict(zip(self._Components, values))
        else:
            values = mix(self._Components, ws=self.MassFractions,
                       T=self.Temperature, P=self.Pressure).MWs
            return OrderedDict(zip(self._Components, values))
    @property
    def Densities(self):
        values = [che(i, T=self.Temperature, P=self.Pressure).rho for i in self._Components]
        return OrderedDict(zip(self._Components, values))
    @property
    def LogPs(self):
        values = mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).logPs
        return OrderedDict(zip(self._Components, values))
    @property
    def Tbs(self):
        values = mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).Tbs
        return OrderedDict(zip(self._Components, values))
    @property
    def Tms(self):
        values = mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).Tms
        return OrderedDict(zip(self._Components, values))
    @property
    def Tflashs(self):
        values = mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).Tflashs
        return OrderedDict(zip(self._Components, values))
    @property
    def Tignitions(self):
        values = mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).Tautoignitions
        return OrderedDict(zip(self._Components, values))
    @property
    def Ttriples(self):
        values = mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).Tts
        return OrderedDict(zip(self._Components, values))
    @property
    def Ptriples(self):
        values = mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).Pts
        return OrderedDict(zip(self._Components, values))
    @property
    def Tcs(self):
        values = mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).Tcs
        return OrderedDict(zip(self._Components, values))
    @property
    def Pcs(self):
        values = mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).Pcs
        return OrderedDict(zip(self._Components, values))
    @property
    def Vcs(self):
        values = mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).Vcs
        return OrderedDict(zip(self._Components, values))
    @property
    def Rhocs(self):
        values = mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).rhocs
        return OrderedDict(zip(self._Components, values))
    @property
    def Zcs(self):
        values = mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).Zcs
        return OrderedDict(zip(self._Components, values))
    @property
    def Hvaps(self):
        values = mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).Hvaps
        return OrderedDict(zip(self._Components, values))
    @property
    def Hfuss(self):
        values = mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).Hfuss
        return OrderedDict(zip(self._Components, values))
    @property
    def Hsubs(self):
        values = mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).Hsubs
        return OrderedDict(zip(self._Components, values))
    @property
    def Hfs(self):
        values = mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).Hfs
        return OrderedDict(zip(self._Components, values))
    @property
    def Hcs(self):
        values = mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).Hcs
        return OrderedDict(zip(self._Components, values))

    @property
    def gammaUNIFAC(self):
        UNIFAC_groups = mix(self.Components, ws=self.MassFractions).UNIFAC_groups
        return UNIFAC(self.Temperature, xs=self.MolarFractions, chemgroups=UNIFAC_groups)
    @property
    def gammaModUNIFAC(self):
        ModUNIFAC_groups = mix(self.Components, ws=self.MassFractions).UNIFAC_groups
        return UNIFAC(self.Temperature, xs=self.MolarFractions, chemgroups=ModUNIFAC_groups)
    @property
    def Ki(self):
        Psat_list = list(self.Psats)
        Ki_list = []
        for i, Psat in enumerate(Psat_list):
            Ki = self.gammaUNIFAC[i] * Psat / self.Pressure
            Ki_list.append(Ki)
        return Ki_list
    @property
    def Psats(self):
        Psats_list = mix(self._Components, ws=self.MassFractions,
                   T=self.Temperature, P=self.Pressure).Psats
        return Psats_list
    @property
    def Tsats(self):
        Psat_list = self.Psats
        Tsats_list = []
        for i, Psat in enumerate(Psat_list):
            def f(T, P=self.Pressure):
                self.Temperature = T
                f = P - self.Psats[i]
                return f
            Tsats_list.append(Newton(f, self.Temperature, 0.01))
        return Tsats_list

    # HIDDEN PROPERTIES #

    @property
    def _w(self):
        for i, c in enumerate(self._Components):
            if c == 'water':
                return i
    @property
    def _CAS(self):
        cas = list(self.CAS)
        if len(cas) == len(self.CAS):
            cas.remove(self.CAS[self._w])
        return cas
    @property
    def _MassFractions(self):
        wi = list(self.MassFractions)
        if len(wi) == len(self.MassFractions):
            wi.remove(self.MassFractions[self._w])
        return wi

    @property
    def _IronOxide(self):
        for i, s in enumerate(self.Components):
            if s == 'fe2o3' or s == 'Fe2O3' or s == self._Fe2O3:
                return True
    @property
    def _Components(self):
        solids = list(self.Components)
        for i, s in enumerate(solids):
            if s == 'fe2o3' or s == 'Fe2O3':
                solids[i] = self._Fe2O3
        return solids

    # FUNCTIONS #

    def Velocity(self, diameter):
        return self.VolumeFlow / Circle(diameter).Area

