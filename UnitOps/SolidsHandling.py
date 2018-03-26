import time
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from collections import Counter, defaultdict, OrderedDict
from Chemicals import *
from .Piping import *


class Mill:
    def __init__(self, power=None, outlet_particle_size=None, dry_grinding=False):
        self.Name = None
        self._d2 = outlet_particle_size
        self._power = power
        self.DryGrinding = dry_grinding

    def Law(self, n):
        if n == 1.0:
            return 'Kicks Law'
        elif n == 2.0:
            return 'Rittingers Law'
        elif n == 1.5:
            return 'Bonds Law'
        else:
            return 'Unknown Law!'

    def Power(self, inlet, n, K):
        if self._power is None:
            if self.DryGrinding is False:
                coef = 1.0
            else:
                coef = 4.0 / 3.0
            if n == 1.0:
                K = K * 3600
                Pot = inlet.MassFlow * K * math.log10(inlet.ParticleSize / self._d2)
            elif n == 2.0:
                K = K * 3.6
                Pot = inlet.MassFlow * K * (1 / self._d2 - 1 / inlet.ParticleSize)
            else:
                K = K * 36
                Pot = inlet.MassFlow * K * (1 / math.sqrt(self._d2) - 1 / math.sqrt(inlet.ParticleSize))
            return Pot * coef
        else:
            return self._power

    def OutletParticleSize(self, inlet, n, K):
        if self._d2 is None:
            if self.DryGrinding is False:
                coef = 1.0
            else:
                coef = 4.0 / 3.0
            E = (self._power / coef) / inlet.MassFlow
            if n == 1.0:
                K = K * 3600
                return inlet.ParticleSize * 10 ** (-E / K)
            elif n == 2.0:
                K = K * 3.6
                return (1 / inlet.ParticleSize + E / K) ** (-1.0)
            else:
                K = K * 36
                return (1 / math.sqrt(inlet.ParticleSize) + E / K) ** (-2.0)
        else:
            return self._d2

    def Outlet(self, inlet, n, K):
        O = Stream(inlet.Components)
        O.wi = inlet.MassFractions
        O.Temperature = inlet.Temperature
        O.Pressure = inlet.Pressure
        O.ParticleSize = self.OutletParticleSize(inlet, n, K)
        return O


class Electrostatic:
    def __init__(self, length, height, electric_field, number_of_plates=2):
        self.Name = None
        self.L = length
        self.H = height
        self.ElectricField = electric_field
        self.NumberOfPlates = number_of_plates
        self.Lambda = 0.066e-6

    @property
    def PlateArea(self):
        return self.L * self.H

    @property
    def NumberOfPlateSides(self):
        return 2 * (self.NumberOfPlates - 1)

    @property
    def CollectingArea(self):
        return self.PlateArea * self.NumberOfPlateSides

    def CunninghamFactor(self, inlet):
        lam = self.Lambda
        dp = inlet.ParticleSize
        Cu = 1 + (lam / dp) * (2.51 + 0.80 * math.exp(-0.55 * dp / lam))
        return Cu

    def DriftSpeed(self, inlet):
        E = self.ElectricField
        q = inlet.ParticleCharge
        Cu = self.CunninghamFactor(inlet)
        mu = inlet.GasViscosity
        dp = inlet.ParticleSize
        uh = (q * E * Cu) / (3 * math.pi * mu * dp)
        return uh

    def Efficiency(self, inlet):
        A = self.CollectingArea
        uh = self.DriftSpeed(inlet)
        Vg = inlet.VolumeFlow
        eff = 1 - math.exp(-A * uh / Vg)
        return eff

    def GasOutlet(self, inlet):
        gasflow = inlet.MassFlow * (1 - inlet.SolidContent)
        gasflows = [gasflow * wi for wi in inlet.GasFractions]
        solidflow = (1 - self.Efficiency(inlet)) * inlet.MassFlow * inlet.SolidContent
        solidflows = [solidflow * wi for wi in inlet.SolidFractions]
        totalflow = gasflow + solidflow
        gf = [gf / totalflow for gf in gasflows]
        sf = [sf / totalflow for sf in solidflows]
        GO = Stream(inlet.Components)
        GO.wi = gf + sf
        GO.Mf = totalflow
        GO.Temperature = inlet.Temperature
        GO.Pressure = inlet.Pressure
        return GO

    def SolidOutlet(self, inlet):
        solidflow = self.Efficiency(inlet) * inlet.MassFlow * inlet.SolidContent
        solidflows = [solidflow * wi for wi in inlet.SolidFractions]
        sf = [sf / solidflow for sf in solidflows]
        SO = Stream(inlet.SolidComponents)
        SO.wi = sf
        SO.Mf = solidflow
        SO.Temperature = inlet.Temperature
        SO.Pressure = inlet.Pressure
        SO.ParticleSize = inlet.ParticleSize
        return SO


class Cyclone:
    def __init__(self):
        self.Name = None
        self.stdVolumeFlow = 0.061944
        self.stdVelocity = 15.0
        self.stdViscosity = 0.000018
        self.stdDensityVariation = 2000.0
        self.stdDiameter = 0.203

    def Diameter(self, inlet):
        return math.sqrt(inlet.VolumeFlow / (self.stdVelocity * 0.5 * 0.2))

    def DensityVariation(self, inlet):
        return inlet.SolidDensity - inlet.GasDensity

    def EffectiveDiameter(self, inlet):
        return inlet.ParticleSize * (((self.Diameter(inlet) / self.stdDiameter) ** 3)
                                     * (self.stdVolumeFlow / inlet.VolumeFlow)
                                     * (self.stdDensityVariation / self.DensityVariation(inlet))
                                     * (inlet.GasViscosity / self.stdViscosity)) ** 0.5

    def Efficiency(self, inlet):
        A = 0.010065
        B = 0.053548
        n = 1.42432
        return 0.01 * (1 / (A + (B / self.EffectiveDiameter(inlet) / 1000) ** n))

    def UnderFlow(self, inlet):
        return self.Efficiency(inlet) * inlet.SolidContent * inlet.MassFlow

    def OverFlow(self, inlet):
        return inlet.MassFlow - self.UnderFlow(inlet)

    def SolidOutlet(self, inlet):
        U = Stream(inlet.SolidComponents)
        U.wi = inlet.SolidFractions
        U.Mf = self.UnderFlow(inlet)
        U.Temperature = inlet.Temperature
        U.Pressure = inlet.Pressure
        U.ParticleSize = inlet.ParticleSize
        return U
        
    def GasOutlet(self, inlet):
        gf = [wfi * inlet.GasContent * (inlet.MassFlow / self.OverFlow(inlet))
              for wfi in inlet.GasFractions]
        sf = [wsi * inlet.SolidContent * (inlet.MassFlow / self.OverFlow(inlet)) * (1 - self.Efficiency(inlet))
              for wsi in inlet.SolidFractions]
        O = Stream(inlet.GasComponents + inlet.SolidComponents)
        O.wi = gf + sf
        O.Mf = self.OverFlow(inlet)
        O.Temperature = inlet.Temperature
        O.Pressure = inlet.Pressure
        O.ParticleSize = inlet.ParticleSize
        return O


class Hydrocyclone(Cyclone):
    def __init__(self, Dc=None, d50=None):
        self.Name = None
        self._Dc = Dc
        self._d50 = d50

    def DensityVariation(self, inlet):
        return inlet.SolidDensity - inlet.LiquidDensity

    def Diameter(self, inlet):
        if self._Dc == None:
            return ((self._d50 / (2.438 * inlet.LiquidViscosity)) * (
            (inlet.VolumeFlow ** 1.2) * (self.DensityVariation(inlet)))) ** (1.0 / 3.0)
        else:
            return self._Dc

    def EffectiveDiameter(self, inlet):
        return 2.438 * (self.Diameter(inlet) ** 3) * inlet.LiquidViscosity / (
        (inlet.VolumeFlow ** 1.2) * (self.DensityVariation(inlet)))

    def Efficiency(self, inlet):
        return 1 - math.exp(-(inlet.ParticleSize / self.EffectiveDiameter(inlet) - 0.115) ** 3)

    def UnderFlow(self, inlet):
        return self.Efficiency(inlet) * inlet.SolidContent * inlet.MassFlow

    def OverFlow(self, inlet):
        return inlet.MassFlow - self.UnderFlow(inlet)

    def SolidOutlet(self, inlet):
        U = Stream(inlet.SolidComponents)
        U.wi = inlet.SolidFractions
        U.Mf = self.UnderFlow(inlet)
        U.Temperature = inlet.Temperature
        U.Pressure = inlet.Pressure
        U.ParticleSize = inlet.ParticleSize
        return U

    def LiquidOutlet(self, inlet):
        lf = [wfi * inlet.LiquidContent * (inlet.MassFlow / self.OverFlow(inlet))
              for wfi in inlet.LiquidFractions]
        sf = [wsi * inlet.SolidContent * (inlet.MassFlow / self.OverFlow(inlet)) * (1 - self.Efficiency(inlet))
              for wsi in inlet.SolidFractions]
        O = Stream(inlet.LiquidComponents + inlet.SolidComponents)
        O.wi = lf + sf
        O.Mf = self.OverFlow(inlet)
        O.Temperature = inlet.Temperature
        O.Pressure = inlet.Pressure
        O.ParticleSize = inlet.ParticleSize
        return O


class Scrubber:
    def __init__(self, number_of_nozzles, nozzle_diameter, throat_diameter):
        self.Name = None
        self.NumberOfNozzles = number_of_nozzles
        self.NozzleDiameter = nozzle_diameter
        self.ThroatDiameter = throat_diameter
        self._SafetyFactor = 1.0  # between 1-2
        self._EquilibriumSlope = 0.1

    @property
    def n(self):
        return 1 - math.exp(-self._SafetyFactor)

    @property
    def ThroatArea(self):
        return Circle(self.ThroatDiameter).Area

    @property
    def NozzleArea(self):
        return self.NumberOfNozzles * Circle(self.NozzleDiameter).Area

    # def _y1(self, gas_inlet):
    #     return gas_inlet.SolidContent
    # def _y2(self, gas_inlet):
    #     return self.AbsorptionEfficiency * self._y1(gas_inlet)
    # def _GasMolarFlux(self, gas_inlet):
    #     return gas_inlet.MolarFlow / self.ThroatArea
    # def _LiquidMolarFlux(self, gas_inlet):
    #     return self._GasMolarFlux(gas_inlet) * self._EquilibriumSlope / (self.n * self._y1(gas_inlet) / (self._y1(gas_inlet) - self._y2(gas_inlet)) - 1)
    # def _LiquidVolumeFlux(self, gas_inlet):
    #     return self._LiquidMolarFlux(gas_inlet) * 18 / 1000 / self.NozzleArea
    # def _GasLiquidRatio(self, gas_inlet):
    #     return gas_inlet.NormalVolumeFlow / (self._LiquidVolumeFlux(gas_inlet) * self.NozzleArea)

    def GasLiquidRatio(self, gas_inlet, liquid_inlet):
        return gas_inlet.NormalVolumeFlow / liquid_inlet.VolumeFlow

    def EffectiveDropletDiameter(self, gas_inlet, liquid_inlet):
        dd = 392.84 * (liquid_inlet.SurfaceTension / (
        liquid_inlet.Velocity(self.NozzleDiameter) * math.sqrt(liquid_inlet.LiquidDensity))) + \
             21.632 * ((liquid_inlet.LiquidViscosity / math.sqrt(
                 liquid_inlet.SurfaceTension * liquid_inlet.LiquidDensity)) ** 0.45) * \
             (1000 * liquid_inlet.VolumeFlow / gas_inlet.VolumeFlow)
        return dd * 1e-6

    def CollectionEfficiency(self, gas_inlet, liquid_inlet):
        A1 = 1.257
        A2 = 0.400
        A3 = 0.550
        Lambda = 1e-9
        C = 1 + (2 * Lambda / gas_inlet.ParticleSize) * (A1 + A2 * math.exp(-A3 * gas_inlet.ParticleSize / Lambda))
        Kp = C * gas_inlet.SolidDensity * gas_inlet.ParticleSize ** 2 * \
             (gas_inlet.Velocity(self.ThroatDiameter) - liquid_inlet.Velocity(self.NozzleDiameter)) \
             / (9 * gas_inlet.GasViscosity * self.EffectiveDropletDiameter(gas_inlet, liquid_inlet))
        return (Kp / (Kp + 0.7)) ** 2

    def OutletTemperature(self, gas_inlet, liquid_inlet):
        return Mixer().OutletTemperature([gas_inlet, liquid_inlet])

    def GasOutlet(self, gas_inlet, liquid_inlet):
        gasflow = gas_inlet.MassFlow * (1 - gas_inlet.SolidContent)
        gasflows = [gasflow * wi for wi in gas_inlet.GasFractions]
        solidflow = (1 - self.CollectionEfficiency(gas_inlet,
                                                   liquid_inlet)) * gas_inlet.MassFlow * gas_inlet.SolidContent
        solidflows = [solidflow * wi for wi in gas_inlet.SolidFractions]
        totalflow = gasflow + solidflow
        gf = [gf / totalflow for gf in gasflows]
        sf = [sf / totalflow for sf in solidflows]
        GO = Stream(gas_inlet.Components)
        GO.wi = gf + sf
        GO.Mf = totalflow
        GO.Temperature = self.OutletTemperature(gas_inlet, liquid_inlet)
        GO.Pressure = gas_inlet.Pressure
        GO.ParticleSize = gas_inlet.ParticleSize
        return GO

    def LiquidOutlet(self, gas_inlet, liquid_inlet):
        liquidflow = liquid_inlet.MassFlow
        liquidflows = [liquidflow * wi for wi in liquid_inlet.LiquidFractions]
        solidflow = self.CollectionEfficiency(gas_inlet, liquid_inlet) * gas_inlet.SolidContent * gas_inlet.MassFlow
        solidflows = [solidflow * wi for wi in gas_inlet.SolidFractions]
        totalflow = liquidflow + solidflow
        lf = [lf / totalflow for lf in liquidflows]
        sf = [sf / totalflow for sf in solidflows]
        LO = Stream(liquid_inlet.Components + gas_inlet.SolidComponents)
        LO.wi = lf + sf
        LO.Mf = totalflow
        LO.Temperature = self.OutletTemperature(gas_inlet, liquid_inlet)
        LO.Pressure = liquid_inlet.Pressure
        LO.ParticleSize = gas_inlet.ParticleSize
        return LO


class Decanter:
    def __init__(self, height, length_width_ratio=4):
        self.Name = None
        self.Height = height
        self.LengthWidthRatio = length_width_ratio

    def SettlingVelocity(self, inlet, porosity):
        xs = inlet.SolidContent
        rhos = inlet.SolidDensity
        rho = inlet.Density
        mu = inlet.Viscosity
        D = inlet.ParticleSize
        ut = ((rhos - rho) * porosity * 9.81 * D ** 2) * math.exp(-4.19 * (1 - porosity)) / (18 * 1000 * mu)
        return ut

    def Area(self, inlet, porosity):
        A = inlet.VolumeFlow / self.SettlingVelocity(inlet, porosity)
        return A

    def Width(self, inlet, porosity):
        W = math.sqrt(self.Area(inlet, porosity) / self.LengthWidthRatio)
        return W

    def Length(self, inlet, porosity):
        L = self.LengthWidthRatio * self.Width(inlet, porosity)
        return L

    def Volume(self, inlet, porosity):
        V = self.Length(inlet, porosity) * self.Width(inlet, porosity) * self.Height
        return V

    def SettlingTime(self, inlet, porosity):
        t = self.Volume(inlet, porosity) / inlet.VolumeFlow
        return t

    def HorizontalVelocity(self, inlet, porosity):
        uh = inlet.VolumeFlow / (self.Width(inlet, porosity) * self.Height)
        return uh

    def HydraulicRadius(self, inlet, porosity):
        W = self.Width(inlet, porosity)
        H = self.Height
        Rh = W * H / (2 * H + W)
        return Rh

    def LiquidOutlet(self, inlet):
        O = Stream(inlet.LiquidComponents)
        O.wi = inlet.LiquidFractions
        O.Mf = inlet.LiquidContent * inlet.MassFlow
        O.Temperature = inlet.Temperature
        O.Pressure = inlet.Pressure
        return O

    def SolidOutlet(self, inlet):
        U = Stream(inlet.SolidComponents)
        U.wi = inlet.SolidFractions
        U.Mf = inlet.SolidContent * inlet.MassFlow
        U.Temperature = inlet.Temperature
        U.Pressure = inlet.Pressure
        U.ParticleSize = inlet.ParticleSize
        return U

    def FlowType(self, inlet, porosity):
        Reynolds = inlet.LiquidDensity * self.HorizontalVelocity(inlet, porosity) * self.HydraulicRadius(inlet,
                                                                                                         porosity) / inlet.Viscosity
        if Reynolds <= 20000.0:
            print('Laminar! Re = ' + str(Reynolds))
            return Reynolds
        else:
            print('Turbulent! Re = ' + str(Reynolds))
            return Reynolds
