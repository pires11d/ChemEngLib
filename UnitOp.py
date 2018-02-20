"""Module that contains all of the most common Unit Operations used in Process Industries"""


import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from Material import thermoFlowStream

from Geometry import *
from Correlation import *


#region BATCH OPERATIONS

class cylTank:
    def __init__(self, tank_diameter, tank_height, cone_angle=0.0):
        self.InitialVolume = 0.0
        self.InitialMixture = None
        self.D = tank_diameter
        self.H = tank_height
        self.Theta = cone_angle

    @property
    def ConeHeight(self):
        return Cone(diameter=self.D, angle=self.Theta).Height
    @property
    def TankHeight(self):
        return self.ConeHeight + self.H
    @property
    def ConeVolume(self):
        return Cone(diameter=self.D, angle=self.Theta).Volume
    @property
    def TankVolume(self):
        Vcyl = Cylinder(self.D, self.H).Volume
        return self.ConeVolume + Vcyl

    def LiquidVolume(self, inlets, outlet_ratio, time):
        V = self.InitialVolume
        for i in inlets:
            V += (i.VolumeFlow - outlet_ratio*i.VolumeFlow) * time
        if V < 0:
            V = 0
        if V > self.TankVolume:
            V = self.TankVolume
        return V

    def LiquidHeight(self, inlets, outlet_ratio, time):
        if self.LiquidVolume(inlets, outlet_ratio, time) < self.ConeVolume:
            h = ((3.0 * self.LiquidVolume(inlets, outlet_ratio, time) * (math.tan(math.radians(self.Theta))) ** 2.0) / math.pi) ** (1.0 / 3.0)
        else:
            h = self.ConeHeight + (self.LiquidVolume(inlets, outlet_ratio, time) - self.ConeVolume) / (0.25 * math.pi * self.D ** 2.0)
        return h

    def LiquidDiameter(self, inlets, outlet_ratio, time):
        if self.LiquidVolume(inlets, outlet_ratio, time) < self.ConeVolume:
            D = 2.0 * self.LiquidHeight(inlets, outlet_ratio, time) / math.tan(math.radians(self.Theta))
        else:
            D = self.D
        return D

    def LiquidVolumes(self, inlets, outlet_ratio_list, t_list, plot=False):
        v_list = []
        for i, t in enumerate(t_list):
            v_list.append(self.LiquidVolume(inlets, outlet_ratio_list[i], t))
        if plot is True:
            plt.plot(t_list, v_list)
            plt.show()
        else:
            return v_list

    def LiquidHeights(self, inlets, outlet_ratio_list, t_list, plot=False):
        h_list = []
        for i, t in enumerate(t_list):
            h_list.append(self.LiquidHeight(inlets, outlet_ratio_list[i], t))
        if plot is True:
            plt.plot(t_list, h_list)
            plt.show()
        else:
            return h_list

    def LiquidDiameters(self, inlets, outlet_ratio_list, t_list, plot=False):
        d_list = []
        for i, t in enumerate(t_list):
            d_list.append(self.LiquidDiameter(inlets, outlet_ratio_list[i], t))
        if plot is True:
            plt.plot(t_list, d_list)
            plt.show()
        else:
            return d_list

    def Animation(self, inlets, outlet_ratio_list, t_list):

        # Some definitions #
        D = self.D
        Hc = self.ConeHeight
        Ht = self.TankHeight
        Hs = self.TankHeight * 1.3
        maximum = max(self.D, Ht)
        h_list = self.LiquidHeights(inlets, outlet_ratio_list, t_list, plot=True)
        d_list = self.LiquidDiameters(inlets, outlet_ratio_list, t_list, plot=True)
        x_list = [tuple([-d / 2, d / 2]) for d in d_list]
        y_list = [tuple([h, h]) for h in h_list]

        # Let the plots begin #
        plt.ion()
        for i, t in enumerate(t_list):
            # Draws cylTank #
            xTk = [0, -D / 2, -D / 2, +D / 2, +D / 2, 0]
            yTk = [0, Hc, Ht, Ht, Hc, 0]
            plt.plot(xTk, yTk, color='black')
            plt.xlim(-maximum * 0.75, +maximum * 0.75)
            plt.ylim(0, maximum * 1.5)
            # Draws Liquid Level#
            if self.Theta > 0.0:
                plt.plot(x_list[i], y_list[i], color='cyan')
            else:
                plt.stackplot(x_list[i], y_list[i], color='cyan')
            # Draws Shower #
            xS = [0]
            yS = [Hs]
            plt.scatter(xS, yS, s=100.0, c='gray')
            # Draws Droplets #
            xD = [0, 0, 0, 0, 0]
            xD1 = np.linspace(0, -0.05, 5)
            xD2 = np.linspace(0, +0.05, 5)
            yD = [0.1 * Ht + 0.9 * Hs, 0.3 * Ht + 0.7 * Hs, 0.5 * Ht + 0.5 * Hs, 0.7 * Ht + 0.3 * Hs,
                  0.9 * Ht + 0.1 * Hs]
            yyD = [Hs, 0.2 * Ht + 0.8 * Hs, 0.4 * Ht + 0.6 * Hs, 0.6 * Ht + 0.4 * Hs, 0.8 * Ht + 0.2 * Hs]
            a = plt.scatter(xD, yD, s=20.0, c='cyan')
            a1 = plt.scatter(xD1, yD, s=10.0, c='cyan')
            a2 = plt.scatter(xD2, yD, s=10.0, c='cyan')
            plt.draw()
            plt.pause(0.001)
            a.remove(), a1.remove(), a2.remove()
            aa = plt.scatter(xD, yyD, s=20.0, c='cyan')
            aa1 = plt.scatter(xD1, yyD, s=10.0, c='cyan')
            aa2 = plt.scatter(xD2, yyD, s=10.0, c='cyan')
            plt.draw()
            plt.pause(0.001)
            aa.remove(), aa1.remove(), aa2.remove()
            if self.Theta == 0.0:
                plt.clf()


class recTank:
    def __init__(self, height, length, width):
        self.InitialVolume = 0.0
        self.InitialMixture = None
        self.H = height
        self.L = length
        self.W = width

    @property
    def TankVolume(self):
        return Cube(self.H, self.L, self.W).Volume

    def LiquidVolume(self, inlets, outlet_ratio, time):
        V = self.InitialVolume
        for i in inlets:
            V += (i.VolumeFlow - outlet_ratio * i.VolumeFlow) * time
        if V < 0:
            V = 0
        if V > self.TankVolume:
            V = self.TankVolume
        return V

    def LiquidHeight(self, inlets, outlet_ratio, time):
        h = self.LiquidVolume(inlets, outlet_ratio, time) / (self.L * self.W)
        return h

    def LiquidVolumes(self, inlets, outlet_ratio_list, t_list, plot=False):
        v_list = []
        for i, t in enumerate(t_list):
            v_list.append(self.LiquidVolume(inlets, outlet_ratio_list[i], t))
        if plot is True:
            plt.plot(t_list, v_list)
            plt.show()
        else:
            return v_list

    def LiquidHeights(self, inlets, outlet_ratio_list, t_list, plot=False):
        h_list = []
        for i, t in enumerate(t_list):
            h_list.append(self.LiquidHeight(inlets, outlet_ratio_list[i], t))
        if plot is True:
            plt.plot(t_list, h_list)
            plt.show()
        else:
            return h_list


class batchDistiller:
    def __init__(self, initial_moles, boiling_rate, pressure=101325.0 ,dx=0.05):
        self.InitialMoles = initial_moles
        self.BoilingRate = boiling_rate
        self.Pressure = pressure
        self.dx = dx

    def InitialVolume(self, binary_mixture):
        Vo = self.InitialMoles*binary_mixture.MolarMass/binary_mixture.Density
        return Vo

    def Temperature(self, binary_mixture):
        return binary_mixture.Temperature

    def Ki(self, binary_mixture):
        binary_mixture.Pressure = self.Pressure
        return binary_mixture.Ki

    def Volatility(self, binary_mixture):
        Ki = max(self.Ki(binary_mixture))
        Kj = min(self.Ki(binary_mixture))
        return Ki/Kj

    def InitialMolarFraction(self, binary_mixture):
        for n, Psat in enumerate(binary_mixture.Psats):
            if Psat == max(binary_mixture.Psats):
                i = n
        xo = binary_mixture.MolarFractions[i]
        return xo

    def LiquidFractionList(self, binary_mixture):
        return np.arange(self.InitialMolarFraction(binary_mixture), 0.0, -self.dx)

    def VaporFractionList(self, binary_mixture):
        alpha = self.Volatility(binary_mixture)
        x_list = self.LiquidFractionList(binary_mixture)
        y_list = []
        for i, x in enumerate(self.LiquidFractionList(binary_mixture)):
            y = alpha*x_list[i] / (1 + x_list[i]*(alpha-1))
            y_list.append(y)
        # plt.plot(x_list, y_list)
        # plt.show()
        return y_list

    def LiquidProfile(self, binary_mixture):
        x_list = self.LiquidFractionList(binary_mixture)
        xo = self.InitialMolarFraction(binary_mixture)
        Wo = self.InitialMoles
        W_list = []
        alpha = self.Volatility(binary_mixture)
        for i, x in enumerate(x_list):
            if x != 1.0 and x != 0.0:
                W_list.append(Wo*math.exp((-1/(alpha-1))*(math.log(xo/x)+alpha*math.log((1-x)/(1-xo)))))
            else:
                W_list.append(0.0)
        # plt.plot(x_list, W_list)
        # plt.show()
        return W_list

    def AverageDistillateComposition(self, binary_mixture):
        wo = self.InitialMoles
        xo = self.InitialMolarFraction(binary_mixture)
        x_list = self.LiquidFractionList(binary_mixture)
        y_list = []
        for i, w in enumerate(self.LiquidProfile(binary_mixture)):
            if x_list[i] == xo:
                y_list.append(self.VaporFractionList(binary_mixture)[0])
            else:
                y_list.append((wo*xo-w*x_list[i])/(wo-w))
        # plt.plot(x_list, y_list)
        # plt.show()
        return y_list

    def TimeList(self, binary_mixture):
        t_list = []
        wo = self.InitialMoles
        D = self.BoilingRate
        for i, w in enumerate(self.LiquidProfile(binary_mixture)):
            t_list.append((wo-w)/D)
        return t_list

    def Txy(self, binary_mixture, x1=True):
        binary_mixture.Pressure = self.Pressure
        if x1 is True:
            j = 0
        else:
            j = 1

        Tsat1 = binary_mixture.Tsats[j]
        Tsat2 = binary_mixture.Tsats[1-j]

        T_list = np.linspace(Tsat1, Tsat2, 20)
        Psats_list = []
        Psat1 = []
        Psat2 = []
        x1 = []
        y1 = []
        for i, T in enumerate(T_list):
            binary_mixture.Temperature = T
            Psats_list.append(binary_mixture.Psats)
            Psat1.append(Psats_list[i][j])
            Psat2.append(Psats_list[i][1-j])
            x1.append((binary_mixture.Pressure - Psat2[i]) / (Psat1[i] - Psat2[i]))
            y1.append(x1[i] * binary_mixture.gammaUNIFAC[1-j] * Psat1[i] / binary_mixture.Pressure)
        plt.plot(x1, T_list, y1, T_list)
        plt.show()


class batchRectifier:
    def __init__(self, initial_moles, boiling_rate, pressure=101325.0, dx=0.05):
        self.InitialMoles = initial_moles
        self.BoilingRate = boiling_rate
        self.Pressure = pressure
        self.dx = dx

#endregion

#region CONTINUOUS OPERATIONS

class Splitter:
    def __init__(self, outlet_fractions):
        self.NumberOfOutlets = len(outlet_fractions)
        self.OutletFractions = outlet_fractions

    def OutletFlows(self, inlet):
        return [inlet.MassFlow * wi for wi in self.OutletFractions]

    def OutletStreams(self, inlet):
        return [thermoFlowStream(components=inlet.Components, wi=inlet.MassFractions,
                           Mf=Fi, T=inlet.Temperature, P=inlet.Pressure, electrolyte=inlet.isElectrolyte) for Fi in self.OutletFlows(inlet)]


class Mixer:
    def __init__(self):
        self._Error = 1e9
        self._Tolerance = 0.1

    def InletFlows(self, inlets):
        return [i.MassFlow for i in inlets]

    def InletVolumeFlows(self, inlets):
        return [i.VolumeFlow for i in inlets]

    def NumberOfInlets(self, inlets):
        return len(inlets)

    def OutletFlow(self, inlets):
        return sum(self.InletFlows(inlets))

    def InletFractions(self, inlets):
        return [F/self.OutletFlow(inlets) for F in self.InletFlows(inlets)]

    def MassBalance(self, inlets):
        dc_list = []
        for i, I in enumerate(inlets):
            dc_list.append(dict(zip(I.Components, I.MassFlows)))
        # scott's method for summing n dictionaries:
        dc = dict(sum((Counter(dci) for dci in dc_list), Counter()))
        return dc

    def OutletComponents(self, inlets):
        return self.MassBalance(inlets).keys()

    def OutletFractions(self, inlets):
        return [mi/self.OutletFlow(inlets) for mi in self.MassBalance(inlets).values()]

    def OutletHeatCapacity(self, inlets):
        Cpo = 0
        for i, I in enumerate(inlets):
            Cpo += self.InletFractions(inlets)[i] * I.HeatCapacity
        return Cpo

    def OutletTemperature(self, inlets):    # TODO: review energy balance
        H_in = 0
        for I in inlets:
            H_in += I.MassFlow * I.HeatCapacity * I.Temperature
        To = H_in / (self.OutletFlow(inlets) * self.OutletHeatCapacity(inlets))
        T = To
        # while self._Error > self._Tolerance:
        #     Cp = thermoFlowStream(components=self.OutletComponents(inlets), wi=self.OutletMassFractions(inlets),
        #                     Mf=self.OutletFlow(inlets), T=T,
        #                     electrolyte=self._isElectrolyte(inlets)).HeatCapacity
        #     H = thermoFlowStream(components=self.OutletComponents(inlets), wi=self.OutletMassFractions(inlets),
        #                     Mf=self.OutletFlow(inlets), T=T,
        #                     electrolyte=self._isElectrolyte(inlets)).Enthalpy*self.OutletFlow(inlets)
        #     T = E_in / (self.OutletFlow(inlets) * Cp)
        #     self._Error = (T - To)**2
        return T

    def OutletStream(self, inlets):
        return thermoFlowStream(components=self.OutletComponents(inlets), wi=self.OutletFractions(inlets),
                          Mf=self.OutletFlow(inlets), T=self.OutletTemperature(inlets),
                          electrolyte=self._isElectrolyte(inlets))

# HIDDEN METHODS

    def _isElectrolyte(self, inlets):
        for I in inlets:
            if I.isElectrolyte is True:
                return True
            else:
                return False


class shellHX:
    def __init__(self, Ntubes, Ltube, Dtube_o, Dtube_i, Npasses, Dshell, P, B, Tc2=None, Th2=None, triangular_pitch=True):
        self.Nt = Ntubes
        self.Lt = Ltube
        self.Dto = Dtube_o
        self.Dti = Dtube_i
        self.Np = Npasses
        self.Ds = Dshell

        self.P = P
        self.B = B
        self.TriangularPitch = triangular_pitch

        self.ral = 0.15

        self._Tc2 = Tc2
        self._Th2 = Th2
        self._error = 10
        self._i = 0

    @property
    def Deq(self):
        if self.TriangularPitch is True:
            coef = 3.44
        else:
            coef = 4.00
        return (coef*(self.P**2.0) - math.pi*(self.Dto**2.0))/(math.pi*self.Dto)
    @property
    def At(self):
        return (math.pi*self.Dti**2)/4
    @property
    def Ai(self):
        return self.At * self.Nt / self.Np
    @property
    def Ao(self):
        return self.Ds * (self.P - self.Dto) * self.B / self.P
    @property
    def Area(self):
        return self.Nt * self.ral * self.Lt


    def ColdFluid(self, fluid1, fluid2):
        if fluid1.Temperature > fluid2.Temperature:
            return fluid2
        else:
            return fluid1
    def HotFluid(self, fluid1, fluid2):
        if fluid1.Temperature > fluid2.Temperature:
            return fluid1
        else:
            return fluid2

    def Q(self, fluid1, fluid2):
        if self._Th2 is None:
            mc = self.ColdFluid(fluid1,fluid2).MassFlow
            Cpc = self.ColdFluid(fluid1,fluid2).HeatCapacity
            Tc1 = self.ColdFluid(fluid1,fluid2).Temperature
            Q = mc * Cpc * (self._Tc2 - Tc1)
        else:
            mh = self.HotFluid(fluid1, fluid2).MassFlow
            Cph = self.HotFluid(fluid1, fluid2).HeatCapacity
            Th1 = self.HotFluid(fluid1, fluid2).Temperature
            Q = mh * Cph * (Th1 - self._Th2)
        return Q

    def Tc2(self, fluid1, fluid2):
        if self._Tc2 is None:
            mc = self.ColdFluid(fluid1,fluid2).MassFlow
            Cpc = self.ColdFluid(fluid1, fluid2).HeatCapacity
            Tc1 = self.ColdFluid(fluid1, fluid2).Temperature
            Tc2 = Tc1 + self.Q(fluid1, fluid2)/(mc*Cpc)
        else:
            Tc2 = self._Tc2
        return Tc2
    def Th2(self, fluid1, fluid2):
        if self._Th2 is None:
            mh = self.HotFluid(fluid1,fluid2).MassFlow
            Cph = self.HotFluid(fluid1, fluid2).HeatCapacity
            Th1 = self.HotFluid(fluid1,fluid2).Temperature
            Th2 = Th1 - self.Q(fluid1,fluid2)/(mh*Cph)
        else:
            Th2 = self._Th2
        return Th2

    def dTc(self, fluid1, fluid2):
        return self.Th2(fluid1, fluid2) - self.ColdFluid(fluid1, fluid2).Temperature
    def dTh(self, fluid1, fluid2):
        return self.HotFluid(fluid1, fluid2).Temperature - self.Tc2(fluid1, fluid2)
    def dTlm(self, fluid1, fluid2):
        dTc = self.dTc(fluid1,fluid2)
        dTh = self.dTh(fluid1, fluid2)
        return (dTh - dTc) / math.log(dTh / dTc)
    def R(self, fluid1, fluid2):
        Tc1 = self.ColdFluid(fluid1, fluid2).Temperature
        Th1 = self.HotFluid(fluid1, fluid2).Temperature
        Tc2 = self.Tc2(fluid1, fluid2)
        Th2 = self.Th2(fluid1, fluid2)
        return (Th1-Th2)/(Tc2-Tc1)
    def S(self, fluid1, fluid2):
        Tc1 = self.ColdFluid(fluid1, fluid2).Temperature
        Th1 = self.HotFluid(fluid1, fluid2).Temperature
        Tc2 = self.Tc2(fluid1, fluid2)
        return (Tc2-Tc1)/(Th1-Tc1)
    def F(self, fluid1, fluid2):
        R = self.R(fluid1, fluid2)
        S = self.S(fluid1, fluid2)
        num = math.sqrt(R**2+1)*math.log((1-S)/(1-R*S))
        den = (R-1)*math.log((2-S*(R+1-math.sqrt(R**2+1)))/(2-S*(R+1+math.sqrt(R**2+1))))
        return num/den
    def dTlm_fixed(self, fluid1, fluid2):
        return self.F(fluid1, fluid2) * self.dTlm(fluid1, fluid2)
    def Ureq(self, fluid1, fluid2):
        U = self.Q(fluid1, fluid2) / (self.Area * self.dTlm_fixed(fluid1, fluid2))
        return U


    def Tcm(self, fluid1, fluid2):
        Tc1 = self.ColdFluid(fluid1, fluid2).Temperature
        Tc2 = self.Tc2(fluid1, fluid2)
        return (Tc1+Tc2)/2
    def Thm(self, fluid1, fluid2):
        Th1 = self.HotFluid(fluid1, fluid2).Temperature
        Th2 = self.Th2(fluid1, fluid2)
        return (Th1+Th2)/2

    def OuterFluid(self, fluid1, fluid2):
        if self.ColdFluid(fluid1,fluid2).MassFlow > self.HotFluid(fluid1,fluid2).MassFlow:
            # self.HotFluid(fluid1,fluid2).Temperature = self._Thm(fluid1, fluid2)
            return self.HotFluid(fluid1,fluid2)
        else:
            # self.ColdFluid(fluid1, fluid2).Temperature = self._Tcm(fluid1, fluid2)
            return self.ColdFluid(fluid1, fluid2)
    def InnerFluid(self, fluid1, fluid2):
        if self.ColdFluid(fluid1,fluid2).MassFlow > self.HotFluid(fluid1,fluid2).MassFlow:
            # self.ColdFluid(fluid1,fluid2).Temperature = self._Tcm(fluid1, fluid2)
            return self.ColdFluid(fluid1,fluid2)
        else:
            # self.HotFluid(fluid1, fluid2).Temperature = self._Thm(fluid1, fluid2)
            return self.HotFluid(fluid1, fluid2)

    def hi(self, fluid1, fluid2):
        i = self.InnerFluid(fluid1, fluid2)
        if i is self.ColdFluid(fluid1, fluid2):
            i.Temperature = self.Tcm(fluid1, fluid2)
        else:
            i.Temperature = self.Thm(fluid1, fluid2)
        ki = i.Conductivity
        Gi = i.MassFlow / self.Ai
        Rei = Gi*self.Dti/i.Viscosity
        Nui = Nu_D(Rei, Pr(i))
        return (Nui*ki/self.Dti)*(self.Dti/self.Dto)
    def ho(self, fluid1, fluid2):
        o = self.OuterFluid(fluid1, fluid2)
        if o is self.ColdFluid(fluid1, fluid2):
            o.Temperature = self.Tcm(fluid1, fluid2)
        else:
            o.Temperature = self.Thm(fluid1, fluid2)
        ko = o.Conductivity
        Go = o.MassFlow / self.Ao
        Reo = Go*self.Deq/o.Viscosity
        Nuo = Nu_shell(Reo, Pr(o))
        return Nuo*ko/self.Deq
    def U(self, fluid1, fluid2):
        hi = self.hi(fluid1, fluid2)
        ho = self.ho(fluid1, fluid2)
        U = (hi * ho) / (hi + ho)
        return U

    # def Qo(self, fluid1, fluid2):
    #     U = self.U(fluid1, fluid2)
    #     A = self.Area
    #     dTlm = self.dTlm2(fluid1, fluid2)
    #     return U*A*dTlm
    # def Tc2real(self, fluid1, fluid2):
    #     Qo = self.Qo(fluid1, fluid2)
    #     mc = self.ColdFluid(fluid1, fluid2).MassFlow
    #     Cpc = self.ColdFluid(fluid1, fluid2).HeatCapacity
    #     Tc1 = min([fluid1.Temperature, fluid2.Temperature])
    #     Tc2 = Tc1 + Qo/(mc*Cpc)
    #     return Tc2
    # def Th2real(self, fluid1, fluid2):
    #     Qo = self.Qo(fluid1, fluid2)
    #     mh = self.HotFluid(fluid1, fluid2).MassFlow
    #     Cph = self.HotFluid(fluid1, fluid2).HeatCapacity
    #     Th1 = max([fluid1.Temperature, fluid2.Temperature])
    #     Th2 = Th1 - Qo/(mh*Cph)
    #     return Th2

    def HeatedFluidStream(self, fluid1, fluid2):
        f = self.ColdFluid(fluid1, fluid2)
        return thermoFlowStream(f.Components, wi=f.MassFractions,
                          Mf=f.MassFlow, T=self.Tc2(fluid1, fluid2))

    def CooledFluidStream(self, fluid1, fluid2):
        f = self.HotFluid(fluid1, fluid2)
        return thermoFlowStream(f.Components, wi=f.MassFractions,
                          Mf=f.MassFlow, T=self.Th2(fluid1, fluid2))


class plateHX:
    def __init__(self, U, Tc2, Th2, phase_change=False):
        self.U = U
        self.Tc2 = Tc2
        self.Th2 = Th2
        self.phase_change = phase_change

        self.Spacing = 3.0e-3
        self.Thickness = 0.75e-3
        self.Length = 1.5
        self.Width = 0.5
        self.Conductivity = 20

    @property
    def PlateArea(self):
        return self.Length * self.Width
    @property
    def CrossSectionalArea(self):
        return self.Spacing * self.Width
    @property
    def De(self):
        return 2 * self.Spacing

    def ColdFluid(self, fluid1, fluid2):
        if fluid1.Temperature > fluid2.Temperature:
            return fluid2
        else:
            return fluid1
    def HotFluid(self, fluid1, fluid2):
        if fluid1.Temperature > fluid2.Temperature:
            return fluid1
        else:
            return fluid2

    def Q(self, fluid1, fluid2):
        if self.phase_change is False:
            mc = self.ColdFluid(fluid1,fluid2).MassFlow
            Cpc = self.ColdFluid(fluid1,fluid2).HeatCapacity
            Tc1 = self.ColdFluid(fluid1,fluid2).Temperature
            Q = mc * Cpc * (self.Tc2 - Tc1)
        else:
            mh = self.HotFluid(fluid1, fluid2).MassFlow
            #Cph = self.HotFluid(fluid1, fluid2).HeatCapacity
            #Th1 = self.HotFluid(fluid1, fluid2).Temperature
            #Q = mh * Cph * (Th1 - self.Th2)
            dHvs = self.HotFluid(fluid1, fluid2).Hvaps.values()
            dHv = sum(dHvs)/len(dHvs)
            Q = mh * dHv
        return Q

    def dTc(self, fluid1, fluid2):
        return self.Th2 - self.ColdFluid(fluid1, fluid2).Temperature
    def dTh(self, fluid1, fluid2):
        return self.HotFluid(fluid1, fluid2).Temperature - self.Tc2
    def dTlm(self, fluid1, fluid2):
        dTc = self.dTc(fluid1,fluid2)
        dTh = self.dTh(fluid1, fluid2)
        return (dTh - dTc) / math.log(dTh / dTc)
    def R(self, fluid1, fluid2):
        Tc1 = self.ColdFluid(fluid1, fluid2).Temperature
        Th1 = self.HotFluid(fluid1, fluid2).Temperature
        Tc2 = self.Tc2
        Th2 = self.Th2
        return (Th1-Th2)/(Tc2-Tc1)
    def S(self, fluid1, fluid2):
        Tc1 = self.ColdFluid(fluid1, fluid2).Temperature
        Th1 = self.HotFluid(fluid1, fluid2).Temperature
        Tc2 = self.Tc2
        return (Tc2-Tc1)/(Th1-Tc1)
    def F(self, fluid1, fluid2):
        R = self.R(fluid1, fluid2)
        S = self.S(fluid1, fluid2)
        num = math.sqrt(R**2+1)*math.log((1-S)/(1-R*S))
        den = (R-1)*math.log((2-S*(R+1-math.sqrt(R**2+1)))/(2-S*(R+1+math.sqrt(R**2+1))))
        return num/den
    def dTlm_fixed(self, fluid1, fluid2):
        try:
            return self.F(fluid1, fluid2) * self.dTlm(fluid1, fluid2)
        except:
            return self.dTlm(fluid1, fluid2)

    def Area(self, fluid1, fluid2):
        A = self.Q(fluid1, fluid2) / (self.U * self.dTlm_fixed(fluid1, fluid2))
        return A

    def Tcm(self, fluid1, fluid2):
        Tc1 = self.ColdFluid(fluid1, fluid2).Temperature
        Tc2 = self.Tc2(fluid1, fluid2)
        return (Tc1+Tc2)/2
    def Thm(self, fluid1, fluid2):
        Th1 = self.HotFluid(fluid1, fluid2).Temperature
        Th2 = self.Th2(fluid1, fluid2)
        return (Th1+Th2)/2

    def h1(self, fluid1, fluid2):
        f1 = fluid1
        k1 = f1.Conductivity
        G1 = f1.MassFlow / self.CrossSectionalArea
        Re1 = G1 * self.De / f1.Viscosity
        Nu1 = Nu_plate(Re1, Pr(f1))
        return Nu1 * k1 / self.De
    def h2(self, fluid1, fluid2):
        f2 = fluid2
        k2 = f2.Conductivity
        G2 = f2.MassFlow / self.CrossSectionalArea
        Re2 = G2 * self.De / f2.Viscosity
        Nu2 = Nu_plate(Re2, Pr(f2))
        return Nu2 * k2 / self.De

    def U_real(self, fluid1, fluid2):
        h1 = self.h1(fluid1, fluid2)
        h2 = self.h2(fluid1, fluid2)
        k = self.Conductivity
        L = self.Thickness
        U = 1 / h1 + L / k  + 1 / h2
        return 1/U

    def HeatedFluidStream(self, fluid1, fluid2):
        f = self.ColdFluid(fluid1, fluid2)
        return thermoFlowStream(f.Components, wi=f.MassFractions,
                          Mf=f.MassFlow, T=self.Tc2(fluid1, fluid2))

    def CooledFluidStream(self, fluid1, fluid2):
        f = self.HotFluid(fluid1, fluid2)
        return thermoFlowStream(f.Components, wi=f.MassFractions,
                          Mf=f.MassFlow, T=self.Th2(fluid1, fluid2))


class Flash:
    def __init__(self, pressure=101325.0):
        self.Pressure = pressure

    def Ki(self, inlet):
        inlet.Pressure = self.Pressure
        return inlet.Ki

    def F0(self, inlet):
        acc = 0
        for i, zi in enumerate(inlet.MolarFractions):
            acc += zi*self.Ki(inlet)[i]
        return acc - 1

    def F1(self, inlet):
        acc = 0
        for i, zi in enumerate(inlet.MolarFractions):
            acc += zi/self.Ki(inlet)[i]
        return 1 - acc

    def VaporFraction(self, inlet):
        if self.F0(inlet) <= 0:
            return 0.0
        elif self.F1(inlet) >= 0:
            return 1.0
        else:
            f = 0.5
            error = 1.0
            while error > 1e-3:
                F = 0
                dF = 0
                for i, zi in enumerate(inlet.MolarFractions):
                    F += zi*(self.Ki(inlet)[i] - 1)/(1+f*(self.Ki(inlet)[i]-1))
                    dF += (-zi*(self.Ki(inlet)[i] - 1)**2.0)/((1+f*(self.Ki(inlet)[i]-1))**2.0)
                fnew = f-(F/dF)
                error = (fnew-f)**2
                f = fnew
            return f

    def VaporFlow(self, inlet):
        return self.VaporFraction(inlet) * inlet.MolarFlow

    def LiquidFlow(self, inlet):
        return (1-self.VaporFraction(inlet)) * inlet.MolarFlow

    def yi(self, inlet):
        yi = []
        for i, zi in enumerate(inlet.MolarFractions):
            yi.append(self.Ki(inlet)[i]*zi / (1 + self.VaporFraction(inlet) * (self.Ki(inlet)[i] - 1)))
        return yi

    def xi(self, inlet):
        xi = []
        for i, zi in enumerate(inlet.MolarFractions):
            xi.append(zi/(1+self.VaporFraction(inlet)*(self.Ki(inlet)[i]-1)))
        return xi

    def VaporOutletStream(self, inlet):
        return thermoFlowStream(inlet.Components, zi=self.yi(inlet), Nf=self.VaporFlow(inlet),
                          T=inlet.Temperature, P=self.Pressure)

    def LiquidOutletStream(self, inlet):
        return thermoFlowStream(inlet.Components, zi=self.xi(inlet), Nf=self.LiquidFlow(inlet),
                          T=inlet.Temperature, P=self.Pressure)

# SOLID-PHASE #

class Mill:
    def __init__(self, power=None, outlet_particle_size=None, dry_grinding=False):
        self._d2 = outlet_particle_size
        self._power = power
        self.DryGrinding = dry_grinding

    def Law(self, n):
        if n == 1.0:
            return 'Kicks Law'
        elif n == 2.0:
            return 'Rittingers Law'
        else:
            return 'Bonds Law'

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

    def OutletStream(self, inlet, n, K):
        O = thermoFlowStream(inlet.Components, wi=inlet.MassFractions, T=inlet.Temperature)
        O.ParticleSize = self.OutletParticleSize(inlet, n, K)
        return O


class Electrostatic:
    def __init__(self, length, height, electric_field, number_of_plates=2):
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

    def GasOutletStream(self, inlet):
        gasflow = inlet.MassFlow * (1 - inlet.SolidContent)
        gasflows = [gasflow * wi for wi in inlet.GasFractions]
        solidflow = (1 - self.Efficiency(inlet)) * inlet.MassFlow * inlet.SolidContent
        solidflows = [solidflow * wi for wi in inlet.SolidFractions]
        totalflow = gasflow + solidflow
        gf = [gf / totalflow for gf in gasflows]
        sf = [sf / totalflow for sf in solidflows]
        return thermoFlowStream(inlet.Components, wi=gf + sf, Mf=totalflow,
                          T=inlet.Temperature, P=inlet.Pressure)

    def SolidOutletStream(self, inlet):
        solidflow = self.Efficiency(inlet) * inlet.MassFlow * inlet.SolidContent
        solidflows = [solidflow * wi for wi in inlet.SolidFractions]
        sf = [sf / solidflow for sf in solidflows]
        return thermoFlowStream(inlet.SolidComponents, wi=sf, Mf=solidflow,
                          T=inlet.Temperature, P=inlet.Pressure)


class Cyclone:
    def __init__(self):
        self.stdVolumeFlow = 0.061944
        self.stdVelocity = 15.0
        self.stdViscosity = 0.000018
        self.stdDensityVariation = 2000.0
        self.stdCycloneDiameter = 0.203

    def CycloneDiameter(self, inlet):
        return math.sqrt(inlet.VolumeFlow / (self.stdVelocity * 0.5 * 0.2))

    def DensityVariation(self, inlet):
        return inlet.SolidDensity - inlet.GasDensity

    def EffectiveDiameter(self, inlet):
        return inlet.ParticleSize * (((self.CycloneDiameter(inlet) / self.stdCycloneDiameter) ** 3)
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

    def UnderFlowStream(self, inlet):
        return thermoFlowStream(components=inlet.SolidComponents, wi=inlet.SolidFractions,
                          Mf=self.UnderFlow(inlet), T=inlet.Temperature, P=inlet.Pressure)

    def OverFlowStream(self, inlet):
        gf = [wfi * inlet.GasContent * (inlet.MassFlow / self.OverFlow(inlet))
              for wfi in inlet.GasFractions]
        sf = [wsi * inlet.SolidContent * (inlet.MassFlow / self.OverFlow(inlet)) * (1 - self.Efficiency(inlet))
              for wsi in inlet.SolidFractions]
        return thermoFlowStream(components=(inlet.GasComponents + inlet.SolidComponents), wi=gf + sf,
                          Mf=self.OverFlow(inlet), T=inlet.Temperature, P=inlet.Pressure)


class Hydrocyclone:
    def __init__(self, Dc=None, d50=None):
        self._Dc = Dc
        self._d50 = d50

    def DensityVariation(self, inlet):
        return inlet.SolidDensity - inlet.LiquidDensity

    def CycloneDiameter(self, inlet):
        if self._Dc == None:
            return ((self._d50 / (2.438 * inlet.LiquidViscosity)) * (
            (inlet.VolumeFlow ** 1.2) * (self.DensityVariation(inlet)))) ** (1.0 / 3.0)
        else:
            return self._Dc

    def EffectiveDiameter(self, inlet):
        return 2.438 * (self.CycloneDiameter(inlet) ** 3) * inlet.LiquidViscosity / (
        (inlet.VolumeFlow ** 1.2) * (self.DensityVariation(inlet)))

    def Efficiency(self, inlet):
        return 1 - math.exp(-(inlet.ParticleSize / self.EffectiveDiameter(inlet) - 0.115) ** 3)

    def UnderFlow(self, inlet):
        return self.Efficiency(inlet) * inlet.SolidContent * inlet.MassFlow

    def OverFlow(self, inlet):
        return inlet.MassFlow - self.UnderFlow(inlet)

    def UnderFlowStream(self, inlet):
        return thermoFlowStream(components=inlet.SolidComponents, wi=inlet.SolidFractions,
                          Mf=self.UnderFlow(inlet), T=inlet.Temperature, P=inlet.Pressure)

    def OverFlowStream(self, inlet):
        lf = [wfi * inlet.LiquidContent * (inlet.MassFlow / self.OverFlow(inlet))
              for wfi in inlet.LiquidFractions]
        sf = [wsi * inlet.SolidContent * (inlet.MassFlow / self.OverFlow(inlet)) * (1 - self.Efficiency(inlet))
              for wsi in inlet.SolidFractions]
        return thermoFlowStream(components=(inlet.LiquidComponents + inlet.SolidComponents), wi=lf + sf,
                          Mf=self.OverFlow(inlet), T=inlet.Temperature, P=inlet.Pressure)


class Venturi:
    def __init__(self, number_of_nozzles, nozzle_diameter, throat_diameter):
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

    def GasOutletStream(self, gas_inlet, liquid_inlet):
        gasflow = gas_inlet.MassFlow * (1 - gas_inlet.SolidContent)
        gasflows = [gasflow * wi for wi in gas_inlet.GasFractions]
        solidflow = (1 - self.CollectionEfficiency(gas_inlet,
                                                   liquid_inlet)) * gas_inlet.MassFlow * gas_inlet.SolidContent
        solidflows = [solidflow * wi for wi in gas_inlet.SolidFractions]
        totalflow = gasflow + solidflow
        gf = [gf / totalflow for gf in gasflows]
        sf = [sf / totalflow for sf in solidflows]
        return thermoFlowStream(gas_inlet.Components, wi=gf + sf, Mf=totalflow,
                          T=self.OutletTemperature(gas_inlet, liquid_inlet), P=gas_inlet.Pressure)

    def LiquidOutletStream(self, gas_inlet, liquid_inlet):
        liquidflow = liquid_inlet.MassFlow
        liquidflows = [liquidflow * wi for wi in liquid_inlet.LiquidFractions]
        solidflow = self.CollectionEfficiency(gas_inlet, liquid_inlet) * gas_inlet.SolidContent * gas_inlet.MassFlow
        solidflows = [solidflow * wi for wi in gas_inlet.SolidFractions]
        totalflow = liquidflow + solidflow
        lf = [lf / totalflow for lf in liquidflows]
        sf = [sf / totalflow for sf in solidflows]
        return thermoFlowStream(liquid_inlet.Components + gas_inlet.SolidComponents, wi=lf + sf, Mf=totalflow,
                          T=self.OutletTemperature(gas_inlet, liquid_inlet))


class Decanter:
    def __init__(self, height, length_width_ratio=4):
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

    def OverFlowStream(self, inlet):
        return thermoFlowStream(inlet.LiquidComponents, wi=[inlet.LiquidComponents],
                          Mf=inlet.LiquidContent * inlet.MassFlow,
                          T=inlet.Temperature, P=inlet.Pressure)

    def UnderFlowStream(self, inlet):
        return thermoFlowStream(inlet.SolidComponents, wi=[inlet.SolidComponents], Mf=inlet.SolidContent * inlet.MassFlow,
                          T=inlet.Temperature, P=inlet.Pressure)

    def FlowType(self, inlet, porosity):
        Reynolds = inlet.LiquidDensity * self.HorizontalVelocity(inlet, porosity) * self.HydraulicRadius(inlet,
                                                                                                         porosity) / inlet.Viscosity
        if Reynolds <= 20000.0:
            print('Laminar! Re = ' + str(Reynolds))
            return Reynolds
        else:
            print('Turbulent! Re = ' + str(Reynolds))
            return Reynolds

            # endregion

#endregion








