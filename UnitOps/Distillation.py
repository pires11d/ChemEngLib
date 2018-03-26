import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
from collections import Counter, defaultdict, OrderedDict
from Chemicals import *
from .Basics import Mixer


class Flash:
    def __init__(self, pressure=101325.0):
        self.Name = None
        self.Pressure = pressure
        self.Inlet = None
        # Drawing #
        self.dt = 0.1
        self.X = 0
        self.Y = 0
        self.Height = 1.5
        self.Width = 0.5

    @property
    def Ki(self):
        self.Inlet.Pressure = self.Pressure
        return self.Inlet.Ki

    @property
    def F0(self):
        acc = 0
        for i, zi in enumerate(self.Inlet.MolarFractions):
            acc += zi*self.Ki[i]
        return acc - 1

    @property
    def F1(self):
        acc = 0
        for i, zi in enumerate(self.Inlet.MolarFractions):
            acc += zi/self.Ki[i]
        return 1 - acc

    @property
    def VaporFraction(self):
        if self.F0 <= 0:
            return 0.0
        elif self.F1 >= 0:
            return 1.0
        else:
            f = 0.5
            error = 1.0
            while error > 1e-3:
                F = 0
                dF = 0
                for i, zi in enumerate(self.Inlet.MolarFractions):
                    F += zi*(self.Ki[i] - 1)/(1+f*(self.Ki[i]-1))
                    dF += (-zi*(self.Ki[i] - 1)**2.0)/((1+f*(self.Ki[i]-1))**2.0)
                fnew = f-(F/dF)
                error = (fnew-f)**2
                f = fnew
            return f

    @property
    def VaporFlow(self):
        return self.VaporFraction * self.Inlet.MolarFlow

    @property
    def LiquidFlow(self):
        return (1-self.VaporFraction) * self.Inlet.MolarFlow

    @property
    def yi(self):
        yi = []
        for i, zi in enumerate(self.Inlet.MolarFractions):
            yi.append(self.Ki[i]*zi / (1 + self.VaporFraction * (self.Ki[i] - 1)))
        return yi

    @property
    def xi(self):
        xi = []
        for i, zi in enumerate(self.Inlet.MolarFractions):
            xi.append(zi/(1+self.VaporFraction*(self.Ki[i]-1)))
        return xi

    @property
    def VaporOutlet(self):
        VO = Stream(self.Inlet.Components)
        VO.zi = self.yi
        VO.Nf = self.VaporFlow
        VO.Temperature = self.Inlet.Temperature
        VO.Pressure = self.Pressure
        VO._phase = 'g'
        return VO

    @property
    def LiquidOutlet(self):
        LO = Stream(self.Inlet.Components)
        LO.zi = self.xi
        LO.Nf = self.LiquidFlow
        LO.Temperature = self.Inlet.Temperature
        LO.Pressure = self.Pressure
        LO._phase = 'l'
        return LO

    @property
    def DrawContour(self):
        X = self.X
        Y = self.Y
        Hc = 0
        H = self.Height
        W = self.Width*1.01
        points = [[X,Y+H],[X,Y+Hc],[X+W/2,Y],[X+W,Y+Hc],[X+W,Y+H]]
        patch = plt.Polygon(points, closed=True, fill=None, lw=2, edgecolor='black')
        return patch

    @property
    def DrawLiquid(self):
        X = self.X
        Y = self.Y
        Hc = 0
        H = self.Height*(1-self.VaporFraction)
        W = self.Width*0.99
        points = [[X,Y+H],[X,Y+Hc],[X+W/2,Y],[X+W,Y+Hc],[X+W,Y+H]]
        patch = plt.Polygon(points, closed=True, fill=True, color=self.Inlet.Color)
        return patch   

    @property
    def DrawVapor(self):
        X = self.X
        Y = self.Y + self.Height*(1-self.VaporFraction)
        Hc = 0
        H = self.Height*(self.VaporFraction)
        W = self.Width*0.99
        points = [[X,Y+H],[X,Y+Hc],[X+W/2,Y],[X+W,Y+Hc],[X+W,Y+H]]
        patch = plt.Polygon(points, closed=True, fill=None, color=self.VaporOutlet.Color, hatch=self.VaporOutlet.Hatch)
        return patch 

    @property
    def DrawInletArrow(self):
        x = self.X
        y = self.Y + self.Height/2
        l = 0.1
        s = 0.05
        patch = plt.arrow(x-l,y,l,0,lw=0.5,head_width=s,head_length=s,fill=True,color='black',length_includes_head=True)
        if self.Inlet.MolarFlow > 0.0:
            patch.set_visible(True)
        else:
            patch.set_visible(False)
        return patch

    @property
    def DrawTopArrow(self):
        x = self.X + self.Width/2
        y = self.Y + self.Height
        l = 0.1
        s = 0.05
        patch = plt.arrow(x,y,0,l,lw=0.5,head_width=s,head_length=s,fill=True,color='black',length_includes_head=True)
        if self.VaporFlow > 0.0:
            patch.set_visible(True)
        else:
            patch.set_visible(False)
        return patch

    @property
    def DrawBottomArrow(self):
        x = self.X + self.Width/2
        y = self.Y
        l = 0.1
        s = 0.05
        patch = plt.arrow(x,y,0,-l,lw=0.5,head_width=s,head_length=s,fill=True,color='black',length_includes_head=True)
        if self.LiquidFlow > 0.0:
            patch.set_visible(True)
        else:
            patch.set_visible(False)
        return patch

    @property
    def FromTop(self):
        f = Flash()
        f.X = self.X
        f.Y = self.Y
        f.Width = self.Width
        f.Height = self.Height
        return f

    @property
    def FromBottom(self):
        f = Flash()
        f.X = self.X
        f.Y = self.Y + self.Height
        f.Width = self.Width
        f.Height = self.Height
        return f     


class Distiller:
    def __init__(self, max_volume, boilup_rate, pressure=101325.0):
        self.Name = None
        self.MaxVolume = max_volume
        self.Boilup = boilup_rate
        self.Pressure = pressure
        self.Mixture = None
        # Drawing #
        self.dt = 0.1
        self.X = 0 
        self.Y = 0
        self.Height = 1.5
        self.Width = 0.5

    @property
    def Ki(self):
        self.Mixture.Pressure = self.Pressure
        return self.Mixture.Ki

    @property
    def Volatility(self):
        Ki = max(self.Ki)
        Kj = min(self.Ki)
        return Ki/Kj

    @property
    def _i(self):
        for i,K in enumerate(self.Ki):
            if K == max(self.Ki):
                flag = i
                break
        return flag

    @property
    def NextMoles(self):
        N = self.Mixture.Moles
        D = self.Boilup
        N = N - D * self.dt
        if N < 0:
            N = 0
        self.Mixture.N = N
        return N

    @property
    def NextMixture(self):
        if self.Mixture.Moles == 0:
            return self.Mixture
        m = self.Mixture
        N = self.Mixture.Moles
        m.N = N
        D = self.Boilup
        x = m.MolarFractions[1-self._i]
        K = self.Ki[1-self._i]
        if N == 0:
            xi = x
        else:
            xi = (N*x - D*x*(1-K)*self.dt) / N
        if xi < 0:
            xi = 0
        if self._i == 0:
            m.zi = [xi, 1-xi]
        elif self._i == 1:
            m.zi = [1-xi, xi]
        self.Mixture = m
        return m

    @property
    def NextTime(self):
        self.NextMoles
        self.NextMixture
        return None

    @property
    def xi(self):
        return self.Mixture.MolarFractions

    @property
    def yi(self):
        alpha = self.Volatility
        xi = self.xi[self._i]
        yi = alpha * xi / (1 + xi * (alpha-1))
        if self._i == 0:
            yi_list = [yi,1-yi]
        elif self._i == 1:
            yi_list = [1-yi,yi]
        return yi_list

    @property
    def Outlet(self):
        O = Stream(self.Mixture.Components)
        O.zi = self.yi
        O.Nf = self.Boilup
        if self.Mixture.Moles == 0.0:
            O.Nf = 0.0
        O.Temperature = self.Mixture.Temperature
        O.Pressure = self.Pressure
        O._phase = 'l'
        return O

    @property
    def ConeHeight(self):
        return 0

    @property
    def LiquidHeight(self):
        return self.Height * self.Mixture.Volume / self.MaxVolume

    @property
    def LiquidWidth(self):
        return self.Width

    @property
    def DrawContour(self):
        X = self.X
        Y = self.Y
        Hc = 0
        H = self.Height
        W = self.Width*1.01
        points = [[X,Y+H],[X,Y+Hc],[X+W/2,Y],[X+W,Y+Hc],[X+W,Y+H]]
        patch = plt.Polygon(points, closed=None, fill=None, lw=2, edgecolor='black')
        return patch

    @property
    def DrawLiquid(self):
        zero = 0.01
        X = self.X
        Y = self.Y
        Hc = self.ConeHeight
        W = self.Width
        h = self.LiquidHeight
        w = (W-self.LiquidWidth)/2
        if self.LiquidWidth < self.Width:
            points = [[X+w+zero,Y+h+zero],[X+W/2,Y+zero],[X+W-w-zero,Y+h+zero]]
        else:
            points = [[X+zero*2,Y+h+zero],[X+zero*2,Y+Hc+zero],[X+W/2,Y+zero],[X+W-zero,Y+Hc+zero],[X+W-zero,Y+h+zero]]
        patch = plt.Polygon(points, closed=True, fill=True,
                            facecolor=self.Mixture.Color,
                            hatch=self.Mixture.Hatch)
        return patch

    @property
    def DrawTopArrow(self):
        x = self.X + self.Width/2
        y = self.Y + self.Height
        l = 0.1
        s = 0.05
        patch = plt.arrow(x,y,0,+l,lw=0.5,head_width=s,head_length=s,fill=True,color='black',length_includes_head=True)
        if self.Outlet.VolumeFlow > 0.0:
            patch.set_visible(True)
        else:
            patch.set_visible(False)
        return patch


class batchDistiller:
    def __init__(self, initial_moles, boiling_rate, pressure=101325.0):
        self.Name = None
        self.InitialMoles = initial_moles
        self.BoilingRate = boiling_rate
        self.Pressure = pressure
        self.dx = 0.05

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