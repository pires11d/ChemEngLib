import math
import matplotlib.pyplot as plt
from matplotlib import animation
from collections import Counter, defaultdict, OrderedDict
from Chemicals import *
from .Piping import Mixer


class shellHX:
    def __init__(self, Ntubes, Ltube, Dtube_o, Dtube_i, Npasses, Dshell, P, B, Tc2=None, Th2=None, triangular_pitch=True):
        self.Name = None
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
            Cpc = self.ColdFluid(fluid1,fluid2).SpecificHeat
            Tc1 = self.ColdFluid(fluid1,fluid2).Temperature
            Q = mc * Cpc * (self._Tc2 - Tc1)
        else:
            mh = self.HotFluid(fluid1, fluid2).MassFlow
            Cph = self.HotFluid(fluid1, fluid2).SpecificHeat
            Th1 = self.HotFluid(fluid1, fluid2).Temperature
            Q = mh * Cph * (Th1 - self._Th2)
        return Q

    def Tc2(self, fluid1, fluid2):
        if self._Tc2 is None:
            mc = self.ColdFluid(fluid1,fluid2).MassFlow
            Cpc = self.ColdFluid(fluid1, fluid2).SpecificHeat
            Tc1 = self.ColdFluid(fluid1, fluid2).Temperature
            Tc2 = Tc1 + self.Q(fluid1, fluid2)/(mc*Cpc)
        else:
            Tc2 = self._Tc2
        return Tc2
    def Th2(self, fluid1, fluid2):
        if self._Th2 is None:
            mh = self.HotFluid(fluid1,fluid2).MassFlow
            Cph = self.HotFluid(fluid1, fluid2).SpecificHeat
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

    def HeatedFluidStream(self, fluid1, fluid2):
        f = self.ColdFluid(fluid1, fluid2)
        hf = Stream(f.Components)
        hf.wi = f.MassFractions
        hf.Mf = f.MassFlow
        hf.Temperature = self.Tc2(fluid1,fluid2)
        hf.Pressure = f.Pressure
        return hf

    def CooledFluidStream(self, fluid1, fluid2):
        f = self.HotFluid(fluid1, fluid2)
        cf = Stream(f.Components)
        cf.wi = f.MassFractions
        cf.Mf = f.MassFlow
        cf.Temperature = self.Th2(fluid1,fluid2)
        cf.Pressure = f.Pressure
        return cf


class plateHX:
    def __init__(self, U, Tc2, Th2, phase_change=False):
        self.Name = None
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
            Cpc = self.ColdFluid(fluid1,fluid2).SpecificHeat
            Tc1 = self.ColdFluid(fluid1,fluid2).Temperature
            Q = mc * Cpc * (self.Tc2 - Tc1)
        else:
            mh = self.HotFluid(fluid1, fluid2).MassFlow
            #Cph = self.HotFluid(fluid1, fluid2).SpecificHeat
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
        hf = Stream(f.Components)
        hf.wi = f.MassFractions
        hf.Mf = f.MassFlow
        hf.Temperature = self.Tc2(fluid1,fluid2)
        hf.Pressure = f.Pressure
        return hf

    def CooledFluidStream(self, fluid1, fluid2):
        f = self.HotFluid(fluid1, fluid2)
        cf = Stream(f.Components)
        cf.wi = f.MassFractions
        cf.Mf = f.MassFlow
        cf.Temperature = self.Th2(fluid1,fluid2)
        cf.Pressure = f.Pressure
        return cf