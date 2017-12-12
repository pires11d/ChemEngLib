class Mass:
    def __init__(self,amount):
        self.amount = amount
    @property
    def g(self):
        return self.amount * 0.001
    @property
    def lb(self):
        return self.amount * 0.45359237
    @property
    def oz(self):
        return self.amount * 0.02834952
    @property
    def ton(self):
        return self.amount * 1000
    @property
    def ton_long(self):
        return self.amount * 1016.047
    @property
    def ton_short(self):
        return self.amount * 907.18474


class Time:
    def __init__(self,amount):
        self.amount = amount
    @property
    def min(self):
        return self.amount*60
    @property
    def h(self):
        return self.amount * 3600


class Temperature:
    def __init__(self,amount):
        self.amount = amount
    @property
    def C(self):
        return self.amount + 273.15
    @property
    def F(self):
        return (self.amount-32)/1.8 + 273.15
    @property
    def R(self):
        return (self.amount-459.67-32)/1.8 + 273.15


class Length:
    def __init__(self,amount):
        self.amount = amount
    @property
    def A(self):
        return self.amount * 1e-10
    @property
    def nm(self):
        return self.amount * 1e-9
    @property
    def um(self):
        return self.amount * 1e-6
    @property
    def mm(self):
        return self.amount * 1e-3
    @property
    def cm(self):
        return self.amount * 1e-2
    @property
    def dm(self):
        return self.amount * 1e-1
    @property
    def km(self):
        return self.amount * 1e3
    @property
    def ft(self):
        return self.amount * 0.3048
    @property
    def inch(self):
        return self.amount * 0.0254


class Area:
    def __init__(self,amount):
        self.amount = amount
    @property
    def mm2(self):
        return self.amount * 1e-6
    @property
    def cm2(self):
        return self.amount * 1e-4
    @property
    def ft2(self):
        return self.amount * 0.09290304
    @property
    def acre(self):
        return self.amount * 4046.873
    @property
    def ha(self):
        return self.amount * 1e4


class Volume:
    def __init__(self,amount):
        self.amount = amount
    @property
    def mm3(self):
        return self.amount * 1e-9
    @property
    def cm3(self):
        return self.amount * 1e-6
    @property
    def mL(self):
        return self.amount * 1e-6
    @property
    def dm3(self):
        return self.amount * 1e-3
    @property
    def L(self):
        return self.amount * 1e-3
    @property
    def ft3(self):
        return self.amount * 0.02831685
    @property
    def gal(self):
        return self.amount * 3.785412e-3


class Force:
    def __init__(self,amount):
        self.amount = amount
    @property
    def kgf(self):
        return self.amount * 9.80665
    @property
    def lbf(self):
        return self.amount * 4.448222


class Pressure:
    def __init__(self,amount):
        self.amount = amount
    @property
    def kPa(self):
        return self.amount * 1e3
    @property
    def bar(self):
        return self.amount * 1e5
    @property
    def mbar(self):
        return self.amount * 1e2
    @property
    def atm(self):
        return self.amount * 101325
    @property
    def psi(self):
        return self.amount * 6.894757e-3
    @property
    def mH2O(self):
        return self.amount * 9806.38
    @property
    def mmH2O(self):
        return self.amount * 9.80638
    @property
    def inH2O(self):
        return self.amount * 0.249082e3
    @property
    def mmHg(self):
        return self.amount * 133.3224
    @property
    def inHg(self):
        return self.amount * 3.38638e3
    @property
    def Torr(self):
        return self.amount * 133.3224


class Energy:
    def __init__(self,amount):
        self.amount = amount
    @property
    def kJ(self):
        return self.amount * 1e3
    @property
    def cal(self):
        return self.amount * 4.184
    @property
    def kcal(self):
        return self.amount * 4184
    @property
    def kWh(self):
        return self.amount * 3.6e6
    @property
    def BTU(self):
        return self.amount * 1055.056
    @property
    def eV(self):
        return self.amount * 1.6e-19


class Power:
    def __init__(self,amount):
        self.amount = amount
    @property
    def kW(self):
        return self.amount * 1e3
    @property
    def HP(self):
        return self.amount * 745.70
    @property
    def BTU_h(self):
        return self.amount * 0.2930711


class Viscosity:
    def __init__(self,amount):
        self.amount = amount
    @property
    def P(self):
        return self.amount * 1e-1
    @property
    def cP(self):
        return self.amount * 1e-3
    @property
    def mNs_m2(self):
        return self.amount * 1e-3


class SurfaceTension:
    def __init__(self,amount):
        self.amount = amount
    @property
    def mN_m(self):
        return self.amount * 1e-3
    @property
    def dyn_cm(self):
        return self.amount * 1e-3
    @property
    def gf_cm(self):
        return self.amount * 0.980665
    @property
    def kgf_cm(self):
        return self.amount * 0.980665e3
    @property
    def lbf_in(self):
        return self.amount * 175.12683699