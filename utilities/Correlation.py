"""Module that contains all transport phenomena correlations (dimensionless numbers)"""

from .Geometry import Circle

def Re_D(stream, diameter):
    v = stream.VolumeFlow / Circle(diameter).Area
    return stream.Density * v * diameter / stream.Viscosity

def Pr(stream):
    return stream.Viscosity * stream.HeatCapacity / stream.Conductivity

def Sc(stream, diffusivity):
    return stream.Viscosity / (diffusivity * stream.Density)

def Bi(stream, h, length):
    return h * length / stream.Conductivity

def Bi_m(diffusivity, kc, length):
    return kc * length / diffusivity

def Fo(stream, length, time):
    return stream.Alpha * time / length**2

def Fo_m(diffusivity, length, time):
    return diffusivity * time / length**2

def Nu_D(Re, Pr, constant_Ts=False):
    if Re > 2300:
        laminar = False
    else:
        laminar = True
    if laminar is False:
        m = 4.0/5.0
        n = 1.0/3.0
        C = 0.027
        return C * Re**m * Pr**n
    else:
        if constant_Ts is False:
            return 4.36
        else:
            return 3.66

def Sh_D(Re, Sc, constant_Cs=False):
    if Re > 2300:
        laminar = False
    else:
        laminar = True
    if laminar is False:
        m = 4.0 / 5.0
        n = 0.4
        C = 0.023
        return C * Re ** m * Sc ** n
    else:
        if constant_Cs is False:
            return 4.36
        else:
            return 3.66

def Nu_shell(Re, Pr):
    m = 0.55
    n = 1.0 / 3.0
    C = 0.36
    return C * Re**m * Pr**n

def Nu_plate(Re, Pr):
    m = 0.65
    n = 0.4
    C = 0.26
    return C * Re**m * Pr**n