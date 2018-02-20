import time
from Converter import *
from Correlation import *
from Numerical import *
from UnitOp import *
from thermoChemical import *


# TODO: Txy with deviation
# TODO: Multicomponent 'batchRectifier'
# TODO: Adsorber


#region FLOW STREAMS

air = thermoFlowStream(['N2', 'O2'], wi=[0.7, 0.3], Vfo=0.007, T=300.0)
gas = thermoFlowStream(['N2', 'H2'], zi=[0.9, 0.1], Vfo=0.007, T=400.0)

water = thermoFlowStream(['water'], wi=[1.0], Mf=1.0, T=25+273)
water2 = thermoFlowStream(['water'], wi=[1.0], Mf=1.0, T=95+273)
water3 = thermoFlowStream(['water'], wi=[1.0], Mf=1.0, T=100+273)

brine = thermoFlowStream(['water', 'NaCl'], wi=[0.9, 0.1], Mf=2.0, electrolyte=True, T=300.0)
acid = thermoFlowStream(['water', 'hcl'], wi=[0.8, 0.1, 0.1], Mf=1.0, electrolyte=True, T=320.0)
liquor = thermoFlowStream(['water', 'hcl', 'FeCl2'], wi=[0.7, 0.2, 0.1], Mf=1.0, electrolyte=True, T=320.0)
feed = thermoFlowStream(['ethanol', 'water'], zi=[0.5, 0.5], Nf=1.0, T=354.0)
mix = thermoFlowStream(['ethanol', 'water'], zi=[0.5, 0.5], Nf=0.0, T=345.0)

oxide = thermoFlowStream(['Fe2O3', 'FeCl2'], wi=[0.9, 0.1], Mf=5.0)
iron = thermoFlowStream(['iron'], wi=[1.0], Mf=Mass(100).ton/3600.0)
iron.ParticleSize = 100e-3

G1 = thermoFlowStream(['N2', 'O2', 'aluminum', 'iron'], wi=[0.7, 0.20, 0.06, 0.04], Vf=0.033)
G1.ParticleSize = 0.5e-6
G1.ParticleCharge = Energy(10).eV

L1 = thermoFlowStream(['water', 'sodium sulphate'], wi=[0.8, 0.2], Vf=0.25)
L1.ParticleSize = 2e-3

GL = thermoFlowStream(['water', 'N2', 'O2'], wi=[0.25, 0.50, 0.25], Mf=1.0)

#endregion

#region UNIT OPERATIONS

t_list = np.linspace(0,200,200)
outlet_ratio_list = [0*t for t in t_list]

cTk = cylTank(tank_diameter=1.0, tank_height=1.0, cone_angle=30.0)
# cTk.LiquidHeights([water],outlet_ratio_list,t_list,plot=True)
# cTk.Animation([water],outlet_ratio_list,t_list)

rTk = recTank(height=1.0, length=4.0, width=1.0)
# rTk.LiquidHeights([water],outlet_ratio_list,t_list,plot=True)

cy = Hydrocyclone(Dc=0.16)
# print "inlet flow ", cy.InletFlow(F1)
# print "underflow ", cy.UnderFlow(F1)
# print "overflow ", cy.OverFlow(F1)
# print "effective diameter", cy.EffectiveDiameter(F1)
# print "efficiency ", cy.Efficiency(F1)
# print "cyclone diameter", cy.CycloneDiameter(F1)
# print cy.OverFlowFluidFractions(F2), cy.OverFlowSolidFractions(F2)
# print cy.UnderthermoFlowStream(F2).MassFractions
# print cy.OverthermoFlowStream(F2).MassFractions

sp = Splitter(outlet_fractions=[0.2, 0.8])
# print sp.OutletFlows(acid)
# print sp.NumberOfOutlets
# Fi = sp.OutletthermoFlowStreams(acid)
# print Fi[0].MassFlow

mi = Mixer()
# Fout = mi.OutletStream([gas, air])
# print Fout.Temperature
# print Fout.GasContent

ve = Venturi(number_of_nozzles=4, nozzle_diameter=0.086, throat_diameter=0.440)
# print ve.GasLiquidRatio(G1, water)
# print ve.EffectiveDropletDiameter(G1, water)
# print ve.CollectionEfficiency(G1, water)
# G2 = ve.GasOutletStream(G1, water)
# print G2.MassFlow
# L2 = ve.LiquidOutletStream(G1, water)
# print L2.MassFlow

m1 = Mill(outlet_particle_size=50e-3)
# print m1.Power(iron, 1.0, K=0.57)
# print m1.Power(iron, 2.0, K=17.0)
# print m1.Power(iron, 1.5, Wi=12.7)

m2 = Mill(power=17000.0)
# print m2.OutletParticleSize(iron, 1.0, K=0.57)
# print m2.OutletParticleSize(iron, 2.0, K=17.0)
# print m2.OutletParticleSize(iron, 1.5, Wi=12.7)
# print m2.OutletStream(iron, 1.5, Wi=12.7).ParticleSize

hx = shellHX(Tc2=313.0, Ntubes=302, Ltube=Length(16).ft,
             Dtube_o=Length(0.062).ft, Dtube_i=Length(0.049).ft,
             Npasses=2.0, Dshell=Length(1.771).ft, P=Length(1.0).inch, B=Length(5.0).inch)
# print "Th2: ", hx.Th2(water, water2)
# print "Tc2: ", hx.Tc2(water, water2)

hx2 = plateHX(1500,55+273,40+273)
#print(hx2.U_real(water,water2))
#print(hx2.Area(water,water2))

hx3 = plateHX(1500,55+273,100+273,True)
#print(hx3.U_real(water,water3))
#print(hx3.Area(water,water3))

dc = Decanter(height=4.5)
# print(dc.SettlingVelocity(L1,0.8))
# print(dc.Width(L1,0.8))
# print(dc.Length(L1, 0.8))
# print(dc.Volume(L1,0.8))
# print(dc.HorizontalVelocity(L1,0.8))
# print(dc.OverthermoFlowStream(L1).Components)
# print(dc.UnderthermoFlowStream(L1).Components)
# print(dc.HydraulicRadius(L1, 0.8))
# dc.FlowType(L1,0.8)

fl = Flash(pressure=Pressure(1.05).bar)
# print (fl.VaporFraction(feed))
# print (fl.Ki(feed))
# print (fl.yi(feed))
# print (fl.xi(feed))

bD = batchDistiller(initial_moles=100, boiling_rate=1, dx=0.025)
# print(bD.InitialVolume(mix))
# print(bD.LiquidFractionList(mix))
# print(bD.InitialMolarFraction(mix))
# print(bD.Volatility(mix))
# print(bD.LiquidProfile(mix))
# print(bD.VaporFractionList(mix))
# print(bD.AverageDistillateComposition(mix))
# print(bD.TimeList(mix))
# bD.Txy(mix)
# plt.plot(bD.TimeList(mix), bD.AverageDistillateComposition(mix))
# plt.show()

esp = Electrostatic(length=6, height=3, electric_field=50000, number_of_plates=2)
# print(esp.CollectingArea)
# print(esp.CunninghamFactor(G1))
# print(esp.DriftSpeed(G1))
# print(esp.Efficiency(G1))
# print(esp.SolidOutletStream(G1).MassFractions, esp.SolidOutletStream(G1).Components)
# print(esp.GasOutletStream(G1).MassFractions, esp.GasOutletStream(G1).Components)

#endregion

#region INTEGRATION EXAMPLES

# def func(x):
#     f = x**2-5*x+6
#     return f
#
# def dFdt(t,y):
#     dF = -0.1*y
#     return dF

# x_list = np.arange(1,100,10)
# y_list = [func(x) for x in x_list]

# y = Euler(dFdt, 5, x_list)
# yy = RungeKutta(dFdt, 5, x_list)
# plt.plot(x_list,y, x_list, yy)
# plt.show()

# plt.plot(x_list, y_list)
# plt.show()
# trap = Trapezoidal(func,x_list)
# simp = Simpson(func,x_list)
# print(trap, simp)

# x = NewtonRaphson(func,298,0.01,1e-3)
# y = QuasiNewton(func,298,0.01,1e-3)
# x = Newt(func,298,0.01)
# y = qNewt(func,298,0.01)

#endregion

#region BATCH DISTILLATION INTEGRATION

# def f(x):
#     alpha = bD.Volatility(mix)
#     y = alpha * x / (1 + x * (alpha - 1))
#     f = -1 / (y - x)
#     return f
#
# def n(f):
#     n_list = []
#     no = bD.InitialMoles
#     n_list.append(no)
#     for i,f in enumerate(f):
#         n = no*math.exp(-f)
#         n_list.append(n)
#     return n_list
#
# print(bD.LiquidFractionList(mix))
# F = Trapezoidal_list(f,bD.LiquidFractionList(mix))
# FF = Simpson_list(f,bD.LiquidFractionList(mix))
# print(n(F))
# print(bD.LiquidProfile(mix))
# plt.plot(bD.LiquidFractionList(mix), n(FF), bD.LiquidFractionList(mix), bD.LiquidProfile(mix))
# plt.show()

#endregion

#region PLOTS

mixx = thermoFlowStream(['heptane', 'hexane'], zi=[0.5, 0.5], Nf=0.0, T=298.0)
mixx2 = thermoFlowStream(['toluene', 'benzene'], zi=[0.5, 0.5], Nf=0.0, T=298.0)
mixx3 = thermoFlowStream(['ethanol', 'water'], zi=[0.5, 0.5], Nf=0.0, T=298.0)
mixx4 = thermoFlowStream(['ethanol', '1-propanol'], zi=[0.5, 0.5], Nf=0.0, T=298.0)
mixx5 = thermoFlowStream(['ethanol', '2-propanol'], zi=[0.5, 0.5], Nf=0.0, T=298.0)
mixx6 = thermoFlowStream(['ethanol', 'butanol'], zi=[0.5, 0.5], Nf=0.0, T=298.0)
mixx7 = thermoFlowStream(['methanol', 'water'], zi=[0.5, 0.5], Nf=0.0, T=298.0)
mixxx = thermoFlowStream(['chloroform', 'ether', 'ethanol'], zi=[0.5, 0.2, 0.3], Nf=0.0, T=298.0)

def Txy(binary_mixture, P=101325, x1=True):
    binary_mixture.Pressure = P
    if x1 is True:
        j = 0
    else:
        j = 1

    Tsat1 = binary_mixture.Tsats[j]
    Tsat2 = binary_mixture.Tsats[1-j]

    T_list = np.linspace(Tsat1, Tsat2, 100)
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
        x1.append((P - Psat2[i]) / (Psat1[i] - Psat2[i]))
        y1.append(x1[i] * binary_mixture.gammaUNIFAC[1-j] * Psat1[i] / binary_mixture.Pressure)

    plt.rc('font', family='Century Gothic')
    plt.plot(x1, T_list, y1, T_list)
    plt.xlabel('x, y')
    plt.ylabel('T (K)')
    plt.show()

def yx(binary_mixture, P=101325, x1=True):
    binary_mixture.Pressure = P
    if x1 is True:
        j = 0
    else:
        j = 1

    Tsat1 = binary_mixture.Tsats[j]
    Tsat2 = binary_mixture.Tsats[1-j]

    T_list = np.linspace(Tsat1, Tsat2, 100)
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
        x1.append((P - Psat2[i]) / (Psat1[i] - Psat2[i]))
        y1.append(x1[i] * binary_mixture.gammaUNIFAC[1-j] * Psat1[i] / binary_mixture.Pressure)

    plt.rc('font', family='Century Gothic')
    plt.plot(x1,y1,x1,x1)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

def alpha(mixture, P=101325):
    x_list = np.linspace(0, 1, 100)
    n = len(mixture.Components) - 1
    mixture.Pressure = P
    Psat_list = mixture.Psats
    gamma_list = []
    alpha_list = []
    for i, x in enumerate(x_list):
        mixture._zi = []
        mixture._zi.append(x)
        for elements in range(n):
            mixture._zi.append((1-x)/n)
        gamma_list.append(mixture.gammaUNIFAC)
        alpha = []
        for j, Psat in enumerate(Psat_list):
            alpha.append(Psat * gamma_list[i][j] / mixture.Pressure)
        alpha_list.append(list(alpha))

    plt.rc('font', family='Century Gothic')
    plt.plot(x_list, alpha_list)
    plt.xlabel('x')
    plt.ylabel('Relative Volatility')
    plt.legend(mixture.Components)
    plt.show()

#alpha(mixxx)
#yx(mixx2)
#Txy(mixx2)

#endregion
