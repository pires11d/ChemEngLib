from chemical import *
from equipment.Tanks import Reactor
import matplotlib.pyplot as plt


#region INPUTS
a = Substance('A')
a.Phase = 'l'
a.Density0 = 900.0
a.Viscosity0 = 1e-5
a.SpecificHeat0 = 3
a.MolarMass = 100e-3

b = Substance("B")
b.Phase = a.Phase
b.Density0 = a.Density0
b.Viscosity0 = a.Viscosity0
b.SpecificHeat0 = a.SpecificHeat0
b.MolarMass = a.MolarMass

c = Substance("C")
c.Phase = a.Phase
c.Density0 = a.Density0
c.Viscosity0 = a.Viscosity0
c.SpecificHeat0 = a.SpecificHeat0
c.MolarMass = a.MolarMass

d = Substance("D")
d.Phase = a.Phase
d.Density0 = a.Density0
d.Viscosity0 = a.Viscosity0
d.SpecificHeat0 = a.SpecificHeat0
d.MolarMass = a.MolarMass

subs = [a,b,c,d]

m = Mixture(subs)
m.zi = [0.5,0.5,0,0]
m.N = 100

r1 = Reaction(m, order=[1,0])
r1.Stoichiometry = [1,2,1,2]
r1.NumberOfReactants = 2
r1.NumberOfProducts = 2
r1.KineticConstant = 1e-1

R = Reactor(r1)
R.dt = 0.5
#endregion

#region PLOTS
def plot(t=10):
    t_list = []
    C_list = []
    time = 0
    while time < t:
        time += R.dt
        t_list.append(time)
        C_list.append(R.Mixture.MolarFractions)
        R.NextMoles

    plt.plot(t_list,C_list)
    plt.show()
#endregion