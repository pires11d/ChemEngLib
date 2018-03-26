from Chemicals import *
from UnitOps.Tanks import Reactor
from _Rinputs import *
import matplotlib.pyplot as plt

t = Reactor()
t.dt = 1
t.Mixture = m

r = Reaction(subs)
r.NumberOfReactants = 2
r.NumberOfProducts = 2
r.Tank = t
r.KineticConstant = 0.2

# print(r.ReactantNames, r.ProductNames)
t_list = []
z_list = []
time = 0
while time < 10:
    time += t.dt
    t_list.append(time)
    z_list.append(t.Mixture.zi)
    dC = r.ConsumptionRate("A")
    Na = m.MolarFractions[0] * m.Moles
    Nc = m.MolarFractions[2] * m.Moles
    Na = Na - dC * t.dt
    Nc = Nc + dC * t.dt
    za = Na/m.Moles
    zc = Nc/m.Moles
    t.Mixture.zi = [za, 0, zc, 0]

plt.plot(t_list,z_list)
plt.show()