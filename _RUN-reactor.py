from Chemicals import *
from UnitOps.Tanks import Reactor
from _Rinputs import *
import matplotlib.pyplot as plt

subs = [a,b,c,d]

m = Mixture(subs)
m.zi = [0.5,0.5,0,0]
m.N = 100

r1 = Reaction(m, order=[1,0])
r1.Stoichiometry = [1,2,1,2]
r1.NumberOfReactants = 2
r1.NumberOfProducts = 2
r1.KineticConstant = 1e-1
# print(r1.Order)

R = Reactor(r1)
R.dt = 0.5

# print(r1.ReactantNames)
# print(r1.ConsumptionRate(a))
# print(r1.ConsumptionRate(b))
# print(r1.ProductNames)
# print(r1.ProductionRate(c))
# print(r1.ProductionRate(d))

t_list = []
C_list = []
time = 0
while time < 10:
    time += R.dt
    t_list.append(time)
    C_list.append(R.Mixture.MolarFractions)
    R.NextMoles

plt.plot(t_list,C_list)
plt.show()