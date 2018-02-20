from Chemical import Substance, Mixture


s = Substance('solvent')
s.Density0 = 680
s.Viscosity0 = 1e-3
s.MolarMass = 40

o = Substance('oil')
o.Density0 = 900
o.Viscosity0 = 1e-5
o.MolarMass = 100

m = Mixture([s,o])
m.wi = [0.5,0.5]
# print(m.MolarMasses)
# print(m.MassFractions())
# print(m.MolarFractions())
