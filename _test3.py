from Chemical import *
from Converter import *
from UnitOp import *

g = Substance('air')
g.Phase = 'g'
g.Density0 = 1.0
g.Viscosity0 = 1e-5
g.SpecificHeat0 = 1
g.MolarMass = 27e-3

s = Substance('solvent')
s.Phase = 'l'
s.Density0 = 680.0
s.Viscosity0 = 1e-3
s.SpecificHeat0 = 2
s.MolarMass = 40e-3

o = Substance('oil')
o.Phase = 'l'
o.Density0 = 900.0
o.Viscosity0 = 1e-5
o.SpecificHeat0 = 3
o.MolarMass = 100e-3

i = Substance('inert')
i.Phase = 's'
i.Density0 = 900.0
i.SpecificHeat0 = 3
i.MolarMass = 100e-3


m = Mixture([s,o,i])
m.wi = [0.2,0.7,0.1]
m.V = 1.0

st = Stream([g,i])
st.wi = [0.8,0.2]
st.Mf = 1.0
st.ParticleSize = 1e-3

st2 = Stream([s,i])
st2.wi = st.wi
st2.Mf = 1.0
st2.ParticleSize = 1e-3

c = Cyclone()
print(c.OverFlow(st),c.UnderFlow(st))

hc = Hydrocyclone(d50=1e-3)
print(hc.OverFlow(st2),hc.UnderFlow(st2))

mm = Mill(power=10000)
print(mm.OutletParticleSize(st,1.5,10))