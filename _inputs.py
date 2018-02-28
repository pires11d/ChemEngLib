from Chemical import *
from Converter import *
from UnitOp import *


#region SUBSTANCES:

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

#endregion

#region MIXTURES:

comp = [o,s]

m0 = Mixture(comp)
m0.wi = [0.8,0.2]
# m0.wi = [1]
m0.V = 1
m0.Temperature = 373

n0 = Mixture(comp)
n0.wi = [0.2,0.8]
# n0.wi = [1]
n0.V = 0.1
n0.Temperature = 273

m = Mixture(comp)
m.N = 0.0
m.wi = [0.5,0.5]
m.Ki = [0.2,1]

n = Mixture(comp)
n.N = 100
n.zi = [0.5,0.5]
n.Ki = [0.2,1] 

#endregion

#region STREAMS:

st0 = Stream(comp)
st0.wi = [0.0,1.0]
# st0.wi = [1]
st0.Vf = VolumeFlow(0).m3_h
st0.Temperature = 373

st = Stream([g,i])
st.wi = [0.8,0.2]
st.Mf = 1
st.ParticleSize = 1e-3

st2 = Stream([s,i])
st2.wi = st.wi
st2.Mf = 1
st2.ParticleSize = 1e-3

F = Stream(comp)
F.wi = [0.2,0.8]
F.Nf = 100
F.Ki = [10,0.2]

#endregion