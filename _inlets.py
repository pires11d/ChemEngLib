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

comp = [s,o]

m = Mixture(comp)
m.V = 0
m.wi = [0.5,0.5]
m.Ki = [10,0.2]

#endregion

#region STREAMS:

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