from Chemical import Substance, Mixture, Stream


s = Substance('solvent')
s.Density0 = 680.0
s.Viscosity0 = 1e-3
s.SpecificHeat0 = 2
s.MolarMass = 40e-3

o = Substance('oil')
o.Density0 = 900.0
o.Viscosity0 = 1e-5
o.SpecificHeat0 = 3
o.MolarMass = 100e-3

m = Mixture([s,o])
m.wi = [0.5,0.5]

st1 = Stream(m)
st1.Mf = 1.0
st2 = Stream(m)
st2.Mf = 2.0

from UnitOp import Hopper, Mixer

dt = 0.1
h = Hopper(100)

while 1<2:
    h.Inlets = [st1,st2]
    h.OutletVolumeFlow = 0.1
    h.Draw(dt)
