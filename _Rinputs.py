from Chemicals import *
from Tools import *

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

