from Chemical import *
from Converter import *
from matplotlib import animation
import matplotlib.pyplot as plt

#region STREAMS

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
m.wi = [0.2,0.8]
m.V = 1.0
m.Temperature = 300

s = Stream([s,o])
s.wi = [0.9,0.1]
s.Vf = VolumeFlow(30.0).m3_h
s.Temperature = 320

#endregion

#region TANK

from UnitOp import *
tk = recTank(1.0,1.0,1.5)
#tk = cylTank(1.5,1.5,30)
tk.Mixture = m
tk.Inlets = [s]
tk.OutletVolumeFlow = VolumeFlow(50.0).m3_h

h = Hopper(initial_angle=30,final_angle=60,Hmin=0.8,Hmax=1.0,r1=3.0,r2=5.5,R=6.0)
h.Mixture = m
h.X = tk.Width * 2

#endregion

#region ANIMATION

# Canvas #
fig = plt.figure()
fig.set_dpi(100)
fig.set_size_inches(8,8)
ax = plt.axes(xlim=(-0.5, 8.0), ylim=(-0.5, 2.5))
# Function Definition #
def init():
    ax.add_patch(tk.DrawContour)
    ax.add_patch(h.DrawContour)
    return tk.DrawContour,h.DrawContour,
def animate(i):
    h.Inlets = [tk.OutletStream]
    tk.NextTime
    h.NextTime
    patch1 = tk.DrawLiquid
    patch2 = h.DrawLiquid
    ax.add_patch(patch1)
    ax.add_patch(patch2)
    return patch1,patch2,
# Function Call #
anim = animation.FuncAnimation(fig,animate,init_func=init,frames=1000,interval=30,blit=True)
plt.show()

#endregion
