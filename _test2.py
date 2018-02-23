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
m.V = 0.1

s = Stream([s,o])
s.wi = [1.0,0.0]
s.Vf = VolumeFlow(120.0).m3_h

#endregion

#region TANK

from UnitOp import *
#tk = recTank(1.0,1.0,1.5)
tk = cylTank(2.0,2.0,30)
tk.Mixture = m
tk.Inlets = [s]
tk.OutletVolumeFlow = VolumeFlow(100.0).m3_h

h = Hopper(initial_angle=30,final_angle=60,Hmin=0.8,Hmax=1.0,r1=3.0,r2=5.5,R=6.0)
h.Mixture = m
h.Inlets = [tk.OutletStream]

#endregion

#region ANIMATION

# Canvas #
fig = plt.figure()
fig.set_dpi(100)
fig.set_size_inches(8,8)
ax = plt.axes(xlim=(-0.5, 10), ylim=(-0.5, 5))
def init():
    ax.add_patch(tk.DrawContour)
    return tk.DrawContour,
def animate(i):
    tk.NextTime
    print(tk.NextVolume)
    # print(tk.OutletStream.MassFractions)
    # print(tk.NextMixture.MassFractions)
    patch = tk.DrawLiquid
    ax.add_patch(patch)
    return patch,
# Function Call #
anim = animation.FuncAnimation(fig,animate,init_func=init,frames=1000,interval=30,blit=True)
plt.show()

#endregion
