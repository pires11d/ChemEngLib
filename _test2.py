from Chemical import *
from Converter import *
from matplotlib import animation
import matplotlib.pyplot as plt

#region COMPONENTS

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

n = Mixture([s,o])
n.wi = [0.8,0.2]
n.V = 0.1

s = Stream([s,o])
s.wi = [1.0,0.0]
s.Vf = VolumeFlow(25).m3_h

#endregion

#region TANKS

from UnitOp import *
tk = recTank(1.0,1.0,1.5)
#tk = cylTank(1.5,1.5,30)
tk.Mixture = m
tk.Inlets = [s]
tk.OutletVolumeFlow = VolumeFlow(50).m3_h

h = Hopper(initial_angle=30,final_angle=40,Hmin=0.8,Hmax=1.0,r1=3.0,r2=5.5,R=6.0)
h.Mixture = n
h.X = tk.Width * 2

p = Pipe(0.1,10.0)
p.From = tk
p.To = h

#endregion

#region ANIMATION

# Canvas #
fig = plt.figure()
fig.set_dpi(100)
fig.set_size_inches(8,8)
ax = plt.axes(xlim=(-0.5, 5.5), ylim=(-0.5, 4.5))
# Function Definition #
def init():
    ax.add_patch(tk.DrawContour)
    ax.add_patch(h.DrawContour)
    return tk.DrawContour,h.DrawContour,
def animate(i):
    tk.NextTime
    p.InletStream = tk.OutletStream
    h.Inlets = [p.OutletStream]
    h.NextTime
    patches = []
    patches.append(tk.DrawLiquid)
    patches.append(h.DrawLiquid)
    patches.append(h.DrawTopArrow)
    patches.append(h.DrawBottomArrow)
    patches.append(tk.DrawTopArrow)
    patches.append(tk.DrawBottomArrow)
    patches.append(p.DrawContour)
    for patch in patches:
        ax.add_patch(patch)
    return patches
# Function Call #
anim = animation.FuncAnimation(fig,animate,init_func=init,frames=1000,interval=30,blit=True)
plt.show()

#endregion
