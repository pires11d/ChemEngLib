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
m.Temperature = 300

s = Stream([s,o])
s.wi = [1.0,0.0]
s.Vf = VolumeFlow(10.0).m3_h
s.Temperature = 320

#endregion

#region TANKS

from UnitOp import *
tk = recTank(1.0,1.0,1.5)
#tk = cylTank(1.5,1.5,30)
tk.Mixture = m
tk.Inlets = [s]
tk.OutletVolumeFlow = VolumeFlow(50.0).m3_h

h = Hopper(initial_angle=30,final_angle=40,Hmin=0.8,Hmax=1.0,r1=3.0,r2=5.5,R=6.0)
h.Mixture = m
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
    patch_list = []
    p1 = tk.DrawLiquid
    p2 = h.DrawLiquid
    p3 = h.DrawTopArrow
    p4 = h.DrawBottomArrow
    p5 = tk.DrawTopArrow
    p6 = tk.DrawBottomArrow
    p7 = p.DrawContour
    patch_list.append(p1)
    patch_list.append(p2)
    patch_list.append(p3)
    patch_list.append(p4)
    patch_list.append(p5)
    patch_list.append(p6)
    patch_list.append(p7)
    for patch in patch_list:
        ax.add_patch(patch)
    return p1,p2,p3,p4,p5,p6,p7,
# Function Call #
anim = animation.FuncAnimation(fig,animate,init_func=init,frames=1000,interval=30,blit=True)
plt.show()

#endregion
