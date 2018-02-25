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
s.Vf = VolumeFlow(0).m3_h

#endregion

#region TANKS

from UnitOp import *

tk = recTank(1.0,1.0,1.5)
tk.Mixture = m
tk.OutletVolumeFlow = VolumeFlow(50).m3_h

h = Hopper(initial_angle=30,final_angle=40,Hmin=0.8,Hmax=1.0,r1=3.0,r2=5.5,R=6.0)
h.Mixture = n
h.OutletVolumeFlow = VolumeFlow(20).m3_h
h.X = tk.Width * 1.5

ctk = cylTank(1.2,1.2,30)
ctk.Mixture = n
ctk.X = h.X + h.Width * 1.5
ctk.OutletVolumeFlow = 0.0

p1 = Pipe(0.1,10.0)
p2 = Pipe(0.1,10.0)
sh1 = Shower()
sh2 = Shower()

p1.From = tk
p1.To = sh1
sh1.From = p1
sh1.To = h
p2.From = h
p2.To = sh2
sh2.From = p2
sh2.To = ctk

# f = Flash()
# f.X = h.X
# f.Y = h.Height * 2

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
    ax.add_patch(ctk.DrawContour)
    ax.add_patch(sh1.DrawContour)
    ax.add_patch(sh2.DrawContour)
    return tk.DrawContour,h.DrawContour,ctk.DrawContour, sh1.DrawContour, sh2.DrawContour

def animate(i):
    if i > 100:
        tk.Inlets.append(ctk.Outlet)
        ctk.OutletVolumeFlow = VolumeFlow(20).m3_h
    # Inlets and Outlets
    tk.Inlets = [s]
    tk.Inlets = [ctk.Outlet]
    p1.Inlet = tk.Outlet
    sh1.Inlet = p1.Outlet
    h.Inlets = [sh1.Outlet]
    p2.Inlet = h.Outlet
    sh2.Inlet = p2.Outlet
    ctk.Inlets = [sh2.Outlet]
    # NextTime function
    tk.NextTime
    h.NextTime
    ctk.NextTime
    # Drawings
    patches = []
    patches.append(tk.DrawLiquid)
    patches.append(tk.DrawTopArrow)
    patches.append(tk.DrawBottomArrow)
    patches.append(h.DrawLiquid)
    patches.append(h.DrawTopArrow)
    patches.append(h.DrawBottomArrow)
    patches.append(ctk.DrawLiquid)
    patches.append(ctk.DrawTopArrow)
    patches.append(ctk.DrawBottomArrow)
    patches.append(p1.DrawContour)
    patches.append(p2.DrawContour)
    patches.append(sh1.DrawLiquid)
    patches.append(sh1.DrawStream)
    patches.append(sh2.DrawLiquid)
    patches.append(sh2.DrawStream)
    for patch in patches:
        ax.add_patch(patch)
    return patches
    #print(tk.NextVolume,h.NextVolume,ctk.NextVolume)

# Function Call #
anim = animation.FuncAnimation(fig,animate,init_func=init,frames=1000,interval=30,blit=True)
plt.show()

#endregion
