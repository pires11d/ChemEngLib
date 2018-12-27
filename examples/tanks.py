from unit_ops import *
from tools import *
from .inputs import *
import matplotlib.pyplot as plt


#region INPUTS
tk = recTank(1.0,1.0,1.5)
tk.Mixture = m0
tk.OutletVolumeFlow = VolumeFlow(100).m3_h

h = Hopper(initial_angle=30,final_angle=40,Hmin=0.8,Hmax=1.0,r1=3.0,r2=5.5,R=6.0)
h.Mixture = n0
h.OutletVolumeFlow = VolumeFlow(50).m3_h
h.X = tk.Width * 1.5

ctk = cylTank(1.2,1.2,30)
ctk.Mixture = n0
ctk.X = h.X + h.Width * 1.5
ctk.OutletVolumeFlow = 0.0

p1 = Pipe(0.1,10.0)
p2 = Pipe(0.1,10.0)
sh1 = Shower()
sh2 = Shower()

p1.FromBottom = tk
p1.ToTop = sh1
sh1.From = p1
sh1.To = h
p2.FromBottom = h
p2.ToTop = sh2
sh2.From = p2
sh2.To = ctk
#endregion

#region ANIMATION

# Canvas
fig = plt.figure()
fig.set_dpi(100)
fig.set_size_inches(5,3)
ax = plt.axes(xlim=(-0.5, 5.5), ylim=(-0.5, 4.5))

def init():
    ax.add_patch(tk.DrawContour)
    ax.add_patch(h.DrawContour)
    ax.add_patch(ctk.DrawContour)
    ax.add_patch(sh1.DrawContour)
    ax.add_patch(sh2.DrawContour)
    return tk.DrawContour,h.DrawContour,ctk.DrawContour, sh1.DrawContour, sh2.DrawContour

def animate(i):
    if i > 50:
        ctk.OutletVolumeFlow = VolumeFlow(40).m3_h
        tk.Inlets = [ctk.Outlet]
    else:
        tk.Inlets = [st0]
    # Inlets and Outlets
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
#endregion

# Function Call
def show_animation():
    anim = animation.FuncAnimation(fig,animate,init_func=init,frames=1000,interval=30,blit=True)
    plt.show()