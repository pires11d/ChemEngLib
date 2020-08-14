from equipment import *
from .inputs import *


#region INPUTS
fl = Flash()
fl.Y = 1.0
p1 = Pipe(0.1,1.0)
p2 = Pipe(0.1,1.0)

t1 = Tank()
t1.Mixture = m
t1.X = fl.X + fl.Width*4
t1.Y = fl.Y + fl.Height
t2 = Tank()
t2.Mixture = m
t2.X = fl.X + fl.Width*4

p1.FromTop = fl
p1.ToTop = t1
p2.FromBottom = fl
p2.ToTop = t2
#endregion

#region ANIMATION:

# Canvas
fig = plt.figure()
fig.set_dpi(100)
fig.set_size_inches(8,8)
ax = plt.axes(xlim=(-0.5, 5.5), ylim=(-0.5, 4.5))

# Animation Function
def init():
    ax.add_patch(fl.DrawContour)
    ax.add_patch(t1.DrawContour)
    ax.add_patch(t2.DrawContour)
    return fl.DrawContour, t1.DrawContour, t2.DrawContour

def draw(DrawObjects):
    for o in DrawObjects:
        ax.add_patch(o)

def animate(i):
    dwi = 0.01*math.sin(i/5)
    F.wi = [F.wi[0]+dwi, F.wi[1]-dwi]
    # Inlets and Outlets
    fl.Inlet = F
    p1.Inlet = fl.VaporOutlet
    p2.Inlet = fl.LiquidOutlet
    t1.Inlets = [p1.Outlet]
    t2.Inlets = [p2.Outlet]
    # print(fl.VaporFlow, fl.LiquidFlow)
    # print("V: "+str(fl.yi))

    # NextTime function
    t1.NextTime
    t2.NextTime

    # Drawings
    patches = []
    patches.append(fl.DrawLiquid)
    patches.append(fl.DrawVapor)
    patches.append(fl.DrawInletArrow)
    patches.append(fl.DrawTopArrow)
    patches.append(fl.DrawBottomArrow)
    patches.append(p1.DrawContour)
    patches.append(p2.DrawContour)
    patches.append(t1.DrawLiquid)
    patches.append(t2.DrawLiquid)
    for patch in patches:
        ax.add_patch(patch)
    return patches

#endregion

# Function Call
def show_animation():
    anim = animation.FuncAnimation(fig,animate,init_func=init,frames=1000,interval=30,blit=True)
    plt.show()