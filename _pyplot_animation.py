import matplotlib.pyplot as plt
from matplotlib import animation

# Canvas #
fig = plt.figure()
fig.set_dpi(100)
fig.set_size_inches(5,5)

ax = plt.axes(xlim=(0, 10), ylim=(0, 10))
patch = plt.Circle((5, -5), 0.75, fc='y')

def init():
    patch.center = (5, 5)
    ax.add_patch(patch)
    return patch,

def animate(i):
    x, y = patch.center
    x = 5 + 0.01 * i
    y = 5 + 0.01 * i
    patch.center = (x, y)
    return patch,

# Animation #
anim = animation.FuncAnimation(fig,animate,init_func=init,frames=1000,interval=10,blit=True)
plt.show()