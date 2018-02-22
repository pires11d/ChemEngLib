from Chemical import *
from matplotlib import animation
import matplotlib.pyplot as plt

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

st1 = Stream([s,o])
st1.wi = [1.0,0.0]
st1.Vf = 1.0
st2 = Stream([s,o])
st2.wi = [0.5,0.5]
st2.Vf = 2.0


from UnitOp import Hopper, Mixer
h = Hopper(100)
h.Volume = h.MaxVolume / 2
h.Inlets = [st1,st2]
h.OutletVolumeFlow = 3.5
#print(h.OutletStream().MassFractions)


# #region ANIMATION

# # Canvas #
# fig = plt.figure()
# fig.set_dpi(100)
# fig.set_size_inches(5,5)
# ax = plt.axes(xlim=(-1, 1), ylim=(0, 2))

# def init():
#     ax.add_patch(h.Contour())
#     return h.Contour(),

# def animate(i):
#     h.NextVolume()
#     patch = h.Liquid()
#     ax.add_patch(patch)
#     return patch,

# # Function Call #
# anim = animation.FuncAnimation(fig,animate,init_func=init,frames=1000,interval=30,blit=True)
# plt.show()

# #endregion
