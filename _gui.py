import time
from Converter import *
from Correlation import *
from Numerical import *
from Material import *
from UnitOp import *


water = FlowStream(['water'], wi=[1.0], Mf=1.0, T=25+273)


#region INTERFACE

from tkinter import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure


# Initializing Parent Interface #
root = Tk()
root.title("cylTank Simulation v1.0")

#region LABELS
title = Label(root, text="cylTank Simulation", font=("Century Gothic", 18))
title.grid(row=0, columnspan=8)
empty = Label(root, text="        ")
empty.grid(row=1, column=3)

label_GP = Label(root, text="Geometry Parameters:", font=("Century Gothic", 11, "bold")).grid(row=2, columnspan=2, sticky=NW)
label_D = Label(root, text="cylTank's Diameter", font=("Century Gothic", 10)).grid(row=2, column=0, sticky=W)
label_H = Label(root, text="cylTank's Height", font=("Century Gothic", 10)).grid(row=3, column=0, sticky=W)
label_theta = Label(root, text="Cone angle", font=("Century Gothic", 10)).grid(row=4, column=0, sticky=W)
ulabel_D = Label(root, text="m", font=("Century Gothic", 10)).grid(row=2, column=2)
ulabel_H = Label(root, text="m", font=("Century Gothic", 10)).grid(row=3, column=2)
ulabel_theta = Label(root, text="°", font=("Century Gothic", 10)).grid(row=4, column=2)

label_OP = Label(root, text="Operation Parameters:", font=("Century Gothic", 11, "bold")).grid(row=6, columnspan=2, sticky=NW)
label_tf = Label(root, text="Timelapse", font=("Century Gothic", 10)).grid(row=6, column=0, sticky=W)
label_Vo = Label(root, text="Initial Volume", font=("Century Gothic", 10)).grid(row=7, column=0, sticky=W)
label_ofr = Label(root, text="Outlet Flow Ratio", font=("Century Gothic", 10)).grid(row=8, column=0, sticky=W)
ulabel_tf = Label(root, text="s", font=("Century Gothic", 10)).grid(row=6, column=2)
ulabel_Vo = Label(root, text="m³", font=("Century Gothic", 10)).grid(row=7, column=2)
ulabel_ofr = Label(root, text="m³/m³", font=("Century Gothic", 10)).grid(row=8, column=2)
#endregion

#region ENTRIES
sv1 = StringVar()
sv1.trace("w", lambda name, index, mode, sv=sv1: onChange(sv1))
sv2 = StringVar()
sv2.trace("w", lambda name, index, mode, sv=sv2: onChange(sv2))
sv3 = StringVar()
sv3.trace("w", lambda name, index, mode, sv=sv3: onChange(sv3))

tb_D = Entry(root, textvariable=sv1)
tb_D.grid(row=2, column=1)
tb_H = Entry(root, textvariable=sv2)
tb_H.grid(row=3, column=1)
tb_angle = Entry(root, textvariable=sv3)
tb_angle.grid(row=4, column=1)
tb_tf = Entry(root)
tb_tf.grid(row=6, column=1)
tb_Vo = Entry(root)
tb_Vo.grid(row=7, column=1)
tb_ofr = Entry(root)
tb_ofr.grid(row=8, column=1)
#endregion

#region INITIAL DATA
tanque = cylTank(0.5, 1.0, 0.0)
D = tanque.D
Hc = tanque.ConeHeight
Ht = tanque.TankHeight
maximum = max(tanque.D, Ht)
#endregion

#region GRAPHS

class Graph(Frame):
    def __init__(self, parent):
        Frame.__init__(self, parent)
        fl = Figure(figsize=(5, 5), dpi=100)

        #region DEFINITIONS

        tb_D.insert(0, tanque.D)
        tb_H.insert(0, tanque.H)
        tb_angle.insert(0, tanque.Theta)
        tb_tf.insert(0, 100)
        tb_Vo.insert(0, tanque.InitialVolume)
        tb_ofr.insert(0, 0.0)

        #endregion

        #region TANK PLOT

        xTk = [0, -D / 2, -D / 2, +D / 2, +D / 2, 0]
        yTk = [0, Hc, Ht, Ht, Hc, 0]
        ax = fl.add_subplot(111)
        ax.plot(xTk, yTk, color='black')
        ax.set_xlim(-maximum * 0.75, +maximum * 0.75)
        ax.set_ylim(0, maximum * 1.5)

        #endregion

        canvas = FigureCanvasTkAgg(fl, self)
        canvas.show()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=BOTTOM, fill=BOTH, expand=True)

class DynamicGraph(Frame):
    def __init__(self, parent):
        Frame.__init__(self, parent)
        fl = Figure(figsize=(5, 5), dpi=100)

        #region DEFINITIONS

        tanque.D = float(tb_D.get())
        tanque.H = float(tb_H.get())
        tanque.Theta = float(tb_angle.get())

        D = tanque.D
        Hc = tanque.ConeHeight
        Ht = tanque.TankHeight
        Hs = tanque.TankHeight * 1.3
        maximum = max(tanque.D, Ht)
        h_list = tanque.LiquidHeights([water], outlet_list, t_list, return_array=True)
        d_list = tanque.LiquidDiameters([water], outlet_list, t_list, return_array=True)
        x_list = [tuple([-d / 2, d / 2]) for d in d_list]
        y_list = [tuple([h, h]) for h in h_list]

        #endregion

        # region STATIC PLOTS

        # cylTank #
        xTk = [0, -D / 2, -D / 2, +D / 2, +D / 2, 0]
        yTk = [0, Hc, Ht, Ht, Hc, 0]
        ax = fl.add_subplot(111)
        ax.plot(xTk, yTk, color='black')
        ax.set_xlim(-maximum * 0.75, +maximum * 0.75)
        ax.set_ylim(0, maximum * 1.5)
        # Shower #
        xS = [0]
        yS = [Hs]
        ax.scatter(xS, yS, s=100.0, c='gray')
        endregion

        #region DYNAMIC PLOTS

        for i, t in enumerate(t_list):
             # liquid level #
             if tanque.theta > 0.0:
                 ax.plot(x_list[i], y_list[i], color='cyan')
             else:
                 ax.stackplot(x_list[i], y_list[i], color='cyan')
            # droplets #
             xd = [0, 0, 0, 0, 0]
             xd1 = np.linspace(0, -0.05, 5)
             xd2 = np.linspace(0, +0.05, 5)
             yd = [0.1 * ht + 0.9 * hs, 0.3 * ht + 0.7 * hs, 0.5 * ht + 0.5 * hs, 0.7 * ht + 0.3 * hs,
                   0.9 * ht + 0.1 * hs]
             yyd = [hs, 0.2 * ht + 0.8 * hs, 0.4 * ht + 0.6 * hs, 0.6 * ht + 0.4 * hs, 0.8 * ht + 0.2 * hs]
             a = ax.scatter(xd, yd, s=20.0, c='cyan')
             a1 = ax.scatter(xd1, yd, s=10.0, c='cyan')
             a2 = ax.scatter(xd2, yd, s=10.0, c='cyan')
             ax.draw()
             ax.pause(0.001)
             a.remove(), a1.remove(), a2.remove()
             aa = ax.scatter(xd, yyd, s=20.0, c='cyan')
             aa1 = ax.scatter(xd1, yyd, s=10.0, c='cyan')
             aa2 = ax.scatter(xd2, yyd, s=10.0, c='cyan')
             ax.draw()
             ax.pause(0.001)
             aa.remove(), aa1.remove(), aa2.remove()

        canvas = FigureCanvasTkAgg(fl, self)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True)

        #endregion

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=BOTTOM, fill=BOTH, expand=True)

#endregion

#region FUNCTIONS:

def onChange(sv):
    DynamicGraph(root).grid(row=2, column=4, rowspan=8)

def Run(event):
    tanque.InitialVolume = float(tb_Vo.get())
    ofr = float(tb_ofr.get())
    t_list = np.arange(0, float(tb_tf.get()), 1)
    outlet_list = [ofr for t in t_list]

    tanque.Animation([water], outlet_list, t_list)
    plt.close()

#endregion

Graph(root).grid(row=2, column=4, rowspan=8)
button_1 = Button(root, text="Run Simulation!")
button_1.bind("<Button-1>", Run)
button_1.grid(row=9, columnspan=3)

root.mainloop()

#endregion