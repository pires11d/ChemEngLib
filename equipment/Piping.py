import math
import matplotlib.pyplot as plt
from matplotlib import animation
from collections import Counter, defaultdict, OrderedDict
from chemical import *
from utilities.Geometry import Circle
from utilities.Correlation import Re_D


class Pipe:
    def __init__(self, diameter, length):
        self.Name = None
        self.Inlet = None
        self.Diameter = diameter
        self.Length = length
        self.InletPressure = 101325.0
        self.FromBottom = None
        self.FromTop = None
        self.FromSide = None
        self.ToBottom = None
        self.ToTop = None
        self.ToSide = None

    @property
    def Height(self):
        # Begin:
        if self.FromBottom != None:
            y0 = self.FromBottom.Y
        elif self.FromTop != None:
            y0 = self.FromTop.Y + self.FromTop.Height
        else:
            y0 = self.FromSide.Y + self.FromSide.Height/2
        # End:
        if self.ToBottom != None:
            y1 = self.ToBottom.Y
        elif self.ToTop != None:
            y1 = self.ToTop.Y + self.ToTop.Height
        else:
            y1 = self.ToSide.Y + self.ToSide.Height/2
        # Difference:
        return y1 - y0

    @property
    def Outlet(self):
        return self.Inlet

    @property
    def InletVelocity(self):
        return self.Inlet.VolumeFlow / Circle(self.Diameter).Area
    
    @property
    def HeadLoss(self):
        v = self.InletVelocity
        L = self.Length
        D = self.Diameter
        f = 64/Re_D(self.Inlet,D)
        g = 9.81
        return (f * L * v**2)/(2 * g * D)

    @property
    def OutletVelocity(self):
        vo = self.InletVelocity
        g = 9.81
        dz = self.Height
        lwf = self.HeadLoss
        return math.sqrt(vo**2 - 2*g*(dz + lwf))

    @property
    def OutletPressure(self):
        Po = self.InletPressure
        g = 9.81
        dz = self.Height
        lwf = self.HeadLoss
        return Po - self.Inlet.Density * g * (dz + lwf)

    @property
    def DrawContour(self):
        # Assigning 'From' connection:
        if self.FromBottom != None:
            xFrom = self.FromBottom.X + self.FromBottom.Width/2
            yFrom = self.FromBottom.Y
            x1 = xFrom
            y1 = yFrom - 0.25
        elif self.FromTop != None:
            xFrom = self.FromTop.X + self.FromTop.Width/2
            yFrom = self.FromTop.Y + self.FromTop.Height
            x1 = xFrom
            y1 = yFrom + 0.25
        else:
            xFrom = self.FromSide.X
            yFrom = self.FromSide.Y + self.FromSide.Height/2
            x1 = xFrom + 0.25
            y1 = yFrom
        # Assigning 'To' connection:
        if self.ToBottom != None:
            xTo = self.ToBottom.X + self.ToBottom.Width/2
            yTo = self.ToBottom.Y
            x2 = xTo
            y2 = yTo - 0.25
        elif self.ToTop != None:
            xTo = self.ToTop.X + self.ToTop.Width/2
            yTo = self.ToTop.Y + self.ToTop.Height
            x2 = xTo
            y2 = yTo + 0.25
        else:
            xTo = self.ToSide.X
            yTo = self.ToSide.Y + self.ToSide.Height/2
            x2 = xTo - 0.25
            y2 = yTo
        xmid = (xFrom+xTo)/2
       # Creating Line
        points = [[xFrom,yFrom],[x1,y1],[xmid,y1],[xmid,y2],[x2,y2],[xTo,yTo]]
        patch = plt.Polygon(points, closed=None, fill=None, lw=1, edgecolor='black')
        return patch


class Shower:
    def __init__(self, diameter=0.1, length=1.0, height = 0.0):
        self.Name = None
        self.Inlet = None
        self.Diameter = diameter
        self.Length = length
        self.Height = height
        self.GeometryFactor = 1.0
        self.From = None
        self.To = None
        # Drawing #
        self.dt = 0.1
        self.Width = 0
    
    @property
    def X(self):
        return self.To.X + self.To.Width/2
        
    @property
    def Y(self):
        return self.To.Y + self.To.Height + 0.25

    @property
    def Outlet(self):
        return self.Inlet

    @property
    def InletVelocity(self):
        return self.Inlet.VolumeFlow / Circle(self.Diameter).Area

    @property
    def DrawContour(self):
        patch = plt.Circle(xy=(self.X,self.Y),radius=self.Diameter/2,lw=2,fill=True,edgecolor='black',facecolor='white')
        return patch

    @property
    def DrawLiquid(self):
        V = self.Outlet.VolumeFlow
        patch = self.DrawContour
        if V > 0.0:
            patch.set_facecolor(self.From.Outlet.Color)
        else:
            patch.set_facecolor('white')
        return patch

    @property
    def DrawStream(self):
        V = self.Outlet.VolumeFlow
        x = self.X
        y = self.Y
        dx = 10 * V
        if dx > 0.5:
            dx = 0.5
        dy = 10 * V
        if dy > 0.125:
            dy = 0.125
        points = [[x,y-0.05],[x-dx,y-dy*2],[x+dx,y-dy*2]]
        patch = plt.Polygon(points,fill=True,color=self.From.Outlet.Color)
        if self.Outlet.VolumeFlow == 0.0:
            patch.set_visible(False)
        return patch
    

class Splitter:
    def __init__(self, outlet_fractions):
        self.Name = None
        self.NumberOfOutlets = len(outlet_fractions)
        self.OutletFractions = outlet_fractions

    def OutletFlows(self, inlet):
        return [inlet.MassFlow * wi for wi in self.OutletFractions]

    def Outlets(self, inlet):
        o_list = [Stream(inlet.Components) for Fi in self.OutletFlows(inlet)]
        for i,o in enumerate(o_list):
            o.Components = inlet.Components
            o.isElectrolyte = inlet.isElectrolyte
            o.wi = inlet.MassFractions
            o.Mf = self.OutletFlows(inlet)[i]
            o.Temperature = inlet.Temperature
            o.Pressure = inlet.Pressure
        return o_list


class Mixer:
    def __init__(self):
        self.Name = None

    def InletFlows(self, inlets):
        return [i.MassFlow for i in inlets]

    def InletVolumeFlows(self, inlets):
        return [i.VolumeFlow for i in inlets]

    def NumberOfInlets(self, inlets):
        return len(inlets)

    def OutletFlow(self, inlets):
        return sum(self.InletFlows(inlets))

    def InletFractions(self, inlets):
        return [F/self.OutletFlow(inlets) for F in self.InletFlows(inlets)]

    def MassBalance(self, inlets):
        # WINDOWS VERSION #
        # dc_list = []
        # for i,I in enumerate(inlets):
        #     items = [(c,I.MassFlows[j]) for j,c in enumerate(I.Components)]
        #     od = OrderedDict(items)
        #     dc_list.append(od)
        #     # dc_list.append(dict(zip(I.Components, I.MassFlows)))
        # dd = OrderedDefaultDict(float)
        # # dd = OrderedDict
        # for dc in dc_list:
        #     for k,v in dc.items():
        #         dd[k] += v
        # # d = OrderedDict(dd)
        # return dd
        dc_list = []
        for i,I in enumerate(inlets):
            items = [(c,I.MassFlows[j]) for j,c in enumerate(I.Components)]
            od = OrderedDict(items)
            dc_list.append(od)
            # dc_list.append(dict(zip(I.Components, I.MassFlows)))
        dd = defaultdict(float)
        for dc in dc_list:
            for k,v in dc.items():
                dd[k] += v
        # d = OrderedDict(dd)
        return dd

    def OutletComponents(self, inlets):
        c_names = self.MassBalance(inlets).keys()
        return c_names

    def OutletFractions(self, inlets):
        return [mi/self.OutletFlow(inlets) for mi in self.MassBalance(inlets).values()]

    def OutletTemperature(self, inlets):
        num = 0
        den = 0
        for I in inlets:
            num += I.MassFlow * I.SpecificHeat * I.Temperature
            den += I.MassFlow * I.SpecificHeat
        T = num / den
        return T

    def Outlet(self, inlets):
        O = Stream(self.OutletComponents(inlets))
        O.isElectrolyte = self._isElectrolyte(inlets)
        O.wi = self.OutletFractions(inlets)
        O.Mf = self.OutletFlow(inlets)
        O.Temperature = self.OutletTemperature(inlets)
        O.isElectrolyte = self._isElectrolyte(inlets)
        # O.Pressure = ...
        return O

    def _isElectrolyte(self, inlets):
        elec = []
        for I in inlets:
            if I.isElectrolyte is True:
                elec.append(I.isElectrolyte)
        return any(elec)
