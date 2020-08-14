import math
import matplotlib.pyplot as plt
from matplotlib import animation
from collections import Counter, defaultdict, OrderedDict
from chemical import *
from utilities.Geometry import *
from .Piping import Mixer


class Tank:
    def __init__(self, max_volume=1.0):
        self.Name = None
        self.MaxVolume = max_volume
        self.Inlets = []
        self.OutletVolumeFlow = 0.0
        self.Mixture = None
        # Drawing #
        self.dt = 0.1
        self.X = 0 
        self.Y = 0
        self.Height = 1.0
        self.Width = 1.0

    @property
    def ConeBottom(self):
        return False

    @property
    def ConeVolume(self):
        return 0.0

    @property
    def ConeHeight(self):
        return 0.0

    @property
    def InletVolumeFlow(self):
        IN = 0
        if len(self.Inlets) > 0:
            for i in self.Inlets:
                IN += i.VolumeFlow
        return IN

    @property
    def NextVolume(self):
        IN = self.InletVolumeFlow
        OUT = self.OutletVolumeFlow
        V = (IN - OUT) * self.dt + self.Mixture.Volume
        if V < 0.0:
            V = 0.0
        self.Mixture.V = V
        return V

    @property
    def InletMolarFlow(self):
        IN = 0
        if len(self.Inlets) > 0:
            for i in self.Inlets:
                IN += i.MolarFlow
        return IN

    @property
    def OutletMolarFlow(self):
        return self.Outlet.MolarFlow

    @property
    def NextMoles(self):
        IN = self.InletMolarFlow
        OUT = self.OutletMolarFlow
        M = (IN - OUT) * self.dt + self.Mixture.Moles
        if M < 0.0:
            M = 0.0
        self.Mixture.M = M
        return M

    @property
    def NextMixture(self):
        m = Mixture(self.Mixture.Components)
        m.V = self.Mixture.Volume
        m.Temperature = self.Mixture.Temperature
        m.Pressure = self.Mixture.Pressure

        # INLET BALANCE #
        IN = Mixer().OutletFlow(self.Inlets) * self.dt
        if IN != 0:
            wi_in = Mixer().OutletFractions(self.Inlets)
            m.Temperature = Mixer().OutletTemperature(self.Inlets)
        else:
            wi_in = self.Mixture.MassFractions
        m.wi = wi_in
        wi_list = []

        if self.Mixture.Volume != 0:
            # MASS BALANCE #
            for i,wi in enumerate(wi_in):
                num = wi * IN + self.Mixture.Mass * self.Mixture.MassFractions[i]
                den = IN + self.Mixture.Mass
                wi = num / den
                wi_list.append(wi)
            m.wi = wi_list
            # ENERGY BALANCE #
            pseudoStream = Stream(self.Mixture.Components)
            pseudoStream.wi = self.Mixture.MassFractions
            pseudoStream.Mf = self.Mixture.Mass / self.dt
            pseudoStream.Temperature = self.Mixture.Temperature

            pseudoInlets = [s for s in self.Inlets]
            pseudoInlets.append(pseudoStream)
            m.Temperature = Mixer().OutletTemperature(pseudoInlets)

        self.Mixture = m
        return m

    @property
    def NextTime(self):
        for i in self.Inlets:
            if i.Vf == None:
                self.NextMoles
                break
        else:
            self.NextVolume
        self.NextMixture
        return None

    @property
    def Outlet(self):
        O = Stream(self.Mixture.Components)
        O.wi = self.Mixture.MassFractions
        if self.Mixture.Volume == 0.0:
            if self.InletVolumeFlow > 0.0:
                O.Vf = self.InletVolumeFlow
            else:
                O.Vf = 0.0
        else:
            O.Vf = self.OutletVolumeFlow
        O.Temperature = self.Mixture.Temperature
        O.Pressure = self.Mixture.Pressure
        return O

    @property
    def LiquidHeight(self):
        return self.Height * self.NextVolume / self.MaxVolume

    @property
    def LiquidWidth(self):
        return self.Width

    @property
    def DrawContour(self):
        X = self.X
        Y = self.Y
        Hc = self.ConeHeight
        H = self.Height
        W = self.Width*1.01
        points = [[X,Y+H],[X,Y+Hc],[X+W/2,Y],[X+W,Y+Hc],[X+W,Y+H]]
        patch = plt.Polygon(points, closed=None, fill=None, lw=2, edgecolor='black')
        return patch

    @property
    def DrawLiquid(self):
        zero = 0.01
        X = self.X
        Y = self.Y
        Hc = self.ConeHeight
        W = self.Width
        h = self.LiquidHeight
        w = (W-self.LiquidWidth)/2
        if self.LiquidWidth < self.Width:
            points = [[X+w+zero,Y+h+zero],[X+W/2,Y+zero],[X+W-w-zero,Y+h+zero]]
        else:
            points = [[X+zero*2,Y+h+zero],[X+zero*2,Y+Hc+zero],[X+W/2,Y+zero],[X+W-zero,Y+Hc+zero],[X+W-zero,Y+h+zero]]
        patch = plt.Polygon(points, closed=True, fill=True,
                            facecolor=self.Mixture.Color,
                            hatch=self.Mixture.Hatch)
        return patch

    @property
    def DrawTopArrow(self):
        x = self.X + self.Width/2
        y = self.Y + self.Height
        l = 0.1
        s = 0.05
        patch = plt.arrow(x,y,0,-l,lw=0.5,head_width=s,head_length=s,fill=True,color='black',length_includes_head=True)
        if self.InletVolumeFlow > 0.0:
            patch.set_visible(True)
        else:
            patch.set_visible(False)
        if self.Mixture.Volume > self.MaxVolume:
            patch.set_color('red')
        return patch

    @property
    def DrawBottomArrow(self):
        x = self.X + self.Width/2 
        y = self.Y
        l = 0.1
        s = 0.05
        patch = plt.arrow(x,y,0,-l,lw=0.5,head_width=s,head_length=s,fill=True,color='black',length_includes_head=True)
        if self.OutletVolumeFlow > 0.0:
            patch.set_visible(True)
        else:
            patch.set_visible(False)
        return patch


class Hopper(Tank):
    def __init__(self, initial_angle, final_angle, Hmin, Hmax, r1, r2, R):
        self.Name = None        
        self.Inlets = []
        self.OutletVolumeFlow = 0.0
        self.Mixture = None
        # Hopper attributes #
        self.InitialAngle = initial_angle
        self.FinalAngle = final_angle
        self.AngularSize = self.FinalAngle - self.InitialAngle
        self.MinAngle = 0
        self.MaxAngle = 360
        self.MinHeight = Hmin
        self.MaxHeight = Hmax
        self.InternalRadius = r1
        self.MiddleRadius = r2
        self.ExternalRadius = R
        # Drawing #
        self.dt = 0.1        
        self.X = 0 
        self.Y = 0
        self.Height = 1.0
        self.Width = self.AngularSize/self.ExternalRadius

    @property
    def ActualHeight(self):
        Hmin = self.MinHeight
        Hmax = self.MaxHeight
        thetaMin = math.radians(self.MinAngle)
        thetaMax = math.radians(self.MaxAngle)
        theta = math.radians(self.InitialAngle)
        return Hmax-(Hmax-Hmin)/(thetaMax-thetaMin)*(thetaMax-theta)

    @property
    def MaxVolume(self):
        PI = math.pi
        r1 = self.InternalRadius
        r2 = self.MiddleRadius
        R = self.ExternalRadius
        Hmin = self.MinHeight
        Hmax = self.MaxHeight
        h_ = self.ActualHeight
        dTheta = math.radians(self.AngularSize)
        Vminor = PI * Hmin * ((3 * r2 * r1 ** 2 - r2 ** 3 - 2 * r1 ** 3)) / (3 * (r1 - r2))
        Vtotal = 0.5 * PI * (h_ - Hmin + 2 * Hmax) * (R ** 2 - r1 ** 2)
        V = Vtotal - Vminor
        V = V * (dTheta) / (2 * PI)
        return V
    

class recTank(Tank):
    def __init__(self, length, width, height):
        self.Name = None
        self.Inlets = []
        self.OutletVolumeFlow = 0.0
        self.Mixture = None
        # Rectangular tank attributes #
        self.Length = length
        # Drawing #
        self.dt = 0.1
        self.X = 0 
        self.Y = 0    
        self.Height = height
        self.Width = width

    @property
    def MaxVolume(self):
        return Cube(self.Height,self.Length,self.Width).Volume

    
class cylTank(Tank):
    def __init__(self, diameter, cylinder_height, cone_angle=0.0):
        self.Name = None
        self.Inlets = []
        self.OutletVolumeFlow = 0.0
        self.Mixture = None
        # Cylindrical tank attributes #
        self.CylinderHeight = cylinder_height
        self.Diameter = diameter
        self.ConeAngle = cone_angle   
        # Drawing #
        self.dt = 0.1
        self.X = 0 
        self.Y = 0  

    @property
    def ConeBottom(self):
        if self.ConeAngle > 0.0:
            return True
        else:
            return False

    @property
    def ConeHeight(self):
        return Cone(diameter=self.Diameter,angle=self.ConeAngle).Height

    @property
    def ConeVolume(self):
        return Cone(diameter=self.Diameter,height=self.ConeHeight).Volume

    @property
    def CylinderVolume(self):
        return Cylinder(diameter=self.Diameter,height=self.CylinderHeight).Volume

    @property
    def MaxVolume(self):
        Vcone = self.ConeVolume
        Vcyl = self.CylinderVolume
        return Vcone + Vcyl

    @property
    def LiquidHeight(self):
        if self.NextVolume <= self.ConeVolume:
            tan = math.tan(math.radians(self.ConeAngle))
            h = (3.0*self.NextVolume*(tan**2)/math.pi)**(1.0/3.0)
        else:
            V = self.NextVolume - self.ConeVolume
            dh = V / Cylinder(self.Diameter,self.CylinderHeight).BaseArea
            h = self.ConeHeight + dh
        return h

    @property
    def LiquidWidth(self):
        if self.NextVolume <= self.ConeVolume:
            d = Cone(height=self.LiquidHeight,angle=self.ConeAngle).Diameter
        else:
            d = self.Diameter
        return d
    
    @property
    def Height(self):
        return self.CylinderHeight + self.ConeHeight

    @property
    def Width(self):
        return self.Diameter


class Reactor(Tank):
    def __init__(self, reaction, max_volume=1.0):
        self.Name = None
        self.Reaction = reaction
        self.MaxVolume = max_volume
        self.Inlets = []
        self.OutletVolumeFlow = 0.0
        self.Mixture = reaction.Mixture
        # Drawing #
        self.dt = 0.1
        self.X = 0 
        self.Y = 0
        self.Height = 1.0
        self.Width = 1.0

    @property
    def NextMoles(self):
        nR = len(self.Reaction.Reactants)
        N = self.Mixture.Moles
        zi = self.Mixture.MolarFractions
        ni = []
        dCr = []
        if N > 0:
            for i,r in enumerate(self.Reaction.Reactants):
                dCr.append(self.Reaction.ConsumptionRate(r) * self.dt)
                N -= dCr[i]
                nr = self.Reaction.Mixture.ComponentMoles[i] - dCr[i]
                ni.append(nr)
            dCp = []
            for i,p in enumerate(self.Reaction.Products):
                dCp.append(self.Reaction.ProductionRate(p) * self.dt)
                N += dCp[i]
                np = self.Mixture.ComponentMoles[i+nR] + dCp[i]
                ni.append(np)
            if N > 0:
                zi = []
                for i,n in enumerate(ni):
                    if n > 0:
                        zi.append(n/N)
                    else:
                        zi.append(0)
            else:
                N = 0
                zi = list([0 for i in self.Mixture.Components]) 
        self.Mixture.N = N
        self.Mixture.zi = zi  
        return self.Mixture.N
        
