
class Reaction:
    def __init__(self, components):
        self.Components = components
        self.ComponentNames = [c.Name for c in self.Components]
        self.Temperature = 298.15
        self.Pressure = 101325.0