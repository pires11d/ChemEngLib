
class Reaction:
    def __init__(self, components):
        self.Components = components
        self.ComponentNames = [c.Name for c in self.Components]
        self.Stoichiometry = [1 for c in self.Components]
        self.Temperature = 298.15
        self.Pressure = 101325.0
        self.NumberOfReactants = 1
        self.NumberOfProducts = 1
        self.isElementary = True
        self.Tank = None
        # Constants #
        self.KineticConstant = None
        self.ActivationEnergy = None
        self.EnthalpyOfReaction = None

    @property
    def Reactants(self):
        nR = self.NumberOfReactants
        R_list = []
        for i,r in enumerate(range(0,nR)):
            R_list.append(self.Components[i])
        return R_list

    @property
    def Products(self):
        nR = self.NumberOfReactants
        nP = self.NumberOfProducts
        n = len(self.Components)
        P_list = []
        for i,p in enumerate(range(0,nP)):
            P_list.append(self.Components[i+nR])
        return P_list

    @property
    def ReactantNames(self):
        return [c.Name for c in self.Reactants]

    @property
    def ProductNames(self):
        return [c.Name for c in self.Products]

    @property
    def Order(self):
        if self.isElementary == True:
            return self.NumberOfReactants
        else:
            return 1

    def ConsumptionRate(self, reactant):
        for i,r in enumerate(self.ReactantNames):
            if r == reactant:
                index = i
        # V = self.Tank.Mixture.Volume
        s = self.Stoichiometry[index]
        wi = self.Tank.Mixture.MolarFractions[index]
        C = self.Tank.Mixture.Moles * wi
        k = self.KineticConstant
        n = s*k*C
        return n

    def ProductionRate(self, product):
        for i,p in enumerate(self.ProductNames):
            if p == product:
                index = i
        # V = self.Tank.Mixture.Volume
        s = self.Stoichiometry[index]
        wi = self.Tank.Mixture.MolarFractions[index]
        C = self.Tank.Mixture.Moles * wi
        k = self.KineticConstant
        n = s*k*C
        return n