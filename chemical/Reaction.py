
class Reaction:
    def __init__(self, mixture, order=[1]):
        self.Mixture = mixture
        self.Components = mixture.Components
        self.ComponentNames = mixture.ComponentNames
        self.Stoichiometry = [1 for c in self.Components]
        self.Temperature = 298.15
        self.Pressure = 101325.0
        self.NumberOfReactants = 1
        self.NumberOfProducts = 1
        self.Order = order
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

    def ConsumptionRate(self, reactant):
        index = None 
        for i,r in enumerate(self.Reactants):
            if r == reactant:
                index = i
        V = self.Mixture.Volume
        s = self.Stoichiometry[index]
        C = 1
        for i,o in enumerate(self.Order):
            C = C * self.Mixture.MolarConcentrations[i] ** o
        k = self.KineticConstant
        rate = s*k*C
        return rate * V

    def ProductionRate(self, product):
        index = None 
        nR = self.NumberOfReactants
        for i,p in enumerate(self.Products):
            if p == product:
                index = i
        V = self.Mixture.Volume
        s = self.Stoichiometry[index+nR]
        C = 1
        for i,o in enumerate(self.Order):
            C = C * self.Mixture.MolarConcentrations[i] ** o
        k = self.KineticConstant
        rate = s*k*C
        return rate * V

    