import numpy as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt
from dataclasses import dataclass, field

def ahrreniusRateConstant(Temperature : float = 273.15)
    pass
@dataclass

class Specie:
    counter : int = field(init= False, default= 0, repr= 0)
    name : str
    molarMass : float
    heatCapacityModel : str = "const" #""polynomial" also available"
    enthalpyFormation : float = 0.0
    entropyFormation : float = 0.0
    massSource : list[str] = field(default_factory= lambda : [])
    heatCapacityValue : float = 900
    heatCapacityCoefficients : list[float] = field(default_factory= lambda : [])
    
    def __post_init__(self):
        type(self).counter += 1
        self.ID = Specie.counter



@dataclass    
class Reaction:
    counter : int = field(init=False, default=0, repr=False)
    name : str
    stochiometricCoefficients : list[float] = field(default_factory= lambda : [])
    speciesID : list[int] = field(default_factory= lambda : [])
    speciesExponent : list[float] = field(default_factory= lambda : [])
    isReversible : bool = True
    ahrreniusPreExponent : float = 1.0
    ahrreniusActivationEnergy : float = 0.0
    species_registry: dict[int, Specie] = field(default_factory=dict, repr=False)

    def __post_init__(self):
        type(self).counter += 1
        self.ID = type(self).counter
        self.entropyChange = 0.0
        self.enthalpyChange = 0.0
        for nu, sid in zip(self.stochiometricCoefficients, self.speciesID):
            sp = self.species_registry[sid]
            self.enthalpyChange += nu * sp.enthalpyFormation
            self.entropyChange += nu * sp.entropyFormation
        

@dataclass
class Mixture:
    pass


@dataclass
class domainSetup:
    lenght : float
    diameter : float

@dataclass
class Zone:
    start : float #m
    end : float #m
    type : str = "null" #"cooling", "reaction", "null"
    heatSource : float = 0 #W/m3
    massSource : list[str] = field(default_factory= lambda : [])

@dataclass
class SolverSetup:
    divisionNumber : float = 120

@dataclass
class BoundaryConditions:
    massFlowrate : float = 50.0
    speciesMassFractions : list[float] = field(default_factory= lambda : [])
