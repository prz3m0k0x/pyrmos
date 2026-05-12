import numpy as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from dataclasses import dataclass, field
import scipy.linalg
from typing import ClassVar, Dict, List


PI = 3.141592653589793
UNIVERSALGASCONSTANT = 8.31446261815324
REFERENCE_TEMPERATURE = 273.15
PRESSURE = 101325.0 #Pa

@dataclass
class Specie:
    counter: ClassVar[int] = 0

    name : str
    molarMass : float
    heatCapacityModel : str = "const" #""polynomial" also available"
    enthalpyFormation : float = 0.0
    entropyFormation : float = 0.0
    massSource : list[str] = field(default_factory= lambda : [])
    heatCapacityValue : float = 900
    heatCapacityCoefficients : list[float] = field(default_factory= lambda : [])
    
    def __post_init__(self):
        self.id = Specie.counter
        type(self).counter += 1
        if self.heatCapacityModel == "polynomial":
            self.heatCapacityValue = np.poly1d(self.heatCapacityCoefficients) #np polynomial object, function of some value (in this case T)

@dataclass
class Reaction:
    counter: ClassVar[int] = 0

    name: str
    stochiometricCoefficients: np.ndarray = field(default_factory=list)
    speciesID: np.ndarray = field(default_factory=list)
    speciesExponent: np.ndarray = field(default_factory=list)
    isReversible: bool = True
    ahrreniusPreExponent: float = 1.0
    ahrreniusActivationEnergy: float = 0.0  # J/mol
    species_registry: Dict[int, "Specie"] = field(default_factory=dict, repr=False)

    entropyChange: float = field(init=False)   # J/mol/K
    enthalpyChange: float = field(init=False)  # J/mol

    def __post_init__(self):
        type(self).counter += 1
        self.id = type(self).counter

        # consistency check
        if not (len(self.stochiometricCoefficients) == len(self.speciesID)):
            raise ValueError("stochiometricCoefficients and speciesID must have same length")

        self.enthalpyChange = 0.0
        self.entropyChange = 0.0

        for nu, sid in zip(self.stochiometricCoefficients, self.speciesID):
            sp = self.species_registry[sid]
            self.enthalpyChange += nu * sp.enthalpyFormation   # J/mol
            self.entropyChange += nu * sp.entropyFormation     # J/mol/K

    def forward_rate_constant(self, T: np.ndarray) -> np.ndarray:
        """Arrhenius forward rate constant k_f(T)."""
        return self.ahrreniusPreExponent * np.exp(-self.ahrreniusActivationEnergy / (UNIVERSALGASCONSTANT * T))

    def equilibrium_constant(self, T: np.ndarray) -> np.ndarray:
        """Equilibrium constant from ΔG = ΔH − TΔS."""
        delta_G = self.enthalpyChange - T * self.entropyChange  # J/mol
        return np.exp(-delta_G / (UNIVERSALGASCONSTANT * T))

    def backward_rate_constant(self, T: np.ndarray) -> np.ndarray:
        """Backward rate constant from k_b = k_f / K_eq."""
        if not self.isReversible:
            return 0.0
        k_f = self.forward_rate_constant(T)
        K_eq = self.equilibrium_constant(T)
        return k_f / K_eq
    

@dataclass
class Mixture:
    densityModel : str = "const" #can also be ideal-incompressible-gas
    densityValue : float = 1.2225 #kg/m3
    species_registry: Dict[int, "Specie"] = field(default_factory=dict, repr=False)

    def equivalentMolarMass(self):
        pass
        
    def idealGasDensity(self):
        pass

    def __post_init__(self):
        
        # if self.densityModel == "ideal-incompressible-gas":



@dataclass
class domainSetup:
    diameter: float
    massFlowRate : float
    inletMassFractions : np.ndarry = field(default_factory=list) 



@dataclass 
class Inlet:
    position: float = 0.0
    massFlowrate: float = 50.0
    speciesMassFractions: List[float] = field(default_factory=list)


@dataclass
class BoundaryConditions:
    counter: int = field(init=False, default=0, repr=False)

    def addBoundary(self):
        self.id = type(self).counter
        type(self).counter += 1


@dataclass
class Zone:
    counter: int = field(init=False, default=0, repr=False)
    lenght: float = 0.0
    type: str = "null"   # "cooling", "reaction", "null"
    heatSource: float = 0.0  # W/m3
    massSource: List[str] = field(default_factory=list)
    reaction : list[int] = field(default_factory=list)

    def __post_init__(self):
        type(self).counter += 1
        self.id = type(self).counter


@dataclass 
class Mesh:
    sizing: float = 0.005  # m
    bias: str = "soft"
    meshZones: List[Zone] = field(default_factory=list)

    # numpy arrays created in __post_init__
    cell_centers: np.ndarray = field(init=False, repr=False)
    cell_sizes:   np.ndarray = field(init=False, repr=False)
    cell_zone_id: np.ndarray = field(init=False, repr=False)
    cell_zone_type: np.ndarray = field(init=False, repr=False)
    cellZoneReaction : np.ndarray = field(init=False, repr=False) 
    lenght: float = field(init=False, default=0.0)

    def __post_init__(self, domain : domainSetup):

        zones = [z for z in self.meshZones if z.lenght > 0.0]
        if not zones:
            # empty mesh
            self.cell_centers = np.array([], dtype=float)
            self.cell_sizes   = np.array([], dtype=float)
            self.cell_volumes = np.array([], dtype=float)
            self.cell_zone_id = np.array([], dtype=int)
            self.cell_zone_type = np.array([], dtype=object)
            self.cellZoneReaction = np.array([], dtype=int)
            self.heatSource = np.array([], dtype=int)
            self.lenght = 0.0
            return

        Lz = np.array([z.lenght for z in zones], dtype=float)   # (Nz,)
        zone_ids = np.array([z.id for z in zones], dtype=int)   # (Nz,)
        zone_types = np.array([z.type for z in zones], dtype=object)


        n_cells = np.maximum(1, np.rint(Lz / self.sizing).astype(int))  # (Nz,)
        dz = Lz / n_cells                                               # (Nz,)

        n_cells_cumsum = np.cumsum(n_cells)
        n_total = int(n_cells_cumsum[-1])

        self.cell_centers = np.empty(n_total, dtype=float)
        self.cell_sizes   = np.empty(n_total, dtype=float)
        self.cell_volumes = np.empty(n_total, dtype=float)
        self.cell_zone_id = np.empty(n_total, dtype=int)
        self.cell_zone_type = np.empty(n_total, dtype=object)
        self.cellZoneReaction = np.empty(n_total, dtype=int)
        self.heatSource = np.empty(n_total, dtype=int)

        z_starts = np.concatenate(([0.0], np.cumsum(Lz[:-1])))

        start = 0
        for i, z in enumerate(zones):
            nc = n_cells[i]
            dz_i = dz[i]
            z0 = z_starts[i]


            idx = slice(start, start + nc)

            # cell centers: z_c = z0 + (k + 0.5) * dz_i, k = 0..nc-1
            k = np.arange(nc, dtype=float)
            self.cell_centers[idx] = z0 + (k + 0.5) * dz_i
            self.cell_sizes[idx]   = dz_i
            self.cell_volumes = self.cell_sizes * (domain.diameter** PI) / 4
            self.cell_zone_id[idx] = z.id
            self.cell_zone_type[idx] = z.type
            self.cellZoneReaction[idx] = np.array[z.reaction]
            self.heatSource[idx] = z.heatSource

            start += nc

        self.lenght = float(z_starts[-1] + Lz[-1])  # total length

zones = [
    Zone(lenght=0.10, type="reaction"),
    Zone(lenght=0.05, type="cooling"),
    Zone(lenght=0.10, type="reaction"),
]

mesh = Mesh(sizing=0.005, meshZones=zones)

print(len(mesh.cell_centers))  # total number of cells
print(mesh.cell_zone_type[:10])  # zone types per cell

class scalarField:
    def __init__(self, variable: str, Mesh: Mesh, type: str = "specie"):
        self.variable = variable
        self.type = type
        self.cellField = np.zeros((1, Mesh.n_cells))
        self.volumetricSources = np.zeros((1, Mesh.n_cells))
    
class globalMatrix:
    def __init__(self, Mixture: Mixture, Mesh: Mesh):
        self.area = (domainSetup.diameter**2 * PI / 4)
        self.density = Mixture.densityValue
        self.velocity = domainSetup.massFlowRate / (self.density * self.area)
        self.massFlux = domainSetup.massFlowRate  # F = m_dot
        self.cellVolume = Mesh.cell_volumes
        self.n_cells = Mesh.n_cells
    
def assemble_sources(Mesh: Mesh,
                     mixture: Mixture,
                     species_list: list[Specie],
                     reactions: dict[int, Reaction],
                     T: np.ndarray,
                     Y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Mesh: mesh object with cell_zone_type, cellZoneReaction
    mixture: gives density, etc.
    species_list: indexed by specie.id
    reactions: dict {reaction_id: Reaction}
    T: (n_cells,) temperature field
    Y: (n_species, n_cells) mass fractions

    returns:
        S_mass: (n_species, n_cells)  # kg/(m^3 s)
        q_rxn: (n_cells,)             # J/(m^3 s)
    """
    rho = mixture.densityValue
    n_species, n_cells = Y.shape
    S_mass = np.zeros((n_species, n_cells))
    q_rxn = np.zeros(n_cells)

    M = np.array([sp.molarMass for sp in species_list])

    reaction_cells = np.argwhere(Mesh.cell_zone_type == "reaction")
    
    for j in range(n_cells):
        if Mesh.cell_zone_type[j] != "reaction":
            continue

        r_id = Mesh.cellZoneReaction[j]
        rxn = reactions[r_id]

        # Build concentrations C_i at this cell
        C_local = np.zeros(len(rxn.speciesID))
        for idx, s_id in enumerate(rxn.speciesID):
            s = s_id  # assuming s_id matches index in species_list/Y
            Y_sj = Y[s, j]
            C_local[idx] = rho * Y_sj / M[s]  # mol/m^3

        Tj = T[j]

        # rate law (for now: irreversible, mass-action)
        k_f = rxn.forward_rate_constant(Tj)
        # product over reactants/products depending on rxn.speciesExponent
        rate_factor = np.prod(C_local ** rxn.speciesExponent)
        r_j = k_f * rate_factor  # mol/(m^3 s)

        # species sources
        for idx, s_id in enumerate(rxn.speciesID):
            s = s_id
            nu_i = rxn.stochiometricCoefficients[idx]  # mol stoich
            S_mass[s, j] += nu_i * r_j * M[s]  # kg/(m^3 s)

        # reaction heat source
        q_rxn[j] += -rxn.enthalpyChange * r_j  # J/(m^3 s)

    return S_mass, q_rxn

def solveEquation(A : np.ndarray, b : np.ndarray):
    return scipy.linalg.solve(a= A, b= b, assume_a= 'triadiagonal')