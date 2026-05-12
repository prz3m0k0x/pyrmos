import numpy as np
import scipy.optimize as optimize
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from dataclasses import dataclass, field
import scipy.linalg
from typing import ClassVar, Dict, List


PI = 3.141592653589793
UNIVERSALGASCONSTANT = 8.31446261815324
T_REF = 273.15
P_REF = 101325.0 #Pa

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
    def enthalpy(self, T):
        T = np.asarray(T, dtype=float)
        if self.heatCapacityModel == "const":

            cp_const = self.heatCapacityValue
            return cp_const * (T - T_REF)

        else:  # polynomial cp(T)
            cp_poly = self.heatCapacityValue  # np.poly1d
            H_poly = cp_poly.integ()          # antiderivative polynomial
            return H_poly(T) - H_poly(T_REF)

@dataclass
class Reaction:
    counter: ClassVar[int] = 0

    name: str
    stochiometricCoefficients: np.ndarray = field(default_factory=list)
    speciesID: np.ndarray = field(default_factory=list)
    speciesExponent: np.ndarray = field(default_factory=list)
    reversedSpecieExponent : np.ndarray = field(default_factory=list)
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
    
    def enthalpyReactionChange(self, T : np.ndarray) -> np.ndarray:
        T =np.asarray(T, dtype= float)
        Hr_T = np.zeros_like(T)

        for nu_i, sid in zip(self.stochiometricCoefficients, self.speciesID):
            sp = self.species_registry[sid]
            dH_i = sp.enthalpy(T)                           # H(T) - H(T_ref)
            H_i_T = sp.enthalpyFormation + dH_i                     # full H_i(T)
            Hr_T += nu_i * H_i_T

        return Hr_T
    
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
        pass
        # if self.densityModel == "ideal-incompressible-gas":



@dataclass
class domainSetup:
    diameter: float
    massFlowRate : float
    inletMassFractions : np.ndarray = field(default_factory=list) 


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

    def __post_init__(self):
        type(self).counter += 1
        self.id = type(self).counter


@dataclass 
class Mesh:
    domain : domainSetup
    sizing: float = 0.005  # m
    bias: str = "soft"
    meshZones: List[Zone] = field(default_factory=list)
    

    # numpy arrays created in __post_init__
    cell_centers: np.ndarray = field(init=False, repr=False)
    cell_sizes:   np.ndarray = field(init=False, repr=False)
    cell_zone_id: np.ndarray = field(init=False, repr=False)
    cell_zone_type: np.ndarray = field(init=False, repr=False)
    heatSource : np.ndarray = field(init=False, repr=False) 
    lenght: float = field(init=False, default=0.0)

    def __post_init__(self):

        zones = [z for z in self.meshZones if z.lenght > 0.0]
        if not zones:
            # empty mesh
            self.cell_centers = np.array([], dtype=float)
            self.cell_sizes   = np.array([], dtype=float)
            self.cell_volumes = np.array([], dtype=float)
            self.cell_zone_id = np.array([], dtype=int)
            self.cell_zone_type = np.array([], dtype=object)
            self.heatSource = np.array([], dtype=float)
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
        self.heatSource = np.empty(n_total, dtype=float)

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
            area = (self.domain.diameter**2) * PI / 4.0
            self.cell_volumes[idx] = self.cell_sizes[idx] * area
            self.cell_zone_id[idx] = z.id
            self.cell_zone_type[idx] = z.type
            self.heatSource[idx] = z.heatSource
            start += nc

        self.n_cells = n_total
        self.lenght = float(z_starts[-1] + Lz[-1])  # total length
        
class scalarField:
    def __init__(self, variable: str, Mesh: Mesh, type: str = "specie"):
        self.variable = variable
        self.type = type
        self.cellField = np.zeros((1, Mesh.n_cells))
        self.volumetricSources = np.zeros((1, Mesh.n_cells))
    
class globalMatrix:
    def __init__(self, domain: domainSetup, mixture: Mixture, Mesh: Mesh):
        self.area = (domain.diameter**2) * PI / 4.0
        self.density = mixture.densityValue
        self.velocity = domain.massFlowRate / (self.density * self.area)
        self.massFlux = domain.massFlowRate
        self.cellVolume = Mesh.cell_volumes
        self.n_cells = Mesh.n_cells
    
def assembleSources(Mesh: Mesh,
                    mixture: Mixture,
                    species_list: list[Specie],
                    reaction: Reaction,
                    T: np.ndarray,
                    Y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Mesh: mesh object with cell_zone_type, cellZoneReaction
    mixture: gives density, etc.
    species_list: list of all species, index == species.id
    reaction: single Reaction active in all 'reaction' cells
    T: (n_cells,) temperature field [K]
    Y: (n_species, n_cells) mass fractions

    returns:
        S_mass: (n_species, n_cells)  # kg/(m^3 s)
        S_q:    (n_cells,)            # J/(m^3 s)
    """

    rho = mixture.densityValue
    n_species, n_cells = Y.shape

    # outputs
    S_mass = np.zeros_like(Y, dtype=float)
    S_q = np.zeros(n_cells, dtype=float)
    # masks
    reaction_mask = (Mesh.cell_zone_type == "reaction")   # (n_cells,)
    cooling_mask = (Mesh.cell_zone_type == "cooling")   # for later

    S_q_cooling = np.zeros(n_cells, dtype=float)
    S_q_cooling[cooling_mask] = np.asarray(Mesh.heatSource)[cooling_mask]

    # if no reaction cells, just return zeros
    if not np.any(reaction_mask):
        return S_mass, S_q + S_q_cooling

    # species molar masses (assume index == specie.id)
    M = np.array([sp.molarMass for sp in species_list], dtype=float)  # (n_species,)

    # indices of species participating in this reaction
    species_idx = np.array(reaction.speciesID, dtype=int)            # (n_rxn_species,)
    M_rxn = M[species_idx]                                           # (n_rxn_species,)

    # restrict T and Y to reaction cells
    T_masked = T[reaction_mask]                                      # (n_reac_cells,)
    Y_rxn = Y[species_idx, :][:, reaction_mask]                      # (n_rxn_species, n_reac_cells)

    # concentrations C_i = rho * Y_i / M_i [mol/m^3]
    C_rxn = rho * Y_rxn / M_rxn[:, np.newaxis]                       # (n_rxn_species, n_reac_cells)

    # Arrhenius rate constants only in reaction cells
    k_f_masked = reaction.forward_rate_constant(T_masked)            # (n_reac_cells,)
    k_r_masked = reaction.backward_rate_constant(T_masked)           # (n_reac_cells,)

    # exponents aligned with species_idx
    alpha = np.asarray(reaction.speciesExponent, dtype=float)[:, np.newaxis]  # (n_rxn_species, 1)

    beta = np.asarray(reaction.reversedSpecieExponent, dtype=float)[:, np.newaxis]

    # forward and backward rates [mol/(m^3 s)]
    rate_forward_masked = k_f_masked * np.prod(C_rxn ** alpha, axis=0)  # (n_reac_cells,)
    # for now, assume irreversible; if reversible, uncomment and use:
    rate_backward_masked = k_r_masked * np.prod(C_rxn ** beta, axis=0)
    rate_masked = rate_forward_masked - rate_backward_masked


    rate = np.zeros(n_cells, dtype=float)
    rate[reaction_mask] = rate_masked

    # species mass sources: S_i = nu_i * rate * M_i  [kg/(m^3 s)]
    for loc, s_id in enumerate(species_idx):
        nu_i = reaction.stochiometricCoefficients[loc]
        S_mass[s_id, reaction_mask] += nu_i * rate_masked * M[s_id]

    Hr_T = reaction.enthalpyReactionChange(T)                        # (n_cells,), J/mol
    S_q[reaction_mask] = -Hr_T[reaction_mask] * rate_masked          # J/(m^3 s)

    return S_mass, S_q + S_q_cooling