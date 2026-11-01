"""H2-only MD educational simulation and trajectory analysis package."""
from .analysis import (fit_self_diffusion, mean_squared_displacement,
                       number_density, radial_distribution)
from .io import load_lammps_dump, load_trajectory, save_trajectory
from .simulation import simulate_h2_centers
from .trajectory import Trajectory
from .atomistic import (AtomicH2Trajectory, simulate_atomistic_h2,
                        load_atomistic_lammps_dump, save_atomistic_trajectory,
                        load_atomistic_trajectory)

__version__ = "0.2.0"
__all__ = ["Trajectory", "simulate_h2_centers", "mean_squared_displacement",
           "fit_self_diffusion", "radial_distribution", "number_density",
           "load_trajectory", "save_trajectory", "load_lammps_dump",
           "AtomicH2Trajectory", "simulate_atomistic_h2",
           "load_atomistic_lammps_dump", "save_atomistic_trajectory",
           "load_atomistic_trajectory"]
