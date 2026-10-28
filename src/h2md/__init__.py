"""H2-only MD educational simulation and trajectory analysis package."""
from .analysis import (fit_self_diffusion, mean_squared_displacement,
                       number_density, radial_distribution)
from .io import load_lammps_dump, load_trajectory, save_trajectory
from .simulation import simulate_h2_centers
from .trajectory import Trajectory

__version__ = "0.1.0"
__all__ = ["Trajectory", "simulate_h2_centers", "mean_squared_displacement",
           "fit_self_diffusion", "radial_distribution", "number_density",
           "load_trajectory", "save_trajectory", "load_lammps_dump"]
