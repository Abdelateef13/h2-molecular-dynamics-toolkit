"""Molecular-center structure and transport analysis."""

from dataclasses import dataclass

import numpy as np

from .trajectory import Trajectory


@dataclass(frozen=True)
class DiffusionFit:
    coefficient: float
    slope: float
    intercept: float
    r_squared: float
    start_time: float
    end_time: float


def mean_squared_displacement(trajectory: Trajectory, max_lag: int | None = None
                              ) -> tuple[np.ndarray, np.ndarray]:
    """Multi-time-origin ensemble MSD of UNWRAPPED molecular centers.

    Cost O(max_lag * frames * molecules). For compact educational trajectories.
    """
    x = trajectory.positions
    if max_lag is None:
        max_lag = max(2, (trajectory.n_frames - 1) // 3)
    if not (1 <= max_lag < trajectory.n_frames):
        raise ValueError("max_lag must be between 1 and n_frames - 1")
    msd = np.zeros(max_lag + 1)
    for lag in range(1, max_lag + 1):
        delta = x[lag:] - x[:-lag]
        msd[lag] = np.mean(np.einsum("...i,...i->...", delta, delta))
    return np.arange(max_lag + 1) * trajectory.frame_dt, msd


def fit_self_diffusion(times: np.ndarray, msd: np.ndarray,
                       fit_start: float = 0.40, fit_end: float = 0.85
                       ) -> DiffusionFit:
    """Fit Einstein MSD(t) = 6*D*t + b in a USER-CHOSEN 3D lag window.

    Fit boundaries are fractions of analyzed lag range. R² alone does not
    prove that a genuine long-time diffusive regime was reached.
    """
    times, msd = np.asarray(times), np.asarray(msd)
    if times.ndim != 1 or msd.shape != times.shape or len(times) < 5:
        raise ValueError("times/msd must be 1D matched arrays (>=5 samples)")
    if not np.all(np.isfinite(times)) or not np.all(np.isfinite(msd)):
        raise ValueError("times and msd must be finite")
    if not np.all(np.diff(times) > 0):
        raise ValueError("times must increase strictly")
    if not 0 <= fit_start < fit_end <= 1:
        raise ValueError("fit fractions must satisfy 0 <= start < end <= 1")
    last = len(times) - 1
    begin = max(1, int(np.ceil(fit_start * last)))
    end = min(len(times), int(np.floor(fit_end * last)) + 1)
    if end - begin < 3:
        raise ValueError("fit window must contain >=3 nonzero-lag points")
    x, y = times[begin:end], msd[begin:end]
    slope, intercept = np.polyfit(x, y, 1)
    fitted = slope * x + intercept
    sst = np.sum((y - y.mean()) ** 2)
    r2 = 1 - np.sum((y - fitted) ** 2) / sst if sst > 0 else float("nan")
    return DiffusionFit(float(slope / 6), float(slope), float(intercept),
                        float(r2), float(x[0]), float(x[-1]))


def radial_distribution(trajectory: Trajectory, bins: int = 70,
                        r_max: float | None = None, frame_stride: int = 5
                        ) -> tuple[np.ndarray, np.ndarray]:
    """g(r) for one species in a cubic periodic cell (minimum-image pairs).

    Counts each unique pair once; normalizes by ideal-gas finite-N number
    of unordered pairs and spherical shell volume.
    """
    if bins < 2 or frame_stride < 1:
        raise ValueError("bins >=2 and frame_stride >=1 required")
    length = trajectory.box_length
    if r_max is None:
        r_max = length / 2.0
    if not (0 < r_max <= length / 2.0):
        raise ValueError("r_max must be in (0, box_length/2]")
    edges = np.linspace(0.0, r_max, bins + 1)
    counts = np.zeros(bins)
    i, j = np.triu_indices(trajectory.n_molecules, 1)
    sampled = trajectory.wrapped_positions[::frame_stride]
    for positions in sampled:
        dr = positions[i] - positions[j]
        dr -= length * np.rint(dr / length)
        r = np.linalg.norm(dr, axis=1)
        counts += np.histogram(r, bins=edges)[0]
    shell_volumes = (4 * np.pi / 3) * (edges[1:] ** 3 - edges[:-1] ** 3)
    pairs = trajectory.n_molecules * (trajectory.n_molecules - 1) / 2
    expected = len(sampled) * pairs * shell_volumes / length**3
    return (edges[:-1] + edges[1:]) / 2.0, counts / expected


def number_density(trajectory: Trajectory) -> float:
    """Number of H2 centers per cubic length unit."""
    return trajectory.n_molecules / trajectory.box_length**3


def diffusion_in_si(fit: DiffusionFit, trajectory: Trajectory) -> float | None:
    """Convert Å²/fs to m²/s; no SI conversion for unspecified LJ units."""
    if trajectory.length_unit == "angstrom" and trajectory.time_unit == "fs":
        return fit.coefficient * 1e-5
    return None


def mass_density_g_cm3(trajectory: Trajectory) -> float | None:
    """Density of H2 molecular centers, only for real Å coordinate units."""
    if trajectory.length_unit != "angstrom":
        return None
    molar_mass = 2.01588  # g/mol, H2
    avogadro = 6.02214076e23  # molecules/mol
    return number_density(trajectory) * molar_mass / avogadro * 1e24
