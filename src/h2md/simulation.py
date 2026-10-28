"""Minimal velocity-Verlet, shifted-force Lennard-Jones MD demonstrator.

Each LJ site represents an H2 molecular *center*, not an individual H atom.
The reduced-unit parameters are illustrative, NOT a validated H2 force field.
"""

import numpy as np

from .trajectory import Trajectory


def _fcc_positions(cells: int, box_length: float) -> np.ndarray:
    basis = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    grid = np.indices((cells, cells, cells)).reshape(3, -1).T
    return ((grid[:, None, :] + basis[None, :, :]).reshape(-1, 3)
            * box_length / cells)


def _forces_and_potential(pos: np.ndarray, length: float, pairs: tuple,
                          cutoff: float) -> tuple[np.ndarray, float]:
    """Minimum-image LJ pairs, potential and radial force both zero at rc."""
    i, j = pairs
    dr = pos[i] - pos[j]
    dr -= length * np.rint(dr / length)
    r2 = np.einsum("ij,ij->i", dr, dr)
    chosen = (r2 < cutoff * cutoff) & (r2 > 1e-14)
    d = dr[chosen]
    r = np.sqrt(r2[chosen])
    inv = 1.0 / r
    inv6 = inv**6
    inv12 = inv6 * inv6
    fc = 24.0 * (2.0 / cutoff**13 - 1.0 / cutoff**7)
    uc = 4.0 * (1.0 / cutoff**12 - 1.0 / cutoff**6)
    radial_force = 24.0 * (2.0 * inv12 - inv6) * inv
    fij = ((radial_force - fc) / r)[:, None] * d
    f = np.zeros_like(pos)
    np.add.at(f, i[chosen], fij)
    np.add.at(f, j[chosen], -fij)
    potential = float(np.sum(4.0 * (inv12 - inv6) - uc + (r - cutoff) * fc))
    return f, potential


def simulate_h2_centers(
    cells: int = 3,
    number_density: float = 0.35,
    target_temperature: float = 2.0,
    equilibration_steps: int = 500,
    production_steps: int = 3000,
    dt: float = 0.004,
    sample_every: int = 10,
    cutoff: float = 2.5,
    seed: int = 42,
) -> Trajectory:
    """Run toy H2-center LJ MD in reduced units (epsilon=sigma=m=kB=1).

    Equilibrate by periodic velocity rescaling, then use NVE for production.
    Parameters are pedagogical; results MUST NOT be quoted as real H2 data.
    """
    if cells < 2 or number_density <= 0 or target_temperature <= 0:
        raise ValueError("cells >= 2 and positive density/temperature required")
    if equilibration_steps < 0 or production_steps < 1 or sample_every < 1 or dt <= 0:
        raise ValueError("invalid simulation steps or timestep")
    n = 4 * cells**3
    length = (n / number_density) ** (1.0 / 3.0)
    if cutoff <= 0 or cutoff >= length / 2.0:
        raise ValueError("cutoff must be positive and strictly less than half the box")
    pos = _fcc_positions(cells, length)
    rng = np.random.default_rng(seed)
    vel = rng.normal(size=(n, 3))
    vel -= vel.mean(axis=0)
    ndof = 3 * n - 3
    temp = (vel * vel).sum() / ndof
    vel *= np.sqrt(target_temperature / temp)
    pairs = np.triu_indices(n, k=1)
    force, pot = _forces_and_potential(pos % length, length, pairs, cutoff)

    def advance(pos: np.ndarray, vel: np.ndarray, force: np.ndarray):
        vel_half = vel + 0.5 * dt * force
        pos = pos + dt * vel_half  # continuous, UNWRAPPED coordinates
        new_force, new_pot = _forces_and_potential(pos % length, length, pairs, cutoff)
        vel = vel_half + 0.5 * dt * new_force
        return pos, vel, new_force, new_pot

    for t in range(equilibration_steps):
        pos, vel, force, pot = advance(pos, vel, force)
        if (t + 1) % 10 == 0:
            temp = (vel * vel).sum() / ndof
            vel *= np.sqrt(target_temperature / temp)

    coordinates, potentials, kinetics, temperatures = [], [], [], []

    def save_frame():
        kin = 0.5 * np.sum(vel * vel)
        coordinates.append(pos.copy())
        potentials.append(pot)
        kinetics.append(kin)
        temperatures.append(2.0 * kin / ndof)

    save_frame()
    for t in range(1, production_steps + 1):
        pos, vel, force, pot = advance(pos, vel, force)
        if t % sample_every == 0:
            save_frame()

    return Trajectory(
        positions=np.stack(coordinates), box_length=length,
        frame_dt=dt * sample_every, length_unit="sigma", time_unit="tau",
        potential_energy=np.array(potentials), kinetic_energy=np.array(kinetics),
        temperature=np.array(temperatures),
    )
