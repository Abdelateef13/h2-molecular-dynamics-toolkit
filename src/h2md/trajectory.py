"""Trajectory container with explicit coordinate and time units."""
from dataclasses import dataclass

import numpy as np


@dataclass
class Trajectory:
    """Unwrapped molecular-center coordinates, shape (frames, molecules, xyz).

    `box_length` is the cubic periodic cell edge in the same length unit
    as `positions`; `frame_dt` is the physical time between saved frames.
    """

    positions: np.ndarray
    box_length: float
    frame_dt: float
    length_unit: str = "sigma"
    time_unit: str = "tau"
    potential_energy: np.ndarray | None = None
    kinetic_energy: np.ndarray | None = None
    temperature: np.ndarray | None = None

    def __post_init__(self) -> None:
        self.positions = np.asarray(self.positions, dtype=np.float64)
        if self.positions.ndim != 3 or self.positions.shape[2] != 3:
            raise ValueError("positions must have shape (n_frames, n_molecules, 3)")
        if self.positions.shape[0] < 2 or self.positions.shape[1] < 2:
            raise ValueError("at least two frames and two molecules are needed")
        if not np.all(np.isfinite(self.positions)):
            raise ValueError("positions must be finite")
        if not np.isfinite(self.box_length) or self.box_length <= 0:
            raise ValueError("box_length must be positive")
        if not np.isfinite(self.frame_dt) or self.frame_dt <= 0:
            raise ValueError("frame_dt must be positive")
        for name in ("potential_energy", "kinetic_energy", "temperature"):
            val = getattr(self, name)
            if val is not None:
                val = np.asarray(val, dtype=np.float64)
                if val.shape != (self.n_frames,):
                    raise ValueError(f"{name} must have one value per frame")
                setattr(self, name, val)

    @property
    def n_frames(self) -> int:
        return self.positions.shape[0]

    @property
    def n_molecules(self) -> int:
        return self.positions.shape[1]

    @property
    def times(self) -> np.ndarray:
        return np.arange(self.n_frames, dtype=float) * self.frame_dt

    @property
    def wrapped_positions(self) -> np.ndarray:
        return np.mod(self.positions, self.box_length)
