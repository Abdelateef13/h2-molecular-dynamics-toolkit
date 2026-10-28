"""Portable NPZ storage and small LAMMPS orthorhombic dump reader."""

from pathlib import Path

import numpy as np

from .trajectory import Trajectory


def save_trajectory(trajectory: Trajectory, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays = dict(positions=trajectory.positions, box_length=trajectory.box_length,
                  frame_dt=trajectory.frame_dt, length_unit=trajectory.length_unit,
                  time_unit=trajectory.time_unit)
    for key in ("potential_energy", "kinetic_energy", "temperature"):
        item = getattr(trajectory, key)
        if item is not None:
            arrays[key] = item
    np.savez_compressed(path, **arrays)


def load_trajectory(path: str | Path) -> Trajectory:
    with np.load(path, allow_pickle=False) as data:
        kwargs = {k: data[k] for k in data.files}
    kwargs["box_length"] = float(kwargs["box_length"])
    kwargs["frame_dt"] = float(kwargs["frame_dt"])
    kwargs["length_unit"] = str(kwargs["length_unit"])
    kwargs["time_unit"] = str(kwargs["time_unit"])
    return Trajectory(**kwargs)


def load_lammps_dump(path: str | Path, timestep_fs: float) -> Trajectory:
    """Read fixed-size cubic orthorhombic LAMMPS `dump custom` trajectory.

    Accepts xu/yu/zu or x/y/z plus ix/iy/iz; molecule centers must each be
    one LAMMPS particle. Requires units real (Å, fs) and equally spaced frames.
    Triclinic or changing-size boxes are intentionally unsupported in v0.1.
    """
    if timestep_fs <= 0:
        raise ValueError("timestep_fs must be positive")
    frames, timesteps = [], []
    origin_ref, edge_ref, ids_ref = None, None, None
    with open(path, encoding="utf-8") as f:
        while True:
            header = f.readline()
            if not header:
                break
            if header.strip() != "ITEM: TIMESTEP":
                raise ValueError("expected ITEM: TIMESTEP in LAMMPS dump")
            step = int(f.readline())
            if f.readline().strip() != "ITEM: NUMBER OF ATOMS":
                raise ValueError("expected ITEM: NUMBER OF ATOMS")
            n = int(f.readline())
            bound_header = f.readline().strip()
            if not bound_header.startswith("ITEM: BOX BOUNDS") or any(
                tag in bound_header.split() for tag in ("xy", "xz", "yz")
            ):
                raise ValueError("only orthorhombic boxes are supported")
            bounds = np.array([[float(a) for a in f.readline().split()]
                               for _ in range(3)])
            if bounds.shape != (3, 2):
                raise ValueError("invalid orthorhombic bounds")
            origin = bounds[:, 0]
            edges = bounds[:, 1] - bounds[:, 0]
            if not np.allclose(edges, edges[0], rtol=1e-6):
                raise ValueError("only cubic boxes supported in v0.1")
            if origin_ref is None:
                origin_ref, edge_ref = origin, float(edges[0])
            if not np.allclose(edges, edge_ref) or not np.allclose(origin, origin_ref):
                raise ValueError("changing box dimensions are unsupported")
            atoms_header = f.readline().strip().split()
            if atoms_header[:2] != ["ITEM:", "ATOMS"]:
                raise ValueError("expected ITEM: ATOMS")
            columns = atoms_header[2:]
            if "id" not in columns:
                raise ValueError("atom id is required to track molecules")
            cartesian = all(v in columns for v in ("xu", "yu", "zu"))
            wrapped = all(v in columns for v in ("x", "y", "z", "ix", "iy", "iz"))
            if not (cartesian or wrapped):
                raise ValueError("requires xu yu zu OR x y z ix iy iz")
            data = np.array([[float(v) for v in f.readline().split()]
                             for _ in range(n)])
            if data.shape != (n, len(columns)):
                raise ValueError("malformed atom line")
            order = np.argsort(data[:, columns.index("id")])
            ids = data[order, columns.index("id")].astype(int)
            if len(np.unique(ids)) != n:
                raise ValueError("duplicate molecule IDs")
            if ids_ref is None:
                ids_ref = ids
            elif not np.array_equal(ids_ref, ids):
                raise ValueError("atom IDs changed across frames")
            if cartesian:
                coordinates = data[order][:, [columns.index(k) for k in ("xu", "yu", "zu")]]
            else:
                xyz = data[order][:, [columns.index(k) for k in ("x", "y", "z")]]
                image = data[order][:, [columns.index(k) for k in ("ix", "iy", "iz")]]
                coordinates = xyz + image * edges
            frames.append(coordinates)
            timesteps.append(step)
    if len(frames) < 2:
        raise ValueError("need two or more dump frames")
    diffs = np.diff(timesteps)
    if np.any(diffs <= 0) or not np.all(diffs == diffs[0]):
        raise ValueError("dump frames must be equally spaced in ascending steps")
    return Trajectory(np.stack(frames), edge_ref, timestep_fs * float(diffs[0]),
                      length_unit="angstrom", time_unit="fs")
