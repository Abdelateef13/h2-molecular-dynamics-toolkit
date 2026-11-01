"""Two-atom H2 molecular dynamics demonstrator (dimensionless teaching model).

Each molecule consists of *two independently integrated H atoms* connected by
an explicit harmonic bond. Intermolecular H--H LJ terms act between different
molecules; intramolecular LJ terms are EXCLUDED.

The parameters are illustrative dimensionless numbers, NOT a validated H2
potential. Classical bonds omit quantum vibration, and classical rotation is
not quantitatively predictive for room-temperature hydrogen.
"""
from dataclasses import dataclass
import numpy as np

from .simulation import _fcc_positions
from .trajectory import Trajectory


@dataclass
class AtomicH2Trajectory:
    """Unwrapped H atom positions, (frames, molecules, 2 atoms, xyz)."""
    atoms: np.ndarray
    box_length: float
    frame_dt: float
    length_unit: str = "sigma"
    time_unit: str = "tau"
    potential_energy: np.ndarray | None = None
    kinetic_energy: np.ndarray | None = None
    temperature: np.ndarray | None = None

    def __post_init__(self):
        self.atoms = np.asarray(self.atoms, dtype=float)
        if self.atoms.ndim != 4 or self.atoms.shape[2:] != (2, 3) or min(self.atoms.shape[:2]) < 2:
            raise ValueError("atoms must have shape (frames>=2, molecules>=2, 2, 3)")
        if not np.isfinite(self.atoms).all() or self.box_length <= 0 or self.frame_dt <= 0:
            raise ValueError("invalid positions/box/timestep")

    @property
    def n_frames(self):
        return self.atoms.shape[0]

    @property
    def n_molecules(self):
        return self.atoms.shape[1]

    @property
    def com(self):
        """Equal-H-mass, *unwrapped* molecular centers of mass."""
        return self.atoms.mean(axis=2)

    def center_trajectory(self):
        return Trajectory(self.com, self.box_length, self.frame_dt,
                          length_unit=self.length_unit, time_unit=self.time_unit,
                          potential_energy=self.potential_energy,
                          kinetic_energy=self.kinetic_energy, temperature=self.temperature)

    def bond_lengths(self):
        """Unwrapped intramolecular H--H separation per frame and molecule."""
        return np.linalg.norm(self.atoms[:, :, 1] - self.atoms[:, :, 0], axis=-1)

    def orientations(self):
        d = self.atoms[:, :, 1] - self.atoms[:, :, 0]
        return d / np.linalg.norm(d, axis=-1)[..., None]


def _atomistic_forces(x: np.ndarray, length: float, i: np.ndarray, j: np.ndarray,
                      cutoff: float, bond_k: float, r0: float):
    """LJ shifted-force intermolecular + harmonic bond. x is (n_mol,2,3)."""
    a = x.reshape(-1, 3)
    dr = a[i] - a[j]
    dr -= length * np.rint(dr / length)
    r2 = np.einsum('ij,ij->i', dr, dr)
    use = (r2 < cutoff**2) & (r2 > 1e-14)
    r = np.sqrt(r2[use]); d = dr[use]
    inv = 1 / r
    inv6 = inv**6; inv12 = inv6**2
    # F_c is scalar radial force: used to shift force and potential smoothly.
    fc = 24 * (2 / cutoff**13 - 1 / cutoff**7)
    uc = 4 * (1 / cutoff**12 - 1 / cutoff**6)
    radial = 24 * (2 * inv12 - inv6) * inv
    fij = ((radial - fc) / r)[:, None] * d
    f = np.zeros_like(a)
    np.add.at(f, i[use], fij)
    np.add.at(f, j[use], -fij)
    potential = float(np.sum(4 * (inv12 - inv6) - uc + (r-cutoff) * fc))

    # Bonded pairs do not interact through the LJ potential.
    d_bond = x[:, 0] - x[:, 1]
    r_bond = np.linalg.norm(d_bond, axis=1)
    if np.any(r_bond < 1e-12):
        raise ValueError("zero H-H bond length")
    bond_force = (-bond_k * (r_bond - r0) / r_bond)[:, None] * d_bond
    f = f.reshape(x.shape)
    f[:, 0] += bond_force
    f[:, 1] -= bond_force
    potential += float(0.5 * bond_k * np.sum((r_bond-r0)**2))
    return f, potential


def simulate_atomistic_h2(cells: int = 2, number_density: float = .24,
                           target_temperature: float = .8,
                           equilibration_steps: int = 400,
                           production_steps: int = 2500,
                           dt: float = .0008, sample_every: int = 10,
                           bond_length: float = .3, bond_k: float = 1200.,
                           cutoff: float = 2.5, seed: int = 42) -> AtomicH2Trajectory:
    """NVE MD for flexible H2 dimers, periodic cubic box, reduced units.

    Two hydrogen atoms per H2. Each atomic mass = 0.5, so H2 mass = 1.
    LJ sigma = epsilon = k_B = 1. NO physically validated parameterization.
    """
    if cells < 2 or number_density <= 0 or target_temperature <= 0:
        raise ValueError('cells>=2 and density/temperature >0 are required')
    if min(dt, cutoff, bond_length, bond_k) <= 0 or sample_every < 1 or production_steps < 1 or equilibration_steps < 0:
        raise ValueError('invalid timesteps or potential parameters')
    n = 4 * cells**3
    length = (n / number_density)**(1/3)
    if cutoff >= length/2:
        raise ValueError('cutoff must be less than half cubic box edge')
    rng = np.random.default_rng(seed)
    center = _fcc_positions(cells, length)
    direction = rng.normal(size=(n,3))
    direction /= np.linalg.norm(direction, axis=1)[:, None]
    x = center[:,None,:] + (.5*bond_length) * np.stack((-direction,direction),axis=1)
    mass = .5
    v = rng.normal(size=x.shape) * np.sqrt(target_temperature/mass)
    v -= np.mean(v.reshape(-1,3), axis=0)
    dof = 6*n - 3  # flexible two-atom molecules, no constraints (remove COM drift)
    v *= np.sqrt(target_temperature / (mass * np.sum(v*v) / dof))
    ai, aj = np.triu_indices(2*n,1)
    intermolecular = (ai//2) != (aj//2)
    ai, aj = ai[intermolecular], aj[intermolecular]
    f, pot = _atomistic_forces(x, length, ai, aj, cutoff, bond_k, bond_length)

    def advance(x,v,f):
        half = v + .5*dt*f/mass
        newx = x + dt*half  # stay UNWRAPPED to preserve continuous atom/bond trajectories
        newf,pot = _atomistic_forces(newx, length, ai, aj, cutoff, bond_k, bond_length)
        newv = half + .5*dt*newf/mass
        return newx,newv,newf,pot

    for step in range(equilibration_steps):
        x,v,f,pot = advance(x,v,f)
        if (step+1) % 10 == 0:
            current_t = mass * np.sum(v*v)/dof
            v *= np.sqrt(target_temperature/current_t)

    atoms, pe, ke, temp = [], [], [], []
    def record():
        kin = .5*mass*np.sum(v*v)
        atoms.append(x.copy()); pe.append(pot); ke.append(kin); temp.append(2*kin/dof)
    record()
    for step in range(1,production_steps+1):
        x,v,f,pot = advance(x,v,f)
        if step % sample_every == 0:
            record()
    return AtomicH2Trajectory(np.stack(atoms), length, dt*sample_every,
                              potential_energy=np.array(pe), kinetic_energy=np.array(ke),
                              temperature=np.array(temp))


def save_atomistic_trajectory(tr: AtomicH2Trajectory, path):
    """Save atomic coordinates, not just molecular centers, in .npz format."""
    from pathlib import Path
    p = Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    arrays = dict(atoms=tr.atoms, box_length=tr.box_length, frame_dt=tr.frame_dt,
                  length_unit=tr.length_unit, time_unit=tr.time_unit)
    for key in ('potential_energy','kinetic_energy','temperature'):
        val = getattr(tr,key)
        if val is not None: arrays[key] = val
    np.savez_compressed(p,**arrays)


def load_atomistic_trajectory(path):
    with np.load(path, allow_pickle=False) as d:
        kwargs = {k:d[k] for k in d.files}
    kwargs.update(box_length=float(kwargs['box_length']), frame_dt=float(kwargs['frame_dt']),
                  length_unit=str(kwargs['length_unit']), time_unit=str(kwargs['time_unit']))
    return AtomicH2Trajectory(**kwargs)


def load_atomistic_lammps_dump(path, timestep_fs: float) -> AtomicH2Trajectory:
    """Load constant-cubic-cell, 2-H-atom H2 LAMMPS dump with molecular IDs.

    Required columns: id mol xu yu zu (or id mol x y z ix iy iz).
    Requires LAMMPS `units real`, x in Å, step time in fs.
    LAMMPS must output consistent unwrapped positions across frames.
    Rejects mixed species, changing dimensions, triclinic boxes, and
    atom lists whose molecule IDs don't define *exactly* two atoms each.
    """
    if timestep_fs <= 0: raise ValueError('timestep_fs > 0 required')
    snapshots=[]; steps=[]; box=None; ids0=None; mol0=None
    with open(path,encoding='utf-8') as f:
        while True:
            header=f.readline()
            if not header: break
            if header.strip() != 'ITEM: TIMESTEP': raise ValueError('expected ITEM: TIMESTEP')
            step=int(f.readline())
            if f.readline().strip() != 'ITEM: NUMBER OF ATOMS': raise ValueError('expected NUMBER OF ATOMS')
            n=int(f.readline())
            bh=f.readline().strip().split()
            if bh[:3]!=['ITEM:','BOX','BOUNDS'] or any(k in bh for k in ('xy','xz','yz')):
                raise ValueError('orthorhombic cubic boxes only')
            bounds=np.array([[float(x) for x in f.readline().split()] for _ in range(3)])
            if bounds.shape != (3,2): raise ValueError('invalid bounds')
            edges=bounds[:,1]-bounds[:,0]
            if not np.allclose(edges,edges[0]): raise ValueError('cubic boxes only')
            if box is None: box=float(edges[0])
            elif not np.allclose(edges,box): raise ValueError('changing cell not supported')
            fields=f.readline().strip().split()
            if fields[:2]!=['ITEM:','ATOMS']: raise ValueError('expected ITEM: ATOMS')
            col=fields[2:]
            if not all(c in col for c in ('id','mol')): raise ValueError('id and mol columns required')
            is_unwrapped=all(c in col for c in ('xu','yu','zu'))
            has_images=all(c in col for c in ('x','y','z','ix','iy','iz'))
            if not (is_unwrapped or has_images):
                raise ValueError('requires xu yu zu or x y z ix iy iz')
            rows=np.array([[float(v) for v in f.readline().split()] for _ in range(n)])
            if rows.shape != (n,len(col)) or n%2: raise ValueError('invalid atom data / odd atom count')
            order=np.argsort(rows[:,col.index('id')])
            rows=rows[order]
            ids=rows[:,col.index('id')].astype(int)
            mol=rows[:,col.index('mol')].astype(int)
            if np.unique(ids).size!=n: raise ValueError('duplicate atom ID')
            if ids0 is None:
                ids0=ids.copy();mol0=mol.copy()
                counts=np.unique(mol,return_counts=True)[1]
                if np.any(counts!=2) or np.any(mol<=0): raise ValueError('exactly 2 H atoms per molecule required')
            elif not np.array_equal(ids,ids0) or not np.array_equal(mol,mol0):
                raise ValueError('atom/molecule mapping changed')
            if is_unwrapped:
                xyz=rows[:,[col.index(c) for c in ('xu','yu','zu')]]
            else:
                xyz=rows[:,[col.index(c) for c in ('x','y','z')]] + edges * rows[:,[col.index(c) for c in ('ix','iy','iz')]]
            # Deterministic grouping by ascending molecular ID, then atom ID.
            by_mol=np.argsort(mol,kind='stable')
            snapshots.append(xyz[by_mol].reshape(-1,2,3))
            steps.append(step)
    if len(snapshots)<2: raise ValueError('at least two frames required')
    differences=np.diff(steps)
    if not np.all(differences>0) or not np.all(differences==differences[0]):
        raise ValueError('timesteps must be equally spaced')
    return AtomicH2Trajectory(np.stack(snapshots),box,timestep_fs*float(differences[0]),
                               length_unit='angstrom',time_unit='fs')
