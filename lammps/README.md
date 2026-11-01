# LAMMPS two-atom H₂ simulation demo

This folder provides a **working-input example** for the existing `h2md atomistic-analyze` command. It uses LAMMPS as the MD engine and the Python package for post-processing. LAMMPS must be installed separately.

- `build_h2_data.py` — deterministically creates `h2_system.data` (32 molecules, 64 atoms by default).
- `h2_system.data` — ready-to-read initial coordinates, molecule IDs and bonds.
- `in.h2_demo` — LAMMPS NVT equilibration + NVE production, trajectory output with atom IDs, molecule IDs and unwrapped Å coordinates.
- `output/` — generated trajectory; ignored by Git.

Run from this directory:

```cmd
python build_h2_data.py
lmp -in in.h2_demo
```

**Units:** `real` -> distances Å, timestep fs, energy kcal/mol, temperature K. The input uses `timestep 0.1` and a dump stride of 100 steps, so the correct toolkit argument is `--timestep-fs 0.1` (not 10). The reader multiplies the timestep by the dump stride to get a 10 fs frame interval.

**LAMMPS potential:** The input uses `bond_style harmonic` with `bond_coeff 1 350.0 0.7414` and `pair_style lj/cut` with `pair_coeff 1 1 0.016 2.9 7.5`. These are **illustrative parameters, NOT a validated H₂ model**. Intramolecular LJ interactions are excluded using `special_bonds lj 0.0 0.0 0.0`. An H₂-specific literature potential, consistent thermodynamic state and convergence studies are needed for physical predictions.

LAMMPS documentation:
- https://docs.lammps.org/read_data.html
- https://docs.lammps.org/dump.html
- https://docs.lammps.org/bond_harmonic.html
- https://docs.lammps.org/Install_windows.html
