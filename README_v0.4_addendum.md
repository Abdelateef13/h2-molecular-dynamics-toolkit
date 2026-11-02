## Version 0.4 — Published H₂ two-site force field (LAMMPS)

A reproducible **H₂-only** simulation workflow using the two-site H₂ nonbonded parameters reported in [Adam et al., International Journal of Hydrogen Energy 78 (2024), 1344–1354](https://doi.org/10.1016/j.ijhydene.2024.06.319), Table 1.

- H-site LJ **σ = 2.7200 Å**, **ε = 0.0199 kcal/mol**, charge = 0; cutoff **12.8 Å**, timestep **1 fs**.
- Explicit 2-site H₂ molecules with a **fixed 0.7414 Å H–H bond (assumption, not in source Table 1)** using LAMMPS SHAKE.
- Short **NPT preparation → NVT production**; the ensemble and demo durations differ from the source paper and must not be treated as a strict reproduction.
- Postprocessing generates center-of-mass MSD/diffusion and RDF, H–H bond statistics, interactive HTML, and an OVITO-compatible trajectory.
- Includes validation checks and tests. See [v0.4 Windows guide](UPGRADE_TO_V0.4_WINDOWS.md), [parameter provenance](docs/PAPER_H2_FORCEFIELD.md), and [LAMMPS input](lammps/paper_h2/in.h2_paper).

```bat
cd lammps\paper_h2
python build_h2_paper.py
lmp -in in.h2_paper
python postprocess.py
```

This is a model *based on published nonbonded parameters*, not a benchmarked prediction of pure-H₂ thermodynamic or transport properties. Run longer simulations, multiple seeds, and independent validation before scientific interpretation.
