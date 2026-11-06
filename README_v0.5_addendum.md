## Version 0.5: H₂ density benchmark and self-diffusion with uncertainty

The new [validation workflow](lammps/paper_h2/validation_v05/README.md) uses three independent LAMMPS seeds, extended NPT equilibration, fixed-volume NVT production, molecular center-of-mass mean-squared displacement and 95% Student-t intervals across independent runs. It benchmarks simulated density against the **CoolProp normal-hydrogen EOS** at 323.15 K and 10 MPa (defaults) and outputs reproducible numerical tables and plots. **No physical density or diffusion values are asserted until the LAMMPS simulations are completed and checked.**

[Windows step-by-step instructions](lammps/paper_h2/validation_v05/WINDOWS_V05.md). The H₂ site LJ parameters come from Adam et al. (2024), DOI:10.1016/j.ijhydene.2024.06.319; 0.7414 Å constrained bond length is a separately declared modeling assumption.
