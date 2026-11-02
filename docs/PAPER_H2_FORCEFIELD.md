# H₂ two-site force field — source provenance and assumptions

**Reference:** Adam, A. M., Bahamon, D., Al Kobaisi, M., & Vega, L. F. (2024). *Molecular dynamics simulations of the interfacial tension and the solubility of brine/H₂/CO₂ systems: Implications for underground hydrogen storage.* International Journal of Hydrogen Energy, **78**, 1344–1354. https://doi.org/10.1016/j.ijhydene.2024.06.319

## What is directly reported in the paper

Table 1, p. 1346: hydrogen is modeled using a **two-site Yang model**, with LJ parameters on each H interaction site: **σ = 2.7200 Å; ε = 0.0199 kcal/mol; q = 0 e**. Section 2.2: **12.8 Å** cutoff, **1.0 fs** timestep, periodic boundaries, and Nosé–Hoover control of temperature/pressure. The article covers H₂ mixtures with brine and CO₂; it does not itself establish accuracy for pure-H₂ diffusion coefficients computed by this tutorial.

## Implementation choices, NOT verbatim force-field values from the article

| Choice | Value | Reason / limitation |
|---|---|---|
| Fixed H–H separation | 0.7414 Å | Representative equilibrium H₂ bond length **assumed by this project**, not specified in paper Table 1; constrained with LAMMPS SHAKE |
| H mass | 1.00794 g/mol | Atomic mass used for each of two sites |
| Harmonic bond coefficient | 1000 kcal mol⁻¹ Å⁻² | Placeholder to define SHAKE's bond geometry; **not fitted or active as a flexible vibrational potential** during constrained MD |
| Initial cubic box | Ideal-gas estimate | Only a density guess for NPT preparation; not an H₂ EOS |
| NPT equilibration | 20 ps demo | **Much shorter** than paper's 2 ns equilibration; non-converged demo |
| NVT production | 50 ps demo | Deliberate change from source paper's NPT sampling, to keep the box fixed for the existing analysis reader |
| Production MSD | Molecular COM | Time-origin-averaged, Einstein D = slope(MSD)/6 in 3D; ensure diffusive limit and independent replicas |
| Cutoff scheme | LJ 12.8 Å, unshifted | Same nominal cutoff; long-range LJ corrections not used |
| Electrostatics | None | Both H interaction sites have q=0, so no Ewald/PPPM is necessary in an H₂-only system |

### Modeling scope

This is **not a strict reproduction** of the original study or proof of high-fidelity pure H₂ thermodynamics. It is a reproducible *H₂-only workflow incorporating the paper's published site parameters*. Before reporting predictive transport properties, validate density/EOS and self-diffusion against suitable pure H₂ literature or experiments, perform 2 ns+ equilibrations and multi-ns production, inspect MSD fit windows, test system-size and timestep effects, and repeat independent runs with uncertainty bars.

### LAMMPS potential

With intramolecular H-site LJ excluded by `special_bonds`, intermolecular H-site pairs interact as

\[U_{LJ}(r) = 4\epsilon [ (\sigma/r)^{12} - (\sigma/r)^{6} ], \qquad r < r_c.\]

The four site-site interactions between each pair of H₂ molecules are retained. SHAKE maintains an H–H distance per molecule. LAMMPS `units real` uses kcal/mol and Å, and pressure is in **atm**, so `P (MPa) × 9.869232667` is used in the input file.

The demo runs **NPT then NVT**; switching to fixed volume is important because the existing `h2md` trajectory parser intentionally rejects frames with changing simulation boxes. The NVT-production density comes from the last equilibrated NPT box, not the initial ideal-gas guess.

### References for implementation

- LAMMPS LJ pair potential: https://docs.lammps.org/pair_lj.html
- LAMMPS constrained bonds (`fix shake`): https://docs.lammps.org/fix_shake.html
- Source paper: https://doi.org/10.1016/j.ijhydene.2024.06.319
