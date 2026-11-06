# Version 0.5 — Independent replicas, density benchmark and H₂ self-diffusion

This workflow extends the **v0.4 paper-based, rigid two-site hydrogen LAMMPS model**. It does not replace the original educational Python demonstrators. It uses **exactly the same nonbonded H₂ site LJ parameters** documented in Adam et al. (2024), Table 1: `sigma = 2.7200 Å`, `epsilon = 0.0199 kcal/mol`, `charge = 0`. The rigid 0.7414 Å H–H bond remains a *separate assumption*. LAMMPS `real` units, time step 1 fs, LJ cutoff 12.8 Å.

**Reference:** Adam, Bahamon, Al Kobaisi & Vega (2024), *Molecular dynamics simulations of the interfacial tension and the solubility of brine/H₂/CO₂ systems: Implications for underground hydrogen storage*, *International Journal of Hydrogen Energy* 78, 1344–1354, DOI: [10.1016/j.ijhydene.2024.06.319](https://doi.org/10.1016/j.ijhydene.2024.06.319).

## Why split into NPT and NVT?

- **NPT equilibration:** let the box volume adapt to target `T` and `P`; collect the time series of `density`, `temperature`, `pressure`, `volume`; estimate the mean density only from the **last 50%** of samples (configurable). Plot the entire trace and inspect whether it is stationary. Real-unit `density` is g/cm³, converted to kg/m³ for comparison.
- **NVT production:** fix the equilibrated box volume so the existing cubic, constant-cell unwrapped `xu yu zu` trajectory reader is valid. Calculate molecular H₂ **center-of-mass self-diffusion**, not site-wise diffusion.
- **Three independent replicas:** initial orientations and Maxwell velocities use different random seeds. Report each estimate separately and use *across-replica* Student-t 95% confidence intervals (very wide with just three runs). Correlated trajectory samples and overlapping time origins **are not independent replicas**.

## Density reference

`CoolProp.CoolProp.PropsSI('Dmass','T',T_K,'P',P_MPa * 1e6,'HEOS::Hydrogen')`

Returns **kg/m³** for **normal H₂** using the CoolProp equation of state based on Leachman et al. (2009). This is a thermodynamic density **reference** (not an experimental self-diffusion reference), and quantum/spin effects in real H₂ are not explicitly represented in the rigid classical two-site model. We show both density at the target T/P and per-run reference density at the *measured* mean T/P.

References: [CoolProp Hydrogen fluid](https://coolprop.org/fluid_properties/fluids/Hydrogen.html) and [CoolProp PropsSI API](https://coolprop.org/coolprop/HighLevelAPI.html).

## Run on Windows

From the repository root (not from this directory), with a virtual environment and LAMMPS available as `lmp`:

```bat
.venv\Scripts\activate.bat
python -m pip install -e ".[dev,validation]"
cd lammps\paper_h2\validation_v05

REM Test the complete workflow first (20 ps NPT + 50 ps NVT per seed)
python run_replicas.py --preset smoke --output runs_smoke --run
python analyze_replicas.py --runs runs_smoke --out analysis_smoke

REM Longer exploratory production (500 ps NPT + 1000 ps NVT per seed)
python run_replicas.py --preset standard --output runs_standard --run
python analyze_replicas.py --runs runs_standard --out analysis_standard

REM Extended protocol (2000 ps NPT + 3000 ps NVT per seed)
python run_replicas.py --preset research --output runs_research --run
python analyze_replicas.py --runs runs_research --out analysis_research
```

Windows detailed guide: [WINDOWS_V05.md](WINDOWS_V05.md).

| Preset | NPT per seed | NVT per seed | Purpose |
|---|---:|---:|---|
| `smoke` | 20 ps | 50 ps | Only verify setup; **not reliable transport science** |
| `standard` | 500 ps | 1 ns | Initial convergence/exploratory study |
| `research` | 2 ns | 3 ns | Stronger sampling; **still requires convergence tests** |

Default: three seeds `2027 2028 2029`, 256 molecules / 512 atoms, 323.15 K and 10 MPa, output trajectory every 0.5 ps. Adjust with `--seeds`, `--cells`, `--temperature-K`, `--pressure-MPa`, `--eq-ps`, `--prod-ps`, `--dump-every`.

`run_replicas.py` creates independent per-seed LAMMPS data and input files. `--run` runs sequentially using LAMMPS executable `lmp`. **Without `--run`, it only prepares the files.** If your executable has another name: `--lmp "C:\\path\\to\\lmp.exe"`. Generated runs are excluded from Git.

## Analysis

The Einstein relation is fitted to molecular-center displacement:

\[ \langle |\mathbf{R}(t_0 + t) - \mathbf{R}(t_0)|^2 \rangle = 6 D t + c. \]

Unwrapped molecular COM trajectories are used, with multiple overlapping start times `t0`, default every ten saved frames. We remove framewise global translational drift. A linear fit is made over **10–50% of the recorded production span**; the alternative **15–40%** window is compared as a basic sensitivity screen. Lag curves are generated to 50% of the production span (to retain averaging over origins).

For unit conversion: with `t` in ps and MSD in Å², `D (m²/s) = (MSD slope / 6) × 10⁻⁸`.

The **mean ± 95% confidence interval** is calculated using Student's t distribution on the `n=3` (or more) *independent run-level D estimates*. A confidence interval that spans below zero is possible for imprecise small-sample estimates, and must not be interpreted as negative physical diffusivity.

## Outputs

`analysis_standard/` or your chosen output:

- `REPORT.md` — human-readable summary, EOS reference, means, 95% CIs and warnings.
- `summary.json` — complete machine-readable results, provenance and analysis choices.
- `replicas.csv` — per-seed pressure, temperature, density, EOS error, diffusion, fit statistics.
- `msd_curves.csv` — numeric MSD curves.
- `density_validation.png` — density vs NPT time for all replicas, reference EOS line.
- `diffusion_msd.png` — multi-origin COM MSD and fit lines for each seed.
- `replica_comparison.png` — density and diffusion by independent seed.

**Do not claim validation solely from a small density error.** Before publishing quantitative D values: establish NPT stationarity, diffusion fit-window stability, finite-size sensitivity (e.g., test 256 vs 500/864 molecules with suitable box), run-length convergence, thermostat influence and additional replicas. The current package does **not** correct finite-size effects or calibrate diffusion against experimental reference values.

## Reproducibility and code review

- The generated `experiment.json` records the provenance, seeds, thermodynamic conditions, sample intervals, and run lengths.
- `fix shake` precedes `fix npt` to avoid the LAMMPS ordering error encountered in v0.4.
- This workflow is compatible by design with LAMMPS 15Jun2023 directives, but its **actual LAMMPS execution must be tested on your machine**; automated tests cover generator logic, parsing, units, formulae and mocked integration only.
- No test datasets should be described as simulation observations.
