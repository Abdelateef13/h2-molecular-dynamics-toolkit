# Upgrade guide — Windows / GitHub Desktop (Version 0.5)

**Keep your existing repository** `C:\Users\ADMIN\Desktop\h2-molecular-dynamics-toolkit` and its `.git` folder. You only need to merge the Version 0.5 add-on files into it.

## 1. Add the files

1. Extract the Version 0.5 add-on ZIP.
2. Copy its `lammps\paper_h2\validation_v05` folder into the **same path** in the existing repository.
3. Copy `tests\test_validation_v05.py` into the existing `tests` folder.
4. Replace `pyproject.toml` and `.gitignore` from the add-on (the full v0.5 distribution also includes the earlier project files).
5. Your old v0.4 scripts, earlier outputs and `.git` directory remain untouched.

## 2. Install Python extras

Open File Explorer in your project folder; click the address bar, enter `cmd`, press Enter.

```bat
.venv\Scripts\activate.bat
python -m pip install -e ".[dev,validation]"
python -m pytest -q
lmp -h
```

If `lmp -h` is not recognized, correct the LAMMPS executable location before proceeding. If you use a different executable name, specify `--lmp "C:\path\to\lmp.exe"` in later commands.

## 3. Smoke test (first)

```bat
cd lammps\paper_h2\validation_v05
python run_replicas.py --preset smoke --output runs_smoke --run
python analyze_replicas.py --runs runs_smoke --out analysis_smoke
```

There are three seeds. `--run` performs actual LAMMPS simulations sequentially. If you omit it, the scripts only generate input files and analysis will say the dumps are missing. **Do not interpret smoke diffusion values scientifically.**

Open `analysis_smoke\REPORT.md` and `analysis_smoke\density_validation.png` and `analysis_smoke\diffusion_msd.png`.

## 4. Long sampling with three independent seeds

```bat
python run_replicas.py --preset standard --output runs_standard --run
python analyze_replicas.py --runs runs_standard --out analysis_standard
```

Standard: 500 ps NPT and 1 ns NVT for each of the 3 seeds. Data are automatically written to `runs_standard\seed_2027\output` and analogous folders. Windows Command Prompt stays attached until simulation completes. The script saves a `lammps.log` under each run.

Research option (more expensive):

```bat
python run_replicas.py --preset research --output runs_research --run
python analyze_replicas.py --runs runs_research --out analysis_research
```

Research preset: 2 ns NPT and 3 ns NVT for each seed. Both are *starting protocols*, not guarantees of convergence.

## 5. How to interpret outputs

Inspect the **whole density trace**: it should settle around a stationary mean rather than drift throughout NPT. Density is compared against normal-H₂ CoolProp EOS for your requested T and P (and per seed at measured T/P). A difference can come from inadequate equilibration, the model, temperature/pressure mismatch, or finite-size effects.

Inspect the **MSD curve**: use only the long-time linear diffusive regime. The dotted fit lines correspond to the default 10–50% of recorded NVT duration. The reported confidence interval is built from independent *seed-level* estimates, not the thousands of highly correlated frames or overlapping origins. With three seeds, uncertainty may be large.

`REPORT.md` flags certain diagnostics, but numerical warning checks cannot by themselves establish physical validity. Use `--fit-start-frac` and `--fit-end-frac` with `analyze_replicas.py` to test window robustness. For example:

```bat
python analyze_replicas.py --runs runs_standard --out analysis_altwindow --fit-start-frac 0.15 --fit-end-frac 0.45
```

## 6. Update GitHub

Open GitHub Desktop -> your existing repo -> review the **new scripts and docs**, not generated large results. Commit: `Add H2 density validation and replica diffusion statistics`. Click **Push origin**. In the README, prominently point to `lammps/paper_h2/validation_v05/README.md`.

**Note:** the original PDF supplies the H₂ site Lennard–Jones parameters, not a reference pure-gas density table or an H₂ self-diffusion validation dataset. CoolProp provides an independent density EOS. Never describe the reference EOS as a value reported in the paper.
