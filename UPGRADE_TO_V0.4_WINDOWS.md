# Upgrade your GitHub H₂ MD Toolkit to v0.4 (Windows)

This is an **additive** upgrade. It adds a folder `lammps/paper_h2`, documentation, and tests. Your v0.3 files remain intact.

## Step 1: Copy the folders

Extract the ZIP. Copy **the contents** of `h2-paper-v04` into the root of the local repository (`h2-molecular-dynamics-toolkit`, the folder with `pyproject.toml`). Merge `lammps/` and `tests/` but do not delete any existing files or `.git`.

## Step 2: Install dependencies

In the repository folder, File Explorer address bar -> `cmd`, Enter:

```cmd
.venv\Scripts\activate.bat
python -m pip install -e ".[dev]"
python -m pytest -q
```

If `.venv` is missing, create it first with `py -m venv .venv`.

## Step 3: Generate the force-field-based H₂ system

```cmd
cd lammps\paper_h2
python build_h2_paper.py
```

Default: 256 molecules / 512 H sites, 323.15 K, 10 MPa target, 1 fs timestep, 20 ps NPT preparation, 50 ps NVT production. The initial box comes from ideal-gas scaling, **not** an H₂ equation of state.

The builder writes `h2_system.data`, `in.h2_paper`, and `model_config.json`. It is safe to re-run; it overwrites those generated files.

## Step 4: Run LAMMPS

```cmd
lmp -in in.h2_paper
```

If `lmp` isn't found, use the actual LAMMPS executable name from your installation, e.g. `lmp.exe`, with an absolute path. The required LAMMPS packages include MOLECULE (bond topology) and RIGID (`fix shake`).

Outputs: `output\h2_paper.lammpstrj` and `output\h2_final.data`. Note that the demo is short; **do not present this D value as converged research data**.

## Step 5: Validate and analyze

```cmd
python postprocess.py
```

After analysis, open:

```text
lammps\paper_h2\output\analysis\h2_interactive.html
```

Additional outputs: MSD, diffusion fit, COM RDF, H–H bond distances, a 3D static image, an XYZ trajectory, and `validation.json`. Bond length near 0.7414 Å and final box size greater than 25.6 Å are checked automatically.

## Optional: Target different conditions

For a 368.15 K, 30 MPa system (start with at least 256 molecules):

```cmd
python build_h2_paper.py --temperature-K 368.15 --pressure-MPa 30 --cells 4
lmp -in in.h2_paper
python postprocess.py
```

For 60 MPa, a larger box is necessary because 12.8 Å must be less than half its periodic edge:

```cmd
python build_h2_paper.py --temperature-K 368.15 --pressure-MPa 60 --cells 5
lmp -in in.h2_paper
python postprocess.py
```

**For longer runs**, change the steps in the builder, keeping the 1 fs timestep:

```cmd
python build_h2_paper.py --eq-steps 2000000 --prod-steps 3000000 --dump-every 2000
```

That requests 2 ns NPT + 3 ns NVT, still **different from the publication's NPT production ensemble**. Longer runs can consume substantial CPU and disk resources.

## Step 6: Publish

Open GitHub Desktop; check the files under Changes, commit with `Add published two-site H2 force field and SHAKE simulation`, and press **Push origin**.

Keep your generated LAMMPS output from being accidentally committed if it's large. A local `.gitignore` in `lammps/paper_h2` excludes generated trajectories and data files; generated `in.h2_paper` and configuration are included as small reproducible examples.
