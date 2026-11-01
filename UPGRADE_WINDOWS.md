# Windows upgrade guide: v0.1 → v0.2

You already published `h2-molecular-dynamics-toolkit` on GitHub. **Do not create another GitHub repository.** You will update its local folder and publish the changes.

## A. Update the existing project files

1. Download and unzip `H2_Molecular_Dynamics_Toolkit_v0.2.zip` somewhere convenient on Desktop.
2. Open your **existing published project's local folder** (the same folder you originally chose in GitHub Desktop).
3. Make a backup copy of that folder, just in case.
4. Open the NEW v0.2 ZIP's extracted `h2-md-v0.2` folder in another File Explorer window.
5. Copy these new/updated items **into the existing repo root**, replacing duplicates when asked:
   - `src` (copy/merge; do not delete the old `src`)
   - `tests` (merge)
   - `assets` (merge)
   - `examples` (merge)
   - `README.md`
   - `pyproject.toml`
   - `UPGRADE_WINDOWS.md`
   - `LICENSE` (same license)
6. Never delete the hidden `.git` directory of the existing repository.

## B. Run the simulation on Windows

1. In the existing repository folder, click the File Explorer address bar, type `cmd`, and press Enter.
2. Activate the existing virtual environment, if you created one: 
   ```cmd
   .venv\Scripts\activate.bat
   ```
   If `.venv` does not exist, first run `py -m venv .venv`, then activate it.
3. Update dependencies in this same folder:
   ```cmd
   python -m pip install -e ".[dev]"
   ```
   An internet connection is needed the first time to install packages.
4. Verify tests:
   ```cmd
   python -m pytest -q
   ```
   Expected: **11 passed**.
5. Run a full two-atom H₂ demo:
   ```cmd
   h2md atomistic-demo --output results/atomistic
   ```
6. Open `results\atomistic\h2_interactive.html` by double-clicking it. Use the Play control, drag to rotate, and scroll to zoom. This page is self-contained; it does not need an internet connection.
7. Open these PNGs for static figures:
   - `results\atomistic\h2_atomistic_3d.png`
   - `results\atomistic\h2_bonds.png`
   - `results\atomistic\h2_analysis.png`
   - `results\atomistic\h2_diagnostics.png`
8. (Optional) Install [OVITO Basic](https://www.ovito.org/) and open `results\atomistic\h2_molecules.xyz` to view the time series. If the H–H bonds are not automatically displayed, apply a *Create bonds* modifier with a small reduced-unit distance cutoff of **0.5σ**. **Be aware that the XYZ coordinates are in reduced units, not Å.** The supplied HTML visualization is an easier first choice.

## C. Publish version 0.2

1. Open **GitHub Desktop** and select your existing `h2-molecular-dynamics-toolkit` repository.
2. Click **Changes**. You should see new Python files and modified README.
3. Type commit summary: `Add two-atom H2 molecular MD and interactive 3D visualization`.
4. Click **Commit to main**.
5. Click **Push origin**.
6. Refresh your repository page and check the **Actions** tab for a passing test run.

If `h2md` is not recognized, run instead:

```cmd
python -m h2md.cli atomistic-demo --output results/atomistic
```

## Scientific interpretation

- `atomistic-demo` models **two H atoms per H₂**, not a single site per molecule.
- H–H bonds are classical harmonic springs, and non-bonded H sites use a simple Lennard–Jones potential (with intramolecular LJ excluded).
- This is **an educational reduced-unit force field**, **not** a validated description of hydrogen at a chosen pressure or temperature. Do not publish its diffusion values as real H₂ transport coefficients.
- Molecular diffusion uses each molecule's **unwrapped center of mass**. Measuring the MSD of individual H atoms is not the same thing.

For actual validated H₂ calculations, use a literature-validated hydrogen intermolecular potential, carefully chosen thermodynamic conditions, sufficient sampling, and independent replication; the v0.2 `atomistic-analyze` command can import two-atom LAMMPS dumps to analyze the result.
