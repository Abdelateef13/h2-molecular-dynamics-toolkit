# Upgrade an existing H₂ MD Toolkit checkout to Version 1.0 — Windows

You do **not** need to rerun your existing LAMMPS simulations.

## 1. Back up and merge

1. Download `H2_Molecular_Dynamics_Toolkit_v1.0_ADDON.zip`.
2. Extract the ZIP on your Desktop.
3. Copy the **contents** of the extracted add-on into your existing `h2-molecular-dynamics-toolkit` folder.
4. In Windows Explorer, choose to **merge directories** and **replace matching files**.
5. Keep your existing `.git` folder, `.venv`, and simulation results (`runs_research`, `runs_500`, `runs_864`, `analysis_research`, `analysis_500`, `analysis_864`). The add-on does not contain or delete those data.

## 2. Reinstall the local package

From a Command Prompt in your existing repository folder:

```cmd
.venv\Scripts\activate.bat
python -m pip install -e ".[dev,validation]"
python -m pytest -q
```

## 3. Generate the real dashboard from your completed simulations

Stay in the **root** folder of the repository and run this one command:

```cmd
h2md-workflow dashboard --analyses lammps\paper_h2\validation_v05\analysis_research lammps\paper_h2\validation_v05\analysis_500 lammps\paper_h2\validation_v05\analysis_864 --out portfolio_report
```

Open `portfolio_report\dashboard.html` in Edge or Chrome. Also inspect `portfolio_report\RESULTS.md` and `portfolio_report\comparison.csv`.

The add-on contains a **preview** dashboard made from the numerical summaries you already provided. Rerunning the command above regenerates this from **your actual analysis folders**.

## 4. Explore other commands

```cmd
h2md-workflow --help
h2md-workflow plan --help
h2md-workflow run --help
h2md-workflow analyze --help
h2md-workflow dashboard --help
```

You can reuse the framework for future temperature-pressure studies by editing `examples\campaigns\h2_baseline.json` and generating a new campaign. Do not run another campaign just to finish presenting this one.

## 5. Add to GitHub Desktop

Open GitHub Desktop and select this repository. Review the changes carefully (especially the `.gitignore`). Commit with:

`Release v1.0: end-to-end H2 MD workflow and interactive dashboard`

Click **Push origin**. On GitHub, pin the repository to your profile if it is not pinned already.

### Important

- Only the original two-site LJ nonbonded parameters were taken from the uploaded Adam et al. paper. The fixed H–H bond is a separate documented modeling assumption.
- Large trajectories should not be committed to Git; keep the lightweight comparison report and sample visualizations.
- The dashboard is an analysis interface, **not a claim of full scientific validation**.
