# H₂ Molecular Dynamics Toolkit — Version 1.0

## Objective

Provide a reusable **end-to-end hydrogen MD workflow** rather than repeatedly adjusting a single number. The framework has five layers:

```text
Experiment configuration (JSON)
              │
              ▼
LAMMPS input generator ──► NPT equilibration ──► NVT production
              │                                    │
              │                                    ▼
              │                         Unwrapped H₂ atom trajectories
              │                                    │
              ▼                                    ▼
   Experiment metadata              Python molecular COM analysis
              │                                    │
              └──────────────────┬─────────────────┘
                                 ▼
                    Automated analysis artifacts
                   MSD, D, RDF, density, CI, QA
                                 │
                                 ▼
                 Campaign dashboard + CSV + report
                                 │
                                 ▼
                    GitHub-ready results and docs
```

## What already exists in this repository

- **Force field**: rigid two-site H₂ with published LJ site parameters from Adam et al. (2024), Table 1. The rigid 0.7414 Å bond is a documented separate approximation.
- **LAMMPS**: script-generated molecules and simulation boxes, NPT then NVT, 16-digit trajectory dumping, SHAKE bond constraints.
- **Analysis**: molecular center-of-mass MSD, multi-time-origin self-diffusion, RDF, density and EOS comparison.
- **Quality**: automated Python tests, source provenance, replica uncertainty, dump precision checks.
- **Visualization**: H₂ molecule animation (existing `h2md atomistic-analyze`), MSD plots and consolidated campaign dashboard.

## Three typical commands (existing run data)

Open Command Prompt in the repository root and install the local project:

```cmd
.venv\Scripts\activate.bat
python -m pip install -e ".[dev,validation]"
```

Then build one dashboard from all three completed system sizes:

```cmd
h2md-workflow dashboard --analyses lammps\paper_h2\validation_v05\analysis_research lammps\paper_h2\validation_v05\analysis_500 lammps\paper_h2\validation_v05\analysis_864 --out portfolio_report
```

Open `portfolio_report\dashboard.html` in Chrome or Edge. No LAMMPS execution, trajectories, or CoolProp calls are needed for this command; it reads the existing `summary.json` files.

Outputs: `dashboard.html` (self-contained interactive comparison), `RESULTS.md` (formatted results text), `comparison.csv` (Excel-friendly), and `comparison.json` (machine-readable).

## Optional: new simulation campaigns driven by one JSON file

The example configuration `examples/campaigns/h2_baseline.json` describes three sizes and common conditions. To generate all inputs without running anything:

```cmd
h2md-workflow plan --config examples\campaigns\h2_baseline.json --campaign campaigns\baseline
```

When ready to execute on a machine with LAMMPS:

```cmd
h2md-workflow run --campaign campaigns\baseline --lmp lmp
h2md-workflow analyze --campaign campaigns\baseline
h2md-workflow dashboard --analyses campaigns\baseline\h2_256\analysis campaigns\baseline\h2_500\analysis campaigns\baseline\h2_864\analysis --out portfolio_report
```

The `run` command skips an existing nonempty trajectory by default. `--force` explicitly allows a rerun and overwrite. Planning refuses to overwrite an existing campaign.

## Scientific boundaries

1. **Do not call the toolkit universally validated**: current EOS comparison at one temperature and pressure is a single-state density check, not a universal force-field validation.
2. **Do not claim finite-size independence** solely because intervals overlap. Record the sample count, uncertainty intervals and potential systematic biases.
3. **The H₂ rigid-bond model is classical**. Quantum rotational/vibrational effects are not reproduced.
4. **Replicate CIs** quantify between-seed statistical variation and not the systematic error of the force field, time step, or sampling choices.
5. **Different T and P states are not directly size comparisons**. The dashboard groups by state so the reader can see conditions.

## Suggested public GitHub layout

```text
README.md                       # short, visual overview + quick start
src/h2md/                      # reusable analysis and workflow Python package
lammps/paper_h2/               # hydrogen model and simulation inputs
examples/campaigns/           # configurable experiment manifest
examples/atomistic-demo/      # sample 3D animation (existing)
tests/                         # automated checks
assets/                        # figures
portfolio_report/              # compact dashboard and summary (optional)
docs/FINAL_FRAMEWORK.md        # architecture, scientific scope
```

Raw LAMMPS trajectories are often large: do not upload them to Git by default. Keep only scripts, small example data, summary CSV files, figures, and concise reports.

## Next feature roadmap (no more parameter tuning required)

1. Publish the existing 256/500/864 molecular-size campaign dashboard.
2. Add an **independent trajectory import interface** for more LAMMPS trajectory styles.
3. Support configurable **density, MSD, RDF, bond-length, and visualization** analysis selection in a single config file.
4. Provide a documented extension point for future molecules (CO₂, CH₄, H₂O) without pretending one force field covers all systems.

A new physical benchmark is optional, not a prerequisite to presenting a complete framework.
