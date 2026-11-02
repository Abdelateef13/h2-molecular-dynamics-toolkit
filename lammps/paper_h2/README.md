# LAMMPS pure H₂ — published two-site site-site LJ (v0.4)

See `../../docs/PAPER_H2_FORCEFIELD.md` for exactly which parameters come from the paper and which are implementation assumptions. `build_h2_paper.py` generates both a periodic H₂ data file and matching LAMMPS input so temperature/pressure are consistent.

```bat
python build_h2_paper.py
lmp -in in.h2_paper
python postprocess.py
```

`postprocess.py` requires the main H₂ package installed with `pip install -e ".[dev]"` from the repository root. The existing parser supports a fixed volume during sampled production, so NPT is used only during preparation.
