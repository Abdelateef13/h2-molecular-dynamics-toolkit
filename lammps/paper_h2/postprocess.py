"""Validate a completed H2 LAMMPS trajectory and generate v0.2 analysis reports.

Run *from lammps/paper_h2* after LAMMPS, with project installed:
    python postprocess.py
"""
from __future__ import annotations
import argparse, json, subprocess, sys
from pathlib import Path
import numpy as np
from h2md.atomistic import load_atomistic_lammps_dump


def diagnose(path: Path, timestep_fs: float=1.0):
    tr=load_atomistic_lammps_dump(path, timestep_fs=timestep_fs)
    b=tr.bond_lengths()
    metrics={
       'frames': tr.n_frames, 'n_h2_molecules':tr.n_molecules,
       'n_h_sites':tr.n_molecules*2,'box_angstrom':float(tr.box_length),
       'bond_mean_angstrom':float(np.mean(b)),
       'bond_std_angstrom':float(np.std(b)),
       'bond_max_abs_error_angstrom':float(np.max(np.abs(b-0.7414))),
       'frame_dt_fs':float(tr.frame_dt),
       'minimum_image_condition_2rc_less_than_L':bool(2*12.8 < tr.box_length),
       'bond_tolerance_condition_max_abs_error_less_than_0p01_A':bool(np.max(np.abs(b-0.7414))<0.01)
    }
    return metrics


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dump', default='output/h2_paper.lammpstrj')
    p.add_argument('--out', default='output/analysis')
    p.add_argument('--timestep-fs', type=float, default=1.0)
    p.add_argument('--no-plots', action='store_true')
    args=p.parse_args(argv)
    path=Path(args.dump)
    if not path.exists():
        p.error(f'LAMMPS trajectory not found: {path} (run lmp -in in.h2_paper first)')
    metrics=diagnose(path,args.timestep_fs)
    out=Path(args.out)
    out.mkdir(parents=True,exist_ok=True)
    (out/'validation.json').write_text(json.dumps(metrics,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(metrics,indent=2))
    if not metrics['minimum_image_condition_2rc_less_than_L'] or not metrics['bond_tolerance_condition_max_abs_error_less_than_0p01_A']:
        raise SystemExit('VALIDATION FAILED: review box size, SHAKE tolerance and input; not plotting.')
    if not args.no_plots:
        subprocess.run([sys.executable,'-m','h2md.cli','atomistic-analyze',
           '--input',str(path),'--timestep-fs',str(args.timestep_fs),
           '--output',str(out)],check=True)
        print('Open',out/'h2_interactive.html')

if __name__=='__main__':main()
