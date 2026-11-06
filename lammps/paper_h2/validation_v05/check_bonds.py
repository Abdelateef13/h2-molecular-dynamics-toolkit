"""Audit saved H-H distances versus dump text precision for existing v0.5 runs.

Example: python check_bonds.py --runs runs_standard
Does not modify completed runs or trajectories. Does not require CoolProp.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

from h2md.atomistic import load_atomistic_lammps_dump
from analyze_replicas import bond_precision_audit, dump_coordinate_precision


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs', default='runs_standard')
    args = parser.parse_args()
    root = Path(args.runs)
    cfg = json.loads((root/'experiment.json').read_text(encoding='utf-8'))
    for seed in cfg['seeds']:
        folder = root / f'seed_{seed}'
        dump = folder/'output'/'h2_production.lammpstrj'
        tr = load_atomistic_lammps_dump(dump, cfg['timestep_fs'])
        digits = dump_coordinate_precision(folder/'in.h2_validation')
        bond = cfg.get('assumptions',{}).get('H_H_bond_A',0.7414)
        stats = bond_precision_audit(tr.atoms, float(bond), digits)
        print(f'Seed {seed} ({tr.n_frames} frames, {tr.n_molecules} molecules):')
        print(json.dumps(stats,indent=2))
    print('A rounding-consistent result is NOT independent verification of SHAKE.')

if __name__ == '__main__':
    main()
