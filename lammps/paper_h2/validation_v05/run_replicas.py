"""Prepare independent long H2 LAMMPS replica simulations; optional sequential execution.

The force-field parameters are imported from the version 0.4 paper-based generator.
No LAMMPS binary is bundled. Run with Python >=3.10.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]  # lammps/paper_h2
sys.path.insert(0, str(ROOT))
from build_h2_paper import (  # noqa: E402
    BOND_A, CUTOFF_A, EPS_KCAL_MOL, H_MASS, MPA_TO_ATM,
    SIGMA_A, TIMESTEP_FS, build_system, write_data,
)

PRESETS_PS = {
    'smoke': (20, 50),
    'standard': (500, 1000),
    'research': (2000, 3000),
}


def render_lammps_input(temp_k: float, pressure_mpa: float, seed: int,
                        eq_steps: int, prod_steps: int,
                        thermo_every: int, dump_every: int) -> str:
    if temp_k <= 0 or pressure_mpa <= 0:
        raise ValueError('temperature and pressure must be positive')
    if eq_steps < 2 * thermo_every or prod_steps < 2 * dump_every:
        raise ValueError('increase simulation steps or decrease output intervals')
    if not (1 <= seed < 900_000_000):
        raise ValueError('seed must be a positive LAMMPS-valid integer')
    return f'''# H2 v0.5: density (NPT) + multi-time-origin molecular diffusion (NVT)
# Site LJ: Adam et al., IJHE 78 (2024) Table 1, DOI 10.1016/j.ijhydene.2024.06.319
# H-H bond length {BOND_A} A is a separate modeling assumption, not in Table 1.
# For valid comparisons, verify NPT convergence and stability of fit windows.
clear
units real
atom_style molecular
boundary p p p
read_data h2_system.data
pair_style lj/cut {CUTOFF_A:.4f}
pair_coeff 1 1 {EPS_KCAL_MOL:.4f} {SIGMA_A:.4f} {CUTOFF_A:.4f}
pair_modify shift no
special_bonds lj 0.0 0.0 0.0
bond_style harmonic
bond_coeff 1 1000.0 {BOND_A:.4f}
neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes

timestep {TIMESTEP_FS:.1f}
thermo {thermo_every}
thermo_style custom step time temp press density pe ke etotal vol lx ly lz
thermo_modify flush yes
velocity all create {temp_k:.6f} {seed} mom yes rot yes dist gaussian
shell mkdir output
# IMPORTANT: SHAKE must be declared before a box-changing fix (NPT).
fix SH all shake 0.00001 200 0 b 1
fix EQ all npt temp {temp_k:.6f} {temp_k:.6f} 100.0 iso {pressure_mpa * MPA_TO_ATM:.8f} {pressure_mpa * MPA_TO_ATM:.8f} 1000.0
variable rho equal density
variable Tk equal temp
variable Patm equal press
variable V equal vol
# columns: timestep, density (g/cm^3), T (K), P (atm), volume (A^3)
fix NPTLOG all ave/time {thermo_every} 1 {thermo_every} v_rho v_Tk v_Patm v_V file output/npt_thermo.dat
run {eq_steps}
unfix NPTLOG
unfix EQ
# Avoid resetting image flags: dump xu/yu/zu coordinates must remain unwrapped.
reset_timestep 0
fix PROD all nvt temp {temp_k:.6f} {temp_k:.6f} 100.0
dump TRAJ all custom {dump_every} output/h2_production.lammpstrj id mol type xu yu zu
dump_modify TRAJ sort id
# %g defaults to six significant digits, which is too imprecise for
# unwrapped coordinates that grow to hundreds/thousands of Angstroms.
dump_modify TRAJ format float %.16g
run {prod_steps}
undump TRAJ
unfix PROD
unfix SH
write_data output/h2_final.data
'''


def prepare(output_dir: Path, preset: str, seeds: list[int], cells: int,
            temp_k: float, pressure_mpa: float,
            eq_ps: float | None = None, prod_ps: float | None = None,
            thermo_every: int = 500, dump_every: int = 500) -> dict:
    if not seeds or len(set(seeds)) != len(seeds):
        raise ValueError('seeds must be nonempty and unique')
    eq_ps = PRESETS_PS[preset][0] if eq_ps is None else eq_ps
    prod_ps = PRESETS_PS[preset][1] if prod_ps is None else prod_ps
    eq_steps = round(1000 * eq_ps / TIMESTEP_FS)
    prod_steps = round(1000 * prod_ps / TIMESTEP_FS)
    if abs(eq_ps*1000/TIMESTEP_FS-eq_steps) > 1e-6 or abs(prod_ps*1000/TIMESTEP_FS-prod_steps)>1e-6:
        raise ValueError('simulation duration must be integer multiple of timestep')
    if cells < 2: raise ValueError('cells must be >=2')
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for seed in seeds:
        path = output_dir / f'seed_{seed}'
        path.mkdir(parents=True,exist_ok=True)
        atoms, edge = build_system(cells=cells,temp_kelvin=temp_k,pressure_mpa=pressure_mpa,seed=seed)
        write_data(path/'h2_system.data',atoms,edge)
        (path/'in.h2_validation').write_text(render_lammps_input(temp_k,pressure_mpa,seed,
            eq_steps,prod_steps,thermo_every,dump_every),encoding='utf-8')
    config = {
        'version':'0.5.0', 'source_doi':'10.1016/j.ijhydene.2024.06.319',
        'model':'two-site LJ H2 with constrained H-H bond',
        'source_parameters': {'site_sigma_A':SIGMA_A,'site_epsilon_kcal_mol':EPS_KCAL_MOL,
                              'cutoff_A':CUTOFF_A,'site_charge_e':0.0},
        'assumptions': {'H_atom_mass_g_mol':H_MASS,'H_H_bond_A':BOND_A,
                        'npt_density':'NPT final fraction, not ideal gas initial density',
                        'diffusion':'NVT center-of-mass Einstein MSD, multiple time origins'},
        'molecules':4*cells**3,'temperature_K':temp_k,'pressure_MPa':pressure_mpa,
        'preset':preset,'seeds':seeds,'timestep_fs':TIMESTEP_FS,
        'equilibration_steps':eq_steps,'production_steps':prod_steps,
        'sample_thermo_every_steps':thermo_every,'dump_every_steps':dump_every,
        'density_reference':'CoolProp HEOS normal Hydrogen (Leachman et al. 2009)',
    }
    (output_dir/'experiment.json').write_text(json.dumps(config,indent=2)+'\n',encoding='utf-8')
    return config


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--preset',choices=PRESETS_PS,default='standard')
    p.add_argument('--output',default='runs_v05')
    p.add_argument('--seeds',type=int,nargs='+',default=[2027,2028,2029])
    p.add_argument('--cells',type=int,default=4)
    p.add_argument('--temperature-K',type=float,default=323.15)
    p.add_argument('--pressure-MPa',type=float,default=10.0)
    p.add_argument('--eq-ps',type=float)
    p.add_argument('--prod-ps',type=float)
    p.add_argument('--thermo-every',type=int,default=500)
    p.add_argument('--dump-every',type=int,default=500)
    p.add_argument('--run',action='store_true',help='run LAMMPS sequentially after generating files')
    p.add_argument('--lmp',default='lmp',help='LAMMPS executable, e.g., lmp.exe')
    args=p.parse_args(argv)
    config=prepare(Path(args.output),args.preset,args.seeds,args.cells,
       args.temperature_K,args.pressure_MPa,args.eq_ps,args.prod_ps,
       args.thermo_every,args.dump_every)
    print(f"Prepared {len(args.seeds)} replicas: {config['molecules']} H2 each; "
          f"{config['equilibration_steps']} NPT steps + {config['production_steps']} NVT steps.")
    if args.run:
        for seed in args.seeds:
            path=Path(args.output).resolve()/f'seed_{seed}'
            print(f'Running seed {seed} in {path}',flush=True)
            with (path/'lammps.log').open('w',encoding='utf-8') as log:
                subprocess.run([args.lmp,'-in','in.h2_validation'],cwd=path,
                               stdout=log,stderr=subprocess.STDOUT,check=True)
            print(f'Finished seed {seed}.',flush=True)
    else:
        print('To execute, rerun this command with --run (or run lmp -in in.h2_validation inside each seed folder).')

if __name__ == '__main__': main()
