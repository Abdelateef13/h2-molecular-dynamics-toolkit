"""Build an H2-only two-site LAMMPS model using Adam et al. (2024) LJ sites.

Source for nonbonded terms: DOI 10.1016/j.ijhydene.2024.06.319,
Table 1, H_H2; cutoff & timestep: Section 2.2.

The fixed intramolecular distance (0.7414 A) is a documented physical
assumption, NOT a bond length reported in the source paper's Table 1.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np

AVOGADRO = 6.02214076e23
R_J_MOL_K = 8.31446261815324
H_MASS = 1.00794  # g/mol, each H atom; modeling choice
SIGMA_A = 2.7200  # Adam et al. (2024), Table 1
EPS_KCAL_MOL = 0.0199  # Adam et al. (2024), Table 1
CHARGE_E = 0.0  # Adam et al. (2024), Table 1
CUTOFF_A = 12.8  # Adam et al. (2024), Section 2.2
BOND_A = 0.7414  # approximate H2 equilibrium geometry; NOT from Table 1
TIMESTEP_FS = 1.0  # Adam et al. (2024), Section 2.2
MPA_TO_ATM = 9.869232667160128


def fcc_centers(cells: int, length: float):
    if cells < 2: raise ValueError('cells >= 2 required')
    if length <= 0: raise ValueError('box length must be positive')
    basis = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]], float)
    spacing = length / cells
    return np.array([((np.array([i,j,k],float)+b+.25)*spacing)
        for i in range(cells) for j in range(cells) for k in range(cells) for b in basis])


def ideal_gas_box_angstrom(n_molecules: int, pressure_mpa: float, temp_kelvin: float):
    """P V = n R T, used ONLY to seed a box prior to NPT equilibration."""
    if n_molecules <= 0 or pressure_mpa <= 0 or temp_kelvin <= 0:
        raise ValueError('n, P and T must be positive')
    volume_m3 = (n_molecules / AVOGADRO) * R_J_MOL_K * temp_kelvin / (pressure_mpa * 1e6)
    return (volume_m3 * 1e30)**(1/3)


def build_system(cells: int=4, temp_kelvin: float=323.15,
                 pressure_mpa: float=10.0, seed: int=2026):
    """Return (Nmolecules,2,3) atom positions and box edge (A)."""
    n = 4*cells**3
    length = ideal_gas_box_angstrom(n, pressure_mpa, temp_kelvin)
    # Reserve margin: box may shrink during equilibration. Check box again
    # post-NPT; L > 2*rc MUST hold in every frame for the minimum image convention.
    if length <= 2*CUTOFF_A + 4.0:
        raise ValueError(f'initial box {length:.2f} A too small for 12.8 A cutoff; increase --cells (try 5)')
    center = fcc_centers(cells, length)
    rng=np.random.default_rng(seed)
    direction = rng.normal(size=(n,3))
    direction /= np.linalg.norm(direction, axis=1)[:,None]
    atomcoords = np.stack([center - .5*BOND_A*direction,
                           center + .5*BOND_A*direction], axis=1)
    assert np.all(atomcoords >= 0) and np.all(atomcoords < length)
    return atomcoords, length


def write_data(path: Path, positions: np.ndarray, length: float):
    path = Path(path)
    n = len(positions)
    lines = [
       'H2 two-site model: H-H LJ parameters Adam et al. (2024) Table 1',
       '', f'{2*n} atoms', f'{n} bonds', '', '1 atom types', '1 bond types',
       '', f'0.0 {length:.12f} xlo xhi', f'0.0 {length:.12f} ylo yhi',
       f'0.0 {length:.12f} zlo zhi', '', 'Masses', '', f'1 {H_MASS:.8f}',
       '', 'Atoms # molecular', ''
    ]
    for m, pair in enumerate(positions, 1):
        for side, pos in enumerate(pair):
            aid = 2*(m-1)+side+1
            lines.append(f'{aid} {m} 1 {pos[0]:.10f} {pos[1]:.10f} {pos[2]:.10f}')
    lines.extend(['', 'Bonds', ''])
    lines.extend(f'{m} 1 {2*m-1} {2*m}' for m in range(1,n+1))
    path.write_text('\n'.join(lines)+'\n',encoding='utf-8')


def render_input(temp_kelvin: float, pressure_mpa: float, eq_steps: int,
                 prod_steps: int, dump_every: int):
    """LAMMPS input for NPT preparation then NVT production at equilibrated box."""
    if eq_steps < 1 or prod_steps < dump_every or dump_every < 1:
        raise ValueError('positive steps required; production >= dump interval')
    atm=pressure_mpa * MPA_TO_ATM
    return f'''# H2 two-site LJ model; Adam et al., IJHE 78 (2024) 1344-1354
# https://doi.org/10.1016/j.ijhydene.2024.06.319  Table 1 + Sec. 2.2
# WARNING: H-H length 0.7414 A is assumed, NOT reported in Table 1.
# NPT equilibration -> NVT diffusion sampling; NOT a complete reproduction
# of the paper's NPT production/mixture simulations.
# From this directory: lmp -in in.h2_paper
clear
units           real
atom_style      molecular
boundary        p p p
read_data       h2_system.data

# LJ site-site H2 parameters from Table 1; q(H2)=0 so no Coulomb style.
pair_style      lj/cut {CUTOFF_A:.4f}
pair_coeff      1 1 {EPS_KCAL_MOL:.4f} {SIGMA_A:.4f} {CUTOFF_A:.4f}
pair_modify     shift no
special_bonds   lj 0.0 0.0 0.0

# Dummy harmonic stiffness provides equilibrium geometry to SHAKE.
# H-H is CONSTRAINED by SHAKE in MD; 1000 is not a fitted bond constant.
bond_style      harmonic
bond_coeff      1 1000.0 {BOND_A:.4f}
neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes

timestep        {TIMESTEP_FS:.1f}         # fs
thermo          500
thermo_style    custom step time temp press density pe ke etotal vol lx ly lz
thermo_modify   flush yes
velocity        all create {temp_kelvin:.4f} 740127 mom yes rot yes dist gaussian

# Nose-Hoover NPT to target P, T. Pressure in 'real' units is atm.
fix             SH all shake 0.00001 200 0 b 1
fix             EQ all npt temp {temp_kelvin:.4f} {temp_kelvin:.4f} 100.0 iso {atm:.7f} {atm:.7f} 1000.0
run             {eq_steps}
unfix           EQ

# Freeze the equilibrated volume by switching NPT to NVT for analysis.
# Existing h2md reader requires a fixed cubic box during saved frames.
reset_timestep  0
fix             PROD all nvt temp {temp_kelvin:.4f} {temp_kelvin:.4f} 100.0
shell           mkdir output
dump            TR all custom {dump_every} output/h2_paper.lammpstrj id mol type xu yu zu
dump_modify     TR sort id
run             {prod_steps}
undump          TR
unfix           PROD
unfix           SH
write_data      output/h2_final.data
'''


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cells', type=int, default=4,help='4*cells^3 H2 molecules (default 256)')
    p.add_argument('--temperature-K', type=float, default=323.15)
    p.add_argument('--pressure-MPa', type=float, default=10.0)
    p.add_argument('--seed', type=int, default=2026)
    p.add_argument('--eq-steps', type=int, default=20000,help='NPT steps (20 ps default)')
    p.add_argument('--prod-steps', type=int, default=50000,help='NVT steps (50 ps default)')
    p.add_argument('--dump-every', type=int, default=500)
    p.add_argument('--output-dir', default='.')
    args=p.parse_args(argv)
    pos, length=build_system(args.cells,args.temperature_K,args.pressure_MPa,args.seed)
    dest=Path(args.output_dir)
    dest.mkdir(parents=True,exist_ok=True)
    write_data(dest/'h2_system.data',pos,length)
    (dest/'in.h2_paper').write_text(render_input(args.temperature_K,args.pressure_MPa,
        args.eq_steps,args.prod_steps,args.dump_every),encoding='utf-8')
    meta = {
       'source_doi':'10.1016/j.ijhydene.2024.06.319',
       'H_site_sigma_angstrom':SIGMA_A,
       'H_site_epsilon_kcal_per_mol':EPS_KCAL_MOL,
       'H_site_charge_e':CHARGE_E,
       'LJ_cutoff_angstrom':CUTOFF_A,
       'bond_length_angstrom_ASSUMED_not_from_paper': BOND_A,
       'h_mass_g_per_mol_ASSUMED': H_MASS,
       'molecules':int(pos.shape[0]),'atoms':int(pos.shape[0]*2),
       'init_box_angstrom': length,
       'init_box_assumption':'ideal gas PV=nRT, not hydrogen EOS',
       'temperature_K':args.temperature_K,'pressure_MPa':args.pressure_MPa,
       'equilibration_ensemble':'NPT','production_ensemble':'NVT',
       'equilibration_steps':args.eq_steps,'production_steps':args.prod_steps,
       'MD_timestep_fs':TIMESTEP_FS,
       'warning':'Short demo runs cannot establish converged transport properties. Validate against experiments and independent replicas.'
    }
    (dest/'model_config.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    print(f'Generated {len(pos)} two-site H2 molecules ({len(pos)*2} H sites)')
    print(f'Initial ideal-gas estimate: L = {length:.3f} A, at {args.pressure_MPa:g} MPa and {args.temperature_K:g} K')
    print('LJ sigma=2.7200 A, epsilon=0.0199 kcal/mol, rc=12.8 A, q=0 (paper Table 1).')
    print('Run from the output directory: lmp -in in.h2_paper')

if __name__=='__main__': main()
