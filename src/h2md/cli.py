"""Command line interface for demo simulation or real trajectory analysis."""
import argparse
import json
from pathlib import Path

import numpy as np

from .analysis import (diffusion_in_si, fit_self_diffusion,
                       mass_density_g_cm3, mean_squared_displacement,
                       number_density, radial_distribution)
from .io import load_lammps_dump, load_trajectory, save_trajectory
from .plotting import plot_analysis
from .simulation import simulate_h2_centers
from .atomistic import (simulate_atomistic_h2, save_atomistic_trajectory,
                        load_atomistic_trajectory, load_atomistic_lammps_dump)
from .visualization import (write_extended_xyz, plot_configuration,
                            write_interactive_html, plot_bond_diagnostics)


def _analyze(traj, out: Path, fit_start: float = .40, fit_end: float = .85,
             model_label: str | None = None):
    out.mkdir(parents=True, exist_ok=True)
    # Enough lags to show the transition, but avoid poorly sampled tail lags.
    max_lag = min(traj.n_frames - 1, max(6, (traj.n_frames - 1) // 3))
    t, msd = mean_squared_displacement(traj, max_lag=max_lag)
    fit = fit_self_diffusion(t, msd, fit_start, fit_end)
    radii, rdf = radial_distribution(traj, bins=65)
    np.savetxt(out / "msd.csv", np.column_stack((t, msd)), delimiter=",",
               comments="", header=f"time_{traj.time_unit},msd_{traj.length_unit}2")
    np.savetxt(out / "rdf.csv", np.column_stack((radii, rdf)), delimiter=",",
               comments="", header=f"radius_{traj.length_unit},g_r")
    plot_analysis(traj, t, msd, fit, radii, rdf, out)
    results = {
        "model": model_label or ("Illustrative single-site LJ H2 centers (NOT quantitatively "
                  "validated H2)" if traj.length_unit == "sigma" else
                  "H2 centers from user-supplied LAMMPS trajectory"),
        "n_molecules": traj.n_molecules,
        "n_frames": traj.n_frames,
        "box_length": traj.box_length,
        "length_unit": traj.length_unit,
        "time_unit": traj.time_unit,
        "frame_dt": traj.frame_dt,
        "number_density_per_length3": number_density(traj),
        "diffusion_coefficient_length2_per_time": fit.coefficient,
        "diffusion_fit_r_squared": fit.r_squared,
        "diffusion_fit_start_time": fit.start_time,
        "diffusion_fit_end_time": fit.end_time,
        "warning": ("A linear MSD fit alone does not establish a diffusive "
                    "regime or convergence. Inspect the chosen fit window, "
                    "trajectory length, and independent runs."),
    }
    si = diffusion_in_si(fit, traj)
    if si is not None:
        results["diffusion_coefficient_m2_s"] = si
    rho = mass_density_g_cm3(traj)
    if rho is not None:
        results["mass_density_g_cm3"] = rho
    if traj.kinetic_energy is not None:
        tot = traj.kinetic_energy + traj.potential_energy
        results["relative_total_energy_change"] = float(
            (tot[-1] - tot[0]) / max(abs(tot[0]), 1e-12))
    (out / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))
    print(f"\nSaved plots, CSV data, and results.json in: {out.resolve()}")


def _atomistic_report(tr, out, fit_start=.40, fit_end=.85):
    """Preserve atomic and molecular trajectory, plus analysis/3D assets."""
    out.mkdir(parents=True, exist_ok=True)
    save_atomistic_trajectory(tr, out / 'atomic_trajectory.npz')
    save_trajectory(tr.center_trajectory(), out / 'com_trajectory.npz')
    write_extended_xyz(tr, out / 'h2_molecules.xyz', stride=3)
    plot_configuration(tr, out / 'h2_atomistic_3d.png')
    plot_bond_diagnostics(tr, out / 'h2_bonds.png')
    write_interactive_html(tr, out / 'h2_interactive.html', max_frames=55)
    model_label = ("Two explicit H atoms per H2, harmonic H-H bonds and intermolecular "
                   "site LJ; UNVALIDATED teaching potential" if tr.length_unit == "sigma"
                   else "H2 molecular COMs reconstructed from two-atom LAMMPS dump")
    _analyze(tr.center_trajectory(),out,fit_start,fit_end,model_label=model_label)
    bond=tr.bond_lengths()
    rfile=out/'results.json'
    result=json.loads(rfile.read_text())
    result['model'] = ('Two explicit H atoms per H2, harmonic H-H bonds, '
                       'intermolecular site LJ; ILLUSTRATIVE unvalidated '
                       'reduced-unit potential' if tr.length_unit=='sigma'
                       else 'H2 molecular COMs reconstructed from 2-atom LAMMPS trajectory')
    result['n_hydrogen_atoms']=2*tr.n_molecules
    result['mean_bond_length_'+tr.length_unit]=float(bond.mean())
    result['bond_length_std_'+tr.length_unit]=float(bond.std())
    result['atomic_forcefield_warning']=('This demo is not a validated H2 force field. '
                'Classical flexible H2 bonds do not reproduce quantum vibrations; '
                'diffusion values are illustrative. Molecular self-diffusion is '
                'calculated from UNWRAPPED H2 centers of mass, not H atom positions.')
    rfile.write_text(json.dumps(result,indent=2)+'\n')
    print('\nAtomistic visualizations: '+str((out/'h2_interactive.html').resolve()))


def main(argv=None):
    parser = argparse.ArgumentParser(prog="h2md", description="H2-only MD tutorial & analysis")
    commands = parser.add_subparsers(dest="command", required=True)
    demo = commands.add_parser("demo", help="run illustrative H2-center LJ simulation + analyze")
    demo.add_argument("--output", default="results", help="output directory")
    demo.add_argument("--steps", type=int, default=3000, help="NVE production steps")
    demo.add_argument("--equil", type=int, default=500, help="thermostatted equilibration steps")
    demo.add_argument("--seed", type=int, default=42)
    demo.add_argument("--fit-start", type=float, default=.40)
    demo.add_argument("--fit-end", type=float, default=.85)
    atom = commands.add_parser('atomistic-demo',
        help='explicit 2-H-atom flexible H2 molecular MD + 3D animation')
    atom.add_argument('--output', default='results/atomistic')
    atom.add_argument('--steps',type=int,default=2500)
    atom.add_argument('--equil',type=int,default=400)
    atom.add_argument('--seed',type=int,default=42)
    atom.add_argument('--fit-start',type=float,default=.40)
    atom.add_argument('--fit-end',type=float,default=.85)
    adump = commands.add_parser('atomistic-analyze',
        help='analyze H2 2-atom/molecule LAMMPS dump + 3D animation')
    adump.add_argument('--input',required=True)
    adump.add_argument('--output',default='results/atomistic-lammps')
    adump.add_argument('--timestep-fs',type=float,
                        help='LAMMPS units real MD timestep (fs); required for dump')
    adump.add_argument('--fit-start',type=float,default=.40)
    adump.add_argument('--fit-end',type=float,default=.85)
    analyze = commands.add_parser("analyze", help="analyze .npz or LAMMPS custom .lammpstrj")
    analyze.add_argument("--input", required=True)
    analyze.add_argument("--output", default="results")
    analyze.add_argument("--timestep-fs", type=float, default=None,
                         help="LAMMPS units real timestep in femtoseconds")
    analyze.add_argument("--fit-start", type=float, default=.40)
    analyze.add_argument("--fit-end", type=float, default=.85)
    args = parser.parse_args(argv)
    if args.command == "demo":
        print("Simulating illustrative single-site H2 centers in reduced LJ units...")
        traj = simulate_h2_centers(production_steps=args.steps,
                                   equilibration_steps=args.equil, seed=args.seed)
        output = Path(args.output)
        save_trajectory(traj, output / "trajectory.npz")
        _analyze(traj, output, args.fit_start, args.fit_end)
    elif args.command == 'atomistic-demo':
        print('Running two-atom H2 illustrative flexible dimer MD (reduced units)...')
        atraj=simulate_atomistic_h2(equilibration_steps=args.equil,
                                   production_steps=args.steps, seed=args.seed)
        _atomistic_report(atraj,Path(args.output),args.fit_start,args.fit_end)
    elif args.command == 'atomistic-analyze':
        p=Path(args.input)
        if p.suffix.lower()=='.npz':
            atraj=load_atomistic_trajectory(p)
        else:
            if args.timestep_fs is None: parser.error('--timestep-fs is required for LAMMPS dump')
            atraj=load_atomistic_lammps_dump(p,timestep_fs=args.timestep_fs)
        _atomistic_report(atraj,Path(args.output),args.fit_start,args.fit_end)
    else:
        path = Path(args.input)
        if path.suffix.lower() == ".npz":
            traj = load_trajectory(path)
        else:
            if args.timestep_fs is None:
                parser.error("--timestep-fs is required for LAMMPS dump input")
            traj = load_lammps_dump(path, timestep_fs=args.timestep_fs)
        _analyze(traj, Path(args.output), args.fit_start, args.fit_end)


if __name__ == "__main__":
    main()
