"""Analyze independent H2 MD runs: NPT density vs CoolProp + NVT diffusion CI.

Install: python -m pip install CoolProp scipy matplotlib numpy
Run: python analyze_replicas.py --runs runs_v05 --out analysis_v05
No demonstration data are fabricated: requires completed LAMMPS trajectories.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import numpy as np

from h2md.atomistic import load_atomistic_lammps_dump

MPA_TO_ATM = 9.869232667160128


def density_window(path: Path, keep_fraction: float = .5):
    """LAMMPS fix ave/time series: step, rho(g/cc), T(K), P(atm), volume(A^3)."""
    if not 0 < keep_fraction <= 1: raise ValueError('keep fraction must be within (0,1]')
    data=np.loadtxt(path,comments='#',ndmin=2)
    if data.shape[1] != 5 or data.shape[0] < 6 or not np.isfinite(data).all():
        raise ValueError('NPT thermo file needs >=6 valid rows and five columns')
    # Pressure is an *instantaneous, signed* virial observable. Individual
    # samples can be zero or negative even when fix npt targets positive P.
    # In contrast, step order, mass density, temperature, and box volume have
    # strict physical constraints. Never discard negative-P samples or take
    # their absolute values: that would bias the reported mean pressure.
    if not np.all(np.diff(data[:,0]) > 0):
        raise ValueError('NPT timesteps must be strictly increasing')
    if np.any(data[:,1] <= 0) or np.any(data[:,2] <= 0) or np.any(data[:,4] <= 0):
        raise ValueError('invalid NPT samples: density, temperature and volume must be > 0')
    n_keep = max(3,int(np.ceil(keep_fraction*len(data))))
    used=data[-n_keep:]
    return {
        'n_total':len(data), 'n_used':len(used),
        'density_kg_m3':float(np.mean(used[:,1])*1000.),
        'temperature_K':float(np.mean(used[:,2])),
        'pressure_MPa':float(np.mean(used[:,3])/MPA_TO_ATM),
        'negative_pressure_samples_used':int(np.count_nonzero(used[:,3] < 0)),
        'density_first_half_kg_m3':float(np.mean(used[:len(used)//2,1])*1000.),
        'density_second_half_kg_m3':float(np.mean(used[len(used)//2:,1])*1000.),
        'volume_A3':float(np.mean(used[:,4])),
        'raw_steps':data[:,0], 'raw_density_kg_m3':data[:,1]*1000.,
    }


def multi_origin_msd(com: np.ndarray, frame_dt_fs: float, lag_fraction: float=.5,
                     n_lags: int=151, origin_stride:int=10, remove_drift:bool=True):
    """Return lag times in ps and molecular COM MSD in A^2.

    Overlapping origins improve the estimator but are *not* independent
    replications. Unwrapped center-of-mass coordinates are required.
    """
    xyz=np.asarray(com,dtype=float)
    if xyz.ndim!=3 or xyz.shape[2]!=3 or xyz.shape[0]<12 or xyz.shape[1]<2:
        raise ValueError('COM data must have >=12 frames, >=2 molecules, xyz')
    if not np.isfinite(xyz).all() or frame_dt_fs<=0 or origin_stride<1:
        raise ValueError('invalid coordinates/time/origin stride')
    if not .1 <= lag_fraction <= .8: raise ValueError('lag fraction must be .1-.8')
    maxlag=int((len(xyz)-1)*lag_fraction)
    lags=np.unique(np.rint(np.linspace(1,maxlag,min(maxlag,n_lags))).astype(int))
    times_ps=lags*frame_dt_fs/1000.0
    values=np.empty(len(lags),dtype=float)
    for k, lag in enumerate(lags):
        start=np.arange(0,len(xyz)-lag,origin_stride)
        delta=xyz[start+lag]-xyz[start]
        if remove_drift:
            # correct per-origin framewise translational drift across molecules
            delta=delta-delta.mean(axis=1,keepdims=True)
        values[k]=np.mean(np.sum(delta*delta,axis=-1))
    return times_ps, values


def fit_diffusion(times_ps, msd_A2, fit_start_fraction=.10,
                  fit_end_fraction=.50, duration_ps=None):
    """Fit MSD = 6Dt + c. D units A^2/ps then converted to m^2/s."""
    times=np.asarray(times_ps,dtype=float)
    msd=np.asarray(msd_A2,dtype=float)
    if len(times)!=len(msd) or len(times)<6 or not np.isfinite(msd).all():
        raise ValueError('invalid MSD values')
    if not 0 < fit_start_fraction < fit_end_fraction <= .8:
        raise ValueError('fit bounds must satisfy 0 < start < end <= .8')
    duration_ps=(float(duration_ps) if duration_ps is not None else float(times[-1])/.5)
    a,b=duration_ps*fit_start_fraction,duration_ps*fit_end_fraction
    mask=(times>=a)&(times<=b)
    if mask.sum()<5:
        raise ValueError('not enough points inside chosen fit window')
    slope,intercept=np.polyfit(times[mask],msd[mask],1)
    if slope<=0: raise ValueError('nonpositive MSD slope; diffusion cannot be reported')
    predicted=slope*times[mask]+intercept
    ss_res=float(np.sum((msd[mask]-predicted)**2))
    ss_tot=float(np.sum((msd[mask]-msd[mask].mean())**2))
    r2=1-ss_res/ss_tot if ss_tot>0 else float('nan')
    return {'D_m2_s':float(slope/6.*1e-8), 'slope_A2_ps':float(slope),
            'intercept_A2':float(intercept),'r2':float(r2),
            'fit_start_ps':float(a),'fit_end_ps':float(b),
            'fit_points':int(np.sum(mask))}


def mean_ci95(values):
    """95% Student-t CI across independent seeds, NOT across MD frames."""
    from scipy.stats import t
    x=np.asarray(values,dtype=float)
    if x.ndim!=1 or len(x)<2 or not np.isfinite(x).all():
        raise ValueError('need at least two finite independent run estimates')
    m=float(x.mean()); s=float(x.std(ddof=1)); se=s/np.sqrt(len(x))
    half=float(t.ppf(.975,len(x)-1)*se)
    return {'mean':m,'sd':s,'n':int(len(x)),'ci95_low':m-half,'ci95_high':m+half,
            'ci95_half_width':half}


def reference_hydrogen_density(temp_k,pressure_mpa):
    try:
        import CoolProp
        from CoolProp.CoolProp import PropsSI
    except ImportError as e:
        raise RuntimeError('CoolProp is REQUIRED for the EOS benchmark: python -m pip install CoolProp') from e
    # CoolProp HEOS pure-normal hydrogen; mass density in kg/m^3
    return float(PropsSI('Dmass','T',float(temp_k),'P',float(pressure_mpa)*1e6,'HEOS::Hydrogen')), CoolProp.__version__


def write_csv(path, columns, rows):
    with Path(path).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(rows)




def dump_coordinate_precision(input_path: Path) -> int:
    """Significant digits for LAMMPS custom dump xu/yu/zu coordinates.

    LAMMPS default for custom dumps is printf %g (six significant digits).
    Detects the explicit %.Ng formatting used in upgraded simulations.
    This only covers the exact generator's dump format, not arbitrary imports.
    """
    content = Path(input_path).read_text(encoding='utf-8')
    match = re.search(r'(?m)^\s*dump_modify\s+TRAJ\s+format\s+float\s+%\.?([0-9]+)g\s*$', content)
    if match:
        return int(match.group(1))
    # For an unmodified legacy v0.5 dump, LAMMPS default is %g.
    if re.search(r'(?m)^\s*dump\s+TRAJ\s+all\s+custom\s+', content):
        return 6
    raise ValueError(f'Cannot identify dump precision from {input_path}; review atom coordinate format')


def bond_precision_audit(atoms, reference_A: float, significant_digits: int,
                         constraint_slack_A: float = 0.0001) -> dict:
    """Check bond residuals against *per-coordinate text-rounding bounds*.

    A printed position with N significant digits has absolute rounding error
    <= 0.5*10**(floor(log10(abs(x)))-N+1) per axis (ordinary printf %g).
    The difference of two rounded positions has component error <= sum of
    individual errors. Norm is 1-Lipschitz; therefore each saved bond-length
    deviation is bounded by the L2 norm of those three summed errors plus
    constraint_slack_A for small MD SHAKE residuals.

    Passing means CONSISTENT with dump rounding, NOT proof SHAKE was exact.
    """
    coords = np.asarray(atoms, dtype=float)
    if coords.ndim != 4 or coords.shape[-2:] != (2, 3):
        raise ValueError('Expected (frames,molecules,2,3) coordinates')
    if not np.isfinite(coords).all() or not 3 <= significant_digits <= 17:
        raise ValueError('Invalid coordinates or dump precision')
    if reference_A <= 0 or constraint_slack_A < 0:
        raise ValueError('Invalid bond reference or constraint slack')

    lengths = np.linalg.norm(coords[:, :, 1, :] - coords[:, :, 0, :], axis=-1)
    deviation = np.abs(lengths - reference_A)
    magnitude = np.abs(coords)
    halfstep = np.zeros_like(magnitude)
    valid = magnitude > 0
    exponent = np.floor(np.log10(magnitude[valid]))
    halfstep[valid] = 0.5 * np.power(10.0, exponent - significant_digits + 1)
    bound = np.linalg.norm(halfstep[:, :, 0, :] + halfstep[:, :, 1, :], axis=-1)
    exceed = deviation - bound - constraint_slack_A
    n_outliers = int(np.count_nonzero(exceed > 1.e-12))
    return {
        'n_bonds_observed': int(deviation.size),
        'max_deviation_A': float(np.max(deviation)),
        'p99_deviation_A': float(np.percentile(deviation, 99)),
        'max_rounding_bound_A': float(np.max(bound)),
        'n_deviations_exceeding_rounding_bound': n_outliers,
        'max_excess_beyond_rounding_bound_A': float(max(0., np.max(exceed))),
        'dump_significant_digits': significant_digits,
        'consistent_with_dump_precision': n_outliers == 0,
    }


def evaluate_runs(runs:Path, out:Path, density_keep:float=.50,
                  fit_start_frac:float=.10,fit_end_frac:float=.50,
                  origin_stride:int=10, plot:bool=True):
    runs,out=Path(runs),Path(out)
    cfg=json.loads((runs/'experiment.json').read_text(encoding='utf-8'))
    out.mkdir(parents=True,exist_ok=True)
    nominal_t,nominal_p=cfg['temperature_K'],cfg['pressure_MPa']
    ref_nominal,cp_version=reference_hydrogen_density(nominal_t,nominal_p)
    rows=[]; curves=[]; thermo=[]; warnings=[]
    for seed in cfg['seeds']:
        r=runs/f'seed_{seed}'
        tp=r/'output'/'npt_thermo.dat'
        dp=r/'output'/'h2_production.lammpstrj'
        if not tp.exists() or not dp.exists():
            raise FileNotFoundError(f'Missing {tp if not tp.exists() else dp}. Run LAMMPS first.')
        d=density_window(tp,density_keep)
        tr=load_atomistic_lammps_dump(dp, timestep_fs=cfg['timestep_fs'])
        if tr.n_molecules!=cfg['molecules']:
            raise ValueError(f'seed {seed}: expected {cfg["molecules"]} H2; got {tr.n_molecules}')
        if tr.box_length <= 2*12.8:
            raise ValueError(f'seed {seed}: box length {tr.box_length:.3f} A <= 2*cutoff')
        bond_A = float(cfg.get('assumptions', {}).get('H_H_bond_A', .7414))
        n_digits = dump_coordinate_precision(r/'in.h2_validation')
        bond_check = bond_precision_audit(tr.atoms, bond_A, n_digits)
        maxbond = bond_check['max_deviation_A']
        if not bond_check['consistent_with_dump_precision']:
            raise ValueError(
                f'seed {seed}: apparent H-H deviation {maxbond:.4f} A exceeds '
                f'the {n_digits}-significant-digit dump rounding bound '
                f'for {bond_check["n_deviations_exceeding_rounding_bound"]} '
                'bond observations. Inspect trajectory, molecule mapping, and SHAKE log.'
            )
        if maxbond > .001:
            warnings.append(
                f'seed {seed}: apparent H-H deviation {maxbond:.4f} A is '
                f'consistent with {n_digits}-significant-digit trajectory rounding '
                '(not independent verification of SHAKE). Use 16-digit dump '
                'coordinates for future runs; retain this limitation in interpretation.'
            )
        # First and last saved frames span production but LAMMPS does not dump step 0
        duration_ps=(tr.n_frames-1)*tr.frame_dt/1000
        times,msd=multi_origin_msd(tr.com,tr.frame_dt,
               origin_stride=origin_stride)
        fit=fit_diffusion(times,msd,fit_start_frac,fit_end_frac,duration_ps)
        # Sensitivity of D to the chosen time window; heuristic, not a proof of convergence.
        alt_start=max(.15,fit_start_frac)
        alt_end=min(.40,fit_end_frac)
        alt = (fit_diffusion(times,msd,alt_start,alt_end,duration_ps)
               if alt_end-alt_start >= .15 else None)
        fit_window_change_pct=(100*abs(alt['D_m2_s']/fit['D_m2_s']-1) if alt else None)
        # Reference at measured mean T,P is meaningful only for a positive
        # absolute pressure. A very short NPT run can have a nonpositive
        # retained mean, despite its positive barostat setpoint. Still report
        # the target-state benchmark, and flag the run for more equilibration.
        if d['pressure_MPa'] > 0:
            ref_at_measured,_=reference_hydrogen_density(d['temperature_K'],d['pressure_MPa'])
            measured_density_error=100*(d['density_kg_m3']/ref_at_measured-1)
        else:
            ref_at_measured=None
            measured_density_error=None
            warnings.append(
                f'seed {seed}: retained mean NPT pressure is {d["pressure_MPa"]:.3g} MPa '
                '(nonpositive); measured-state EOS comparison omitted. '
                'Extend NPT equilibration and inspect pressure convergence'
            )
        row={
           'seed':seed,'production_ps_recorded':round(duration_ps,3),
           'npt_samples_used':d['n_used'],'temperature_K':d['temperature_K'],
           'pressure_MPa':d['pressure_MPa'],
           'md_density_kg_m3':d['density_kg_m3'],
           'ref_density_measured_TP_kg_m3':ref_at_measured,
           'density_error_percent_at_measured_TP':measured_density_error,
           'negative_pressure_samples_used':d['negative_pressure_samples_used'],
           'D_m2_s':fit['D_m2_s'],'msd_R2':fit['r2'],
           'D_fit_start_ps':fit['fit_start_ps'],'D_fit_end_ps':fit['fit_end_ps'],
           'D_fit_window_sensitivity_percent':fit_window_change_pct,
           'bond_max_deviation_A':maxbond,'bond_dump_precision_digits':n_digits,
           'bond_max_rounding_bound_A':bond_check['max_rounding_bound_A'],
           'bond_rounding_outliers':bond_check['n_deviations_exceeding_rounding_bound'],
           'nvt_box_A':tr.box_length,
           }
        rows.append(row)
        curves.append((seed,times,msd,fit))
        thermo.append((seed,d))
        # A basic plateau screen, not a rigorous test of equilibrium.
        split=abs(d['density_second_half_kg_m3']-d['density_first_half_kg_m3'])/d['density_kg_m3']
        if split>.03:
            warnings.append(f'seed {seed}: retained NPT density half-means differ by {100*split:.1f}%: extend NPT and inspect the trace')
        if abs(d['pressure_MPa']/nominal_p-1)>.10:
            warnings.append(f'seed {seed}: NPT mean pressure differs >10% from target; inspect equilibration')
        if abs(d['temperature_K']/nominal_t-1)>.03:
            warnings.append(f'seed {seed}: NPT mean temperature differs >3% from target; inspect equilibration')
        if fit['r2'] < .98:
            warnings.append(f'seed {seed}: MSD fit R²={fit["r2"]:.3f}; inspect diffusive window')
        if fit_window_change_pct is not None and fit_window_change_pct > 20:
            warnings.append(f'seed {seed}: diffusion differs by {fit_window_change_pct:.1f}% in an alternative fit window; extend production')
    dstat=mean_ci95([r['md_density_kg_m3'] for r in rows])
    Dstat=mean_ci95([r['D_m2_s'] for r in rows])
    if cfg['preset']=='smoke':
        warnings.append('Smoke mode is only a software check; 50 ps production is too short to establish reliable diffusion')
    if len(rows)<3: warnings.append('At least 3 independent replicas recommended; 5+ preferable for precision')
    if Dstat['mean']>0 and Dstat['ci95_half_width']/Dstat['mean']>.30:
        warnings.append('Relative diffusion CI half-width >30%; add replicas and/or extend production')
    aggregate={
      'source_doi':cfg['source_doi'],'reference':'CoolProp HEOS::Hydrogen (normal H2; Leachman EOS)',
      'coolprop_version':cp_version,'temperature_K_target':nominal_t,
      'pressure_MPa_target':nominal_p,'density_reference_target_kg_m3':ref_nominal,
      'density_kg_m3':dstat, 'diffusion_m2_s':Dstat,
      'density_relative_error_percent_at_target':100*(dstat['mean']/ref_nominal-1),
      'density_keep_fraction':density_keep,
      'MSD_fit_fraction_of_recorded_production':[fit_start_frac,fit_end_frac],
      'MSD_origin_stride_frames':origin_stride,'replicas':rows,'warnings':warnings,
      'bond_precision_method':'Output-precision-aware roundoff envelope; passing does not independently prove SHAKE constraint accuracy.',
      'interpretation':'95% CI uses independent-seed estimates, not autocorrelated frames or overlapping time origins. Results must be checked for equilibration, finite-size, fit-window, and force-field errors.'
    }
    (out/'summary.json').write_text(json.dumps(aggregate,indent=2)+'\n',encoding='utf-8')
    write_csv(out/'replicas.csv',list(rows[0].keys()),rows)
    with (out/'msd_curves.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['seed','lag_ps','msd_A2'])
        for seed,times,msd,fit in curves:
            w.writerows((seed,float(t),float(m)) for t,m in zip(times,msd))
    densityline=f"{dstat['mean']:.5g} [{dstat['ci95_low']:.5g}, {dstat['ci95_high']:.5g}] kg/m³"
    dline=f"{Dstat['mean']:.5g} [{Dstat['ci95_low']:.5g}, {Dstat['ci95_high']:.5g}] m²/s"
    lines=[
        '# H₂ MD validation report — Version 0.5','',
        f"**Conditions:** {nominal_t:g} K, {nominal_p:g} MPa, {cfg['molecules']} H₂ molecules, {len(rows)} independent seeds.",
        f"**Run lengths:** {cfg['equilibration_steps']*cfg['timestep_fs']/1000:g} ps NPT + {cfg['production_steps']*cfg['timestep_fs']/1000:g} ps NVT per seed.",
        f"**Normal-H₂ reference density (CoolProp/Leachman EOS):** {ref_nominal:.5g} kg/m³.",
        f"**NPT density mean (95% CI between seeds):** {densityline}.",
        f"**Density deviation from EOS at target state:** {aggregate['density_relative_error_percent_at_target']:+.2f}%.",
        f"**H₂ molecular COM self-diffusion mean (95% CI between seeds):** {dline}.",
        f"**MSD fit interval:** {fit_start_frac:.0%}–{fit_end_frac:.0%} of recorded production duration, with overlapping time origins.",
        '',
        '**Important:** Density is an EOS benchmark, not experimental validation of diffusion. The t-based interval quantifies between-seed statistical variation only; it does not include force-field systematic errors or finite-size and fit-window biases. Negative instantaneous NPT pressures are retained; only nonpositive *mean* pressure prevents a measured-state EOS comparison.','',
        '**Bond-quality note:** The audit compares bond length deviations to a conservative estimate of the coordinate text-rounding envelope. If the text uses only six significant digits, the saved coordinates cannot prove SHAKE tolerance. See per-seed `bond_dump_precision_digits`, `bond_max_deviation_A`, and `bond_rounding_outliers` in `replicas.csv`.','',
        '## Checks and warnings','',
    ]
    lines += [f'- {w}' for w in warnings] if warnings else ['- No automated screening warnings. This is not a convergence guarantee.']
    lines += ['', '## Model provenance','',
              '- LJ site parameters: Adam et al. (2024), Table 1, DOI:10.1016/j.ijhydene.2024.06.319.',
              '- H–H constraint length 0.7414 Å is a separate assumption, not reported in Table 1.',
              '- Reference equation of state: Leachman et al. (2009), via CoolProp `HEOS::Hydrogen`.',
              '- Normal hydrogen (EOS) differs from the classical rigid dimer model in treatment of molecular quantum effects.','']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    if plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(8,4.8))
        for seed,d in thermo:
            ax.plot(d['raw_steps']*cfg['timestep_fs']/1000, d['raw_density_kg_m3'],
                    lw=.9,alpha=.8,label=f'Seed {seed}')
        ax.axhline(ref_nominal,c='black',lw=1.5,ls='--',label='CoolProp EOS (target T, P)')
        ax.set(xlabel='NPT time (ps)',ylabel='Density (kg/m³)',title='H₂ density equilibration and EOS benchmark')
        ax.legend();fig.tight_layout();fig.savefig(out/'density_validation.png',dpi=180);plt.close(fig)
        fig,ax=plt.subplots(figsize=(8,4.8))
        for seed,times,msd,fit in curves:
            line=ax.plot(times,msd,lw=1.6,label=f'Seed {seed}')[0]
            x=np.array([fit['fit_start_ps'],fit['fit_end_ps']])
            y=fit['slope_A2_ps']*x+fit['intercept_A2']
            ax.plot(x,y,ls='--',c=line.get_color(),lw=2)
        ax.set(xlabel='Lag time (ps)',ylabel='COM MSD (Å²)',title='Multi-origin molecular MSD and linear fit')
        ax.legend();fig.tight_layout();fig.savefig(out/'diffusion_msd.png',dpi=180);plt.close(fig)
        fig,axes=plt.subplots(1,2,figsize=(10,4))
        seeds=[str(r['seed']) for r in rows]
        axes[0].plot(seeds,[r['md_density_kg_m3'] for r in rows],'o',ms=9)
        axes[0].axhline(ref_nominal,c='black',ls='--',lw=1.4,label='EOS')
        axes[0].set(xlabel='Independent seed',ylabel='Density (kg/m³)',title='Density by run')
        axes[0].legend()
        axes[1].plot(seeds,[r['D_m2_s'] for r in rows],'o',ms=9)
        axes[1].axhline(Dstat['mean'],c='black',ls='--',lw=1.4,label='Replica mean')
        axes[1].set(xlabel='Independent seed',ylabel='D (m²/s)',title='COM diffusion by run')
        axes[1].legend();fig.tight_layout();fig.savefig(out/'replica_comparison.png',dpi=180);plt.close(fig)
    return aggregate


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',default='runs_v05')
    p.add_argument('--out',default='analysis_v05')
    p.add_argument('--density-keep-fraction',type=float,default=.5)
    p.add_argument('--fit-start-frac',type=float,default=.10)
    p.add_argument('--fit-end-frac',type=float,default=.50)
    p.add_argument('--origin-stride',type=int,default=10)
    p.add_argument('--no-plots',action='store_true')
    a=p.parse_args(argv)
    result=evaluate_runs(Path(a.runs),Path(a.out),a.density_keep_fraction,
             a.fit_start_frac,a.fit_end_frac,a.origin_stride,not a.no_plots)
    print(json.dumps({k:result[k] for k in ['density_reference_target_kg_m3','density_kg_m3',
          'density_relative_error_percent_at_target','diffusion_m2_s','warnings']},indent=2))
    print(f'Wrote {a.out}/REPORT.md and figures.')

if __name__=='__main__':main()
