"""End-to-end orchestration and comparative reporting for H2 LAMMPS campaigns.

Designed as a thin, testable layer over the existing v0.5 engine. It does not
change force-field parameters or silently overwrite existing simulation data.
"""
from __future__ import annotations

import argparse
import csv
import html
import importlib.util
import json
import re
import subprocess
from pathlib import Path

SCHEMA_VERSION = '1.0'


def _engine_file(name: str) -> Path:
    """Locate the existing engine in an editable repository checkout."""
    path = Path(__file__).resolve().parents[2] / 'lammps' / 'paper_h2' / 'validation_v05' / name
    if not path.is_file():
        raise FileNotFoundError(f'Missing {path}; install from the complete GitHub checkout, not the wheel.')
    return path


def _load_engine(name: str):
    file = _engine_file(name)
    spec = importlib.util.spec_from_file_location(f'h2md_workflow_{file.stem}', file)
    if spec is None or spec.loader is None:
        raise RuntimeError(f'Cannot import {file}')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _read_json(path: Path) -> dict:
    data = json.loads(path.read_text(encoding='utf-8'))
    if not isinstance(data, dict):
        raise ValueError(f'{path} must hold a JSON object')
    return data


def _safe_id(value: str) -> str:
    if not isinstance(value, str) or not re.fullmatch(r'[a-zA-Z0-9][a-zA-Z0-9_-]{0,63}', value):
        raise ValueError('Experiment ids must be simple letters, digits, dash, underscore (max 64 chars)')
    return value


def plan(config_path: Path, campaign: Path) -> list[dict]:
    """Create all LAMMPS inputs, never run LAMMPS in this step."""
    config = _read_json(Path(config_path))
    defaults = config.get('defaults', {})
    experiments = config.get('experiments')
    if not isinstance(defaults, dict) or not isinstance(experiments, list) or not experiments:
        raise ValueError('Expected defaults object and non-empty experiments array')
    campaign = Path(campaign)
    if (campaign / 'campaign.json').exists():
        raise FileExistsError(f'Campaign already exists: {campaign} (will not overwrite)')
    if campaign.exists() and any(campaign.iterdir()):
        raise FileExistsError(f'Campaign directory not empty: {campaign}')
    engine = _load_engine('run_replicas.py')
    rows = []
    ids = set()
    for exp in experiments:
        if not isinstance(exp, dict):
            raise ValueError('Each experiment must be an object')
        opts = defaults | exp
        eid = _safe_id(exp.get('id'))
        if eid in ids:
            raise ValueError(f'Duplicate experiment ID: {eid}')
        ids.add(eid)
        cells = int(opts.get('cells', 4))
        seeds = opts.get('seeds', [2027, 2028, 2029])
        if not isinstance(seeds, list) or not seeds or any(type(seed) is not int for seed in seeds):
            raise ValueError('seeds must be a nonempty list of integer seeds')
        preset = opts.get('preset', 'standard')
        if preset not in engine.PRESETS_PS:
            raise ValueError(f'Unsupported preset: {preset}')
        rows.append(dict(id=eid, cells=cells, seeds=seeds, preset=preset,
                         temperature_K=float(opts.get('temperature_K', 323.15)),
                         pressure_MPa=float(opts.get('pressure_MPa', 10.0)),
                         equilibration_ps=opts.get('equilibration_ps'),
                         production_ps=opts.get('production_ps'),
                         thermo_every=int(opts.get('thermo_every', 500)),
                         dump_every=int(opts.get('dump_every', 1000))))
    campaign.mkdir(parents=True, exist_ok=True)
    for row in rows:
        run_dir = campaign / row['id'] / 'runs'
        engine.prepare(run_dir, row['preset'], row['seeds'], row['cells'],
                       row['temperature_K'], row['pressure_MPa'],
                       row['equilibration_ps'], row['production_ps'],
                       row['thermo_every'], row['dump_every'])
    manifest = {'version': SCHEMA_VERSION, 'source_config': str(Path(config_path).resolve()),
                'experiments': rows}
    (campaign / 'campaign.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    return rows


def _manifest(campaign: Path) -> dict:
    return _read_json(Path(campaign) / 'campaign.json')


def execute(campaign: Path, lmp: str = 'lmp', force: bool = False) -> list[str]:
    """Run sequentially; do not replace successful trajectories unless --force is explicit."""
    campaign = Path(campaign)
    manifest = _manifest(campaign)
    done = []
    for exp in manifest['experiments']:
        eid = _safe_id(exp['id'])
        for seed in exp['seeds']:
            root = campaign / eid / 'runs' / f'seed_{seed}'
            dump = root / 'output' / 'h2_production.lammpstrj'
            if dump.is_file() and dump.stat().st_size > 0 and not force:
                done.append(f'SKIPPED existing {eid}/seed_{seed}')
                continue
            with (root / 'lammps.log').open('w', encoding='utf-8') as log:
                subprocess.run([lmp, '-in', 'in.h2_validation'], cwd=root,
                               stdout=log, stderr=subprocess.STDOUT, check=True)
            if not dump.is_file() or dump.stat().st_size == 0:
                raise RuntimeError(f'LAMMPS returned successfully but trajectory missing: {dump}')
            done.append(f'COMPLETED {eid}/seed_{seed}')
    return done


def analyze_campaign(campaign: Path, plot: bool = True) -> list[Path]:
    campaign = Path(campaign)
    engine = _load_engine('analyze_replicas.py')
    outputs = []
    for exp in _manifest(campaign)['experiments']:
        eid = _safe_id(exp['id'])
        runs = campaign / eid / 'runs'
        out = campaign / eid / 'analysis'
        engine.evaluate_runs(runs, out, plot=plot)
        outputs.append(out)
    return outputs


def _stats(summary: dict, key: str) -> tuple[float, float, int]:
    item = summary[key]
    for attr in ('mean', 'ci95_half_width', 'n'):
        if attr not in item:
            raise ValueError(f'Missing {key}.{attr}')
    mean = float(item['mean'])
    half = float(item['ci95_half_width'])
    n = int(item['n'])
    if not (mean > 0 and half >= 0 and n >= 2):
        raise ValueError(f'Invalid summary stats for {key}')
    return mean, half, n


def collect(analyses: list[Path]) -> list[dict]:
    """Collect existing v0.5 analyses; no raw trajectory or EOS dependency required."""
    if not analyses:
        raise ValueError('Pass at least one analysis folder containing summary.json')
    rows = []
    for item in analyses:
        source = Path(item)
        path = source / 'summary.json' if source.is_dir() else source
        if not path.is_file():
            raise FileNotFoundError(f'Analysis summary missing: {path}')
        d = _read_json(path)
        density, drho, n = _stats(d, 'density_kg_m3')
        diffusion, dD, dn = _stats(d, 'diffusion_m2_s')
        molecules = int(d.get('molecules', d.get('n_molecules', 0)) or 0)
        if not molecules:
            # v0.5 has molecule count in experiment.json beside analysis
            for candidate in (path.parent.parent / 'runs' / 'experiment.json',
                              path.parent.parent / 'experiment.json',
                              path.parent.parent / f'runs_{path.parent.name.removeprefix("analysis_")}' / 'experiment.json'):
                if candidate.exists():
                    molecules = int(_read_json(candidate)['molecules'])
                    break
        if not molecules:
            # The v0.5 summary has replicate-level metadata but omits N;
            # the input directory's semantic name may identify the run size.
            labels = {'analysis_research': 256, 'analysis_500': 500, 'analysis_864': 864}
            molecules = labels.get(path.parent.name, 0)
        if molecules <= 0:
            raise ValueError(f'Cannot establish molecule count for {path}. Place beside runs/experiment.json.')
        target_rho = float(d['density_reference_target_kg_m3'])
        T = float(d.get('temperature_K_target', 0))
        P = float(d.get('pressure_MPa_target', 0))
        row = {'name': path.parent.name, 'molecules': molecules, 'replicas': n,
               'temperature_K': T, 'pressure_MPa': P,
               'density_kg_m3': density, 'density_ci95_half_width': drho,
               'reference_density_kg_m3': target_rho,
               'density_error_percent': float(d['density_relative_error_percent_at_target']),
               'diffusion_m2_s': diffusion, 'diffusion_ci95_half_width': dD,
               'diffusion_1e6_m2_s': diffusion * 1e6,
               'diffusion_ci95_1e6_m2_s': dD * 1e6,
               'n_warnings': len(d.get('warnings', [])),
               'source_file': str(path)}
        if dn != n:
            raise ValueError(f'Density and diffusion replica count disagree: {path}')
        rows.append(row)
    return sorted(rows, key=lambda row: (row['temperature_K'], row['pressure_MPa'], row['molecules']))


def _write_markdown(rows: list[dict], target: Path):
    lines = ['# H₂ molecular dynamics — campaign overview', '',
             'This is a comparative *simulation* summary, not an experimental validation of self-diffusion.',
             'Confidence intervals represent **between-replica statistical uncertainty** only.', '',
             '| H₂ molecules | T (K) | P (MPa) | Replicas | Density (kg/m³) | EOS error | D (10⁻⁶ m²/s) | Warnings |',
             '|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        lines.append(f"| {r['molecules']} | {r['temperature_K']:g} | {r['pressure_MPa']:g} | {r['replicas']} | "
                     f"{r['density_kg_m3']:.4f} ± {r['density_ci95_half_width']:.4f} | "
                     f"{r['density_error_percent']:+.3f}% | "
                     f"{r['diffusion_1e6_m2_s']:.4f} ± {r['diffusion_ci95_1e6_m2_s']:.4f} | {r['n_warnings']} |")
    lines += ['', '## Reproducibility and interpretation', '',
              '- Density is from the equilibrated portion of the NPT stage; H₂ self-diffusion is from molecular COM positions during NVT production.',
              '- Diffusion CI does not include force-field bias, finite-size corrections, or fitting-window systematic error.',
              '- Compare sizes only at the same thermodynamic state and protocol; overlapping CIs do not prove size independence.',
              '- Model: two-site LJ site parameters from Adam et al., *International Journal of Hydrogen Energy* (2024), DOI:10.1016/j.ijhydene.2024.06.319. Rigid H–H bond is a separate implementation assumption.',
              '- Source files and exact per-run diagnostics are in `summary.json`, `replicas.csv` and each analysis `REPORT.md`.', '']
    target.write_text('\n'.join(lines), encoding='utf-8')


def _create_dashboard(rows: list[dict], output: Path):
    import plotly.graph_objects as go
    import plotly.io as pio
    groups = sorted(set((r['temperature_K'],r['pressure_MPa']) for r in rows))
    fD = go.Figure(); fRho = go.Figure()
    for T, P in groups:
        subset = [r for r in rows if r['temperature_K']==T and r['pressure_MPa']==P]
        legend = f'{T:g} K · {P:g} MPa'
        hover = [f"{r['molecules']} molecules · {r['replicas']} replicas" for r in subset]
        fD.add_trace(go.Scatter(x=[r['molecules'] for r in subset], y=[r['diffusion_1e6_m2_s'] for r in subset],
                                error_y=dict(type='data',array=[r['diffusion_ci95_1e6_m2_s'] for r in subset], visible=True),
                                name=legend, text=hover, hovertemplate='%{text}<br>D=%{y:.4f} ×10⁻⁶ m²/s<extra></extra>',
                                mode='markers+lines'))
        fRho.add_trace(go.Scatter(x=[r['molecules'] for r in subset], y=[r['density_error_percent'] for r in subset],
                                  name=legend, text=hover, mode='markers+lines',
                                  hovertemplate='%{text}<br>EOS deviation=%{y:+.3f}%<extra></extra>'))
    fD.update_layout(title='H₂ self-diffusion · 95% replica confidence intervals',
                      xaxis_title='H₂ molecules', yaxis_title='D (10⁻⁶ m²/s)',
                      template='plotly_white', height=410, margin=dict(l=55,r=20,t=70,b=55))
    fRho.update_layout(title='Hydrogen density · deviation from reference EOS',
                        xaxis_title='H₂ molecules', yaxis_title='Density deviation (%)',
                        template='plotly_white', height=410, margin=dict(l=55,r=20,t=70,b=55))
    fRho.add_hline(y=0, line_dash='dash', line_color='gray')
    figures = (pio.to_html(fD, full_html=False, include_plotlyjs=True, config={'responsive':True}),
               pio.to_html(fRho, full_html=False, include_plotlyjs=False, config={'responsive':True}))
    trs = []
    for r in rows:
        trs.append('<tr>' + ''.join(f'<td>{html.escape(str(v))}</td>' for v in (
            r['molecules'], f"{r['temperature_K']:.2f}", f"{r['pressure_MPa']:.2f}", r['replicas'],
            f"{r['density_kg_m3']:.4f} ± {r['density_ci95_half_width']:.4f}",
            f"{r['density_error_percent']:+.3f}%",
            f"{r['diffusion_1e6_m2_s']:.4f} ± {r['diffusion_ci95_1e6_m2_s']:.4f}",
            r['n_warnings'])) + '</tr>')
    content = '''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>H₂ MD | Analysis Dashboard</title>
<style>body{margin:0;background:#f3f6fb;font-family:Inter,Segoe UI,Arial,sans-serif;color:#172237}.wrap{max-width:1180px;margin:0 auto;padding:34px 22px}header{background:#172d4b;color:white;padding:35px 0}header h1{font-size:32px;margin:0 0 6px}header p{opacity:.82;margin:0}.eyebrow{text-transform:uppercase;letter-spacing:.12em;font-size:12px;color:#91c6ee}.grid{display:grid;grid-template-columns:1fr 1fr;gap:18px}.panel{background:white;border:1px solid #e5ebf3;border-radius:14px;padding:18px;box-shadow:0 5px 20px #1a34510b;margin-bottom:18px}.chart{min-width:0;overflow:hidden}table{width:100%;border-collapse:collapse;font-size:14px;white-space:nowrap}td,th{padding:12px;border-bottom:1px solid #e5ebf3;text-align:right}th:first-child,td:first-child{text-align:left}th{background:#f8fafc}footer{font-size:13px;color:#5b6b7f;line-height:1.7}.scroll{overflow:auto}@media(max-width:850px){.grid{grid-template-columns:1fr}}</style></head><body>
<header><div class="wrap"><div class="eyebrow">Molecular simulation · Reproducible analytics</div><h1>H₂ Molecular Dynamics Toolkit</h1><p>Campaign summary · Density benchmark · Self-diffusion · Uncertainty</p></div></header>
<main class="wrap"><section class="grid"><div class="panel chart">__FIG1__</div><div class="panel chart">__FIG2__</div></section><section class="panel"><h2>Simulation comparison</h2><div class="scroll"><table><thead><tr><th>H₂ molecules</th><th>T (K)</th><th>P (MPa)</th><th>Runs</th><th>Density (kg/m³) ±95% CI</th><th>EOS error</th><th>D (10⁻⁶ m²/s) ±95% CI</th><th>Warnings</th></tr></thead><tbody>__ROWS__</tbody></table></div></section>
<section class="panel"><h2>Scope and limitations</h2><footer>H₂ molecular centers of mass are used for self-diffusion. Reference density comes from the CoolProp normal-hydrogen EOS. Error bars represent 95% Student-t intervals across independent seed-level estimates only; they do not include force-field errors or finite-size and fitting-window bias. The original force-field paper reports nonbonded two-site LJ parameters; the rigid H–H constraint is a separately documented modeling choice. This dashboard is generated from previously analyzed LAMMPS runs, not from invented data.</footer></section></main></body></html>'''
    output.write_text(content.replace('__FIG1__', figures[0]).replace('__FIG2__', figures[1]).replace('__ROWS__', '\n'.join(trs)), encoding='utf-8')


def dashboard(analyses: list[Path], out: Path) -> list[dict]:
    rows = collect(analyses)
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    with (out/'comparison.csv').open('w', newline='', encoding='utf-8') as file:
        w=csv.DictWriter(file, fieldnames=list(rows[0])); w.writeheader();w.writerows(rows)
    (out/'comparison.json').write_text(json.dumps(rows,indent=2)+'\n', encoding='utf-8')
    _write_markdown(rows,out/'RESULTS.md')
    _create_dashboard(rows,out/'dashboard.html')
    return rows


def main(argv=None):
    parser=argparse.ArgumentParser(prog='h2md-workflow',description='Plan, run, analyze, and present H2 LAMMPS campaigns')
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('plan',help='generate a multi-experiment campaign from JSON')
    p.add_argument('--config',required=True);p.add_argument('--campaign',required=True)
    r=sub.add_parser('run',help='execute prepared LAMMPS jobs sequentially')
    r.add_argument('--campaign',required=True);r.add_argument('--lmp',default='lmp');r.add_argument('--force',action='store_true')
    a=sub.add_parser('analyze',help='analyze all completed campaign runs')
    a.add_argument('--campaign',required=True);a.add_argument('--no-plots',action='store_true')
    d=sub.add_parser('dashboard',help='compare any existing v0.5 analysis results; no LAMMPS required')
    d.add_argument('--analyses',nargs='+',required=True);d.add_argument('--out',default='portfolio_report')
    args=parser.parse_args(argv)
    if args.command=='plan':
        rows=plan(Path(args.config),Path(args.campaign))
        print(f'Prepared {len(rows)} experiment(s) in {args.campaign}; no LAMMPS jobs executed.')
    elif args.command=='run':
        for msg in execute(Path(args.campaign),args.lmp,args.force):print(msg)
    elif args.command=='analyze':
        outputs=analyze_campaign(Path(args.campaign),not args.no_plots)
        print('Analyzed: '+', '.join(str(p) for p in outputs))
    else:
        rows=dashboard([Path(p) for p in args.analyses],Path(args.out))
        print(f'Generated {args.out}/dashboard.html, RESULTS.md, comparison.csv ({len(rows)} experiments)')

if __name__ == '__main__':main()
