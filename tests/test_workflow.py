import json
from pathlib import Path

import pytest

from h2md.workflow import collect, dashboard, plan, _safe_id


def summary(n=256, D=1.75e-6):
    return {
        'molecules': n,
        'temperature_K_target': 323.15,
        'pressure_MPa_target': 10.0,
        'density_reference_target_kg_m3': 7.100073004,
        'density_kg_m3': {'mean': 7.1, 'ci95_half_width': 0.03, 'n': 3},
        'diffusion_m2_s': {'mean': D, 'ci95_half_width': 0.1e-6, 'n': 3},
        'density_relative_error_percent_at_target': -0.001,
        'warnings': []
    }


def test_no_path_traversal():
    for candidate in ('../secret', 'hello/world', '', 'bad space', '.hidden'):
        with pytest.raises(ValueError):
            _safe_id(candidate)


def test_collect_and_dashboard(tmp_path):
    items=[]
    for n in (256, 500, 864):
        p=tmp_path/f'analysis_{n}';p.mkdir()
        (p/'summary.json').write_text(json.dumps(summary(n=n)))
        items.append(p)
    rows=collect(items)
    assert [r['molecules'] for r in rows] == [256,500,864]
    out=tmp_path/'export'
    dashboard(items,out)
    assert '<!doctype html>' in (out/'dashboard.html').read_text(encoding="utf-8")
    assert '256' in (out/'RESULTS.md').read_text(encoding="utf-8")
    assert (out/'comparison.csv').exists()


def test_missing_summary(tmp_path):
    with pytest.raises(FileNotFoundError):
        collect([tmp_path/'missing'])


def test_incompatible_stats(tmp_path):
    p=tmp_path/'a';p.mkdir()
    d=summary();d['diffusion_m2_s']['n']=4
    (p/'summary.json').write_text(json.dumps(d))
    with pytest.raises(ValueError, match='replica count'):
        collect([p])


def test_plan_produces_all_inputs(tmp_path):
    cfg={'defaults':{'preset':'smoke','temperature_K':323.15,'pressure_MPa':10,
                     'seeds':[2027,2028],'dump_every':1000},
         'experiments':[{'id':'h2_256','cells':4}]}
    input_file=tmp_path/'config.json';input_file.write_text(json.dumps(cfg))
    campaign=tmp_path/'campaign'
    experiments=plan(input_file,campaign)
    assert len(experiments)==1
    assert (campaign/'h2_256/runs/seed_2027/in.h2_validation').exists()
    assert (campaign/'h2_256/runs/seed_2028/h2_system.data').exists()
    with pytest.raises(FileExistsError):
        plan(input_file,campaign)


def test_bad_duplicate_ids_rejected(tmp_path):
    cfg={'experiments':[{'id':'a','cells':4},{'id':'a','cells':5}]}
    p=tmp_path/'config.json';p.write_text(json.dumps(cfg))
    with pytest.raises(ValueError,match='Duplicate'):
        plan(p,tmp_path/'campaign')
