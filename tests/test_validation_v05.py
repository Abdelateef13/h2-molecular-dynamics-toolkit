"""Numerical and integration tests for independent-seed H2 density/diffusion workflow."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

V05=Path(__file__).resolve().parents[1]/'lammps'/'paper_h2'/'validation_v05'
sys.path.insert(0,str(V05))
import analyze_replicas as ana  # noqa: E402
import run_replicas as run  # noqa: E402


def test_sim_script_sha_ke_order_and_real_units():
    x=run.render_lammps_input(323.15,10,2027,20000,50000,500,500)
    assert x.index('fix SH all shake') < x.index('fix EQ all npt')
    assert 'units real' in x and '98.69232' in x
    assert 'id mol type xu yu zu' in x
    assert 'lj/cut 12.8000' in x and '0.0199 2.7200' in x
    assert 'output/npt_thermo.dat' in x


def test_setup_replica_inputs_have_independent_seeds(tmp_path):
    cfg=run.prepare(tmp_path,'smoke',[1001,1002,1003],4,323.15,10.0)
    assert cfg['molecules']==256 and cfg['production_steps']==50000
    assert len(set((tmp_path/f'seed_{s}'/'in.h2_validation').read_text() for s in cfg['seeds']))==3
    assert len(set((tmp_path/f'seed_{s}'/'h2_system.data').read_text() for s in cfg['seeds']))==3
    assert (tmp_path/'experiment.json').exists()


def test_msd_conversion_units_and_linear_fit():
    t=np.linspace(.1,500,401)
    desired=2.5e-7
    slope=desired*6/1e-8
    m=10+slope*t
    result=ana.fit_diffusion(t,m,.1,.5,1000)
    assert result['D_m2_s']==pytest.approx(desired,rel=1e-12)
    assert result['r2']==pytest.approx(1.)


def test_msd_drifts_are_removed_when_requested():
    n=41
    # Identical global translation of all molecules: zero relative diffusion.
    frames=np.zeros((n,4,3))
    frames+=np.arange(n)[:,None,None]*np.array([3.0,1.0,0])[None,None,:]
    t,msd=ana.multi_origin_msd(frames,500,origin_stride=2)
    assert np.allclose(msd,0,atol=1e-12)


def test_density_npt_window(tmp_path):
    p=tmp_path/'npt_thermo.dat'
    rows=np.array([[i*500, .006 + .0001*i, 322.+i*.1, 98.69, 99999.] for i in range(1,11)])
    np.savetxt(p,rows,header='TimeStep v_rho v_Tk v_Patm v_V')
    res=ana.density_window(p,.5)
    assert res['n_used']==5
    assert res['density_kg_m3']==pytest.approx(rows[-5:,1].mean()*1000)
    assert res['pressure_MPa']==pytest.approx(98.69/ana.MPA_TO_ATM)


def test_independent_t_confidence_interval():
    a=ana.mean_ci95([1.,2.,3.])
    assert a['mean']==pytest.approx(2.)
    assert a['ci95_half_width']==pytest.approx(4.302652729749*1/np.sqrt(3))


def test_needs_two_replicas_for_ci():
    with pytest.raises(ValueError,match='at least two'):
        ana.mean_ci95([1.])


def test_integration_report_on_synthetic_trajectories(tmp_path,monkeypatch):
    # Synthetic data only for testing the software pipeline, NEVER a physical output.
    cfg=run.prepare(tmp_path/'runs','smoke',[101,102,103],4,323.15,10.)
    rng=np.random.default_rng(42)
    for seed in cfg['seeds']:
        rep=tmp_path/'runs'/f'seed_{seed}'/'output'
        rep.mkdir(parents=True,exist_ok=True)
        np.savetxt(rep/'npt_thermo.dat',[[i*500,.0075+.00002*np.sin(i),323.15,98.6923267,114234.] for i in range(1,21)],header='TimeStep v_rho v_Tk v_Patm v_V')
        # 256 molecules, 512 atoms; dump 100 frames sampled every 500 fs
        initial=rng.uniform(1,49,(256,3))
        displacement=np.cumsum(rng.normal(scale=.2,size=(100,256,3)),axis=0)
        with (rep/'h2_production.lammpstrj').open('w') as f:
            for frame in range(100):
                f.write(f'ITEM: TIMESTEP\n{(frame+1)*500}\nITEM: NUMBER OF ATOMS\n512\nITEM: BOX BOUNDS pp pp pp\n0 55\n0 55\n0 55\nITEM: ATOMS id mol type xu yu zu\n')
                for m in range(256):
                    center=initial[m]+displacement[frame,m]
                    for site in range(2):
                        xyz=center+np.array([.7414*(site-.5),0.,0.])
                        f.write(f'{m*2+site+1} {m+1} 1 {xyz[0]} {xyz[1]} {xyz[2]}\n')
    monkeypatch.setattr(ana,'reference_hydrogen_density',lambda T,P:(7.52,'mock'))
    output=tmp_path/'analysis'
    aggregate=ana.evaluate_runs(tmp_path/'runs',output,plot=True)
    assert len(aggregate['replicas'])==3
    assert aggregate['density_reference_target_kg_m3']==7.52
    assert aggregate['diffusion_m2_s']['mean']>0
    assert (output/'REPORT.md').exists()
    assert (output/'replicas.csv').exists()
    assert (output/'msd_curves.csv').exists()
    assert (output/'density_validation.png').exists()
    assert (output/'diffusion_msd.png').exists()
    assert (output/'replica_comparison.png').exists()


def test_density_accepts_negative_instantaneous_virial_pressure(tmp_path):
    """Instantaneous NPT pressure need not be positive at every sample."""
    p = tmp_path / 'npt_thermo.dat'
    # Retained last six samples average 100 atm but contain negative samples.
    pressures = [50., -300., 300., 50., -400., 400., 150., -150., 300., 300., -250., 250.]
    rows = [[(i + 1) * 500, .008, 323.15, pressures[i], 100000.]
            for i in range(len(pressures))]
    np.savetxt(p, rows, header='TimeStep v_rho v_Tk v_Patm v_V')
    result = ana.density_window(p, .5)
    assert result['n_used'] == 6
    assert result['pressure_MPa'] == pytest.approx(np.mean(pressures[-6:]) / ana.MPA_TO_ATM)
    assert result['negative_pressure_samples_used'] == 2


def test_density_still_rejects_nonphysical_volume(tmp_path):
    p = tmp_path / 'npt_thermo.dat'
    data = np.array([[(i+1)*500, .008, 323.15, -100., 100000.]
                     for i in range(10)])
    data[-1, 4] = -1.
    np.savetxt(p, data)
    with pytest.raises(ValueError, match='density, temperature and volume'):
        ana.density_window(p)


def test_report_handles_nonpositive_mean_pressure_without_fake_eos(tmp_path, monkeypatch):
    # A regression check for a short/unequilibrated NPT segment. We do not
    # treat negative pressures as physical target-pressure equilibrium.
    cfg=run.prepare(tmp_path/'runs', 'smoke', [11,12,13], 4, 323.15, 10.)
    rng=np.random.default_rng(8)
    for seed in cfg['seeds']:
        rep=tmp_path/'runs'/f'seed_{seed}'/'output'
        rep.mkdir(parents=True,exist_ok=True)
        np.savetxt(rep/'npt_thermo.dat',
                   [[(i+1)*500,.008,323.15,-50.,20000.] for i in range(20)])
        xyz=rng.uniform(0,50,(256,3))
        with (rep/'h2_production.lammpstrj').open('w') as f:
            for frame in range(90):
                xyz += rng.normal(scale=.12,size=(256,3))
                f.write(f'ITEM: TIMESTEP\n{(frame+1)*500}\nITEM: NUMBER OF ATOMS\n512\nITEM: BOX BOUNDS pp pp pp\n0 55\n0 55\n0 55\nITEM: ATOMS id mol type xu yu zu\n')
                for m,center in enumerate(xyz):
                    for atom in range(2):
                        pos=center + np.array([.7414*(atom-.5),0,0])
                        f.write(f'{2*m+atom+1} {m+1} 1 {pos[0]} {pos[1]} {pos[2]}\n')
    def eos(t,p):
        assert p>0, 'EOS must not receive a nonpositive mean pressure'
        return (8.0, 'mock')
    monkeypatch.setattr(ana,'reference_hydrogen_density',eos)
    result=ana.evaluate_runs(tmp_path/'runs',tmp_path/'analysis',plot=False)
    assert len(result['replicas']) == 3
    assert all(x['ref_density_measured_TP_kg_m3'] is None for x in result['replicas'])
    assert all(x['density_error_percent_at_measured_TP'] is None for x in result['replicas'])
    assert any('nonpositive' in w for w in result['warnings'])
    assert 'null' in (tmp_path/'analysis'/'summary.json').read_text()


def test_dump_precision_detects_legacy_and_new_format(tmp_path):
    f=tmp_path/'in.h2_validation'
    f.write_text('dump TRAJ all custom 500 output/h2.lammpstrj id mol type xu yu zu\n'
                 'dump_modify TRAJ sort id\n')
    assert ana.dump_coordinate_precision(f)==6
    f.write_text(f.read_text()+'dump_modify TRAJ format float %.16g\n')
    assert ana.dump_coordinate_precision(f)==16


def test_new_generator_writes_high_precision():
    x=run.render_lammps_input(323.15,10,2027,20000,50000,500,500)
    assert 'dump_modify TRAJ format float %.16g' in x


def test_bond_roundoff_at_large_unwrapped_coordinates_is_accepted():
    rng=np.random.default_rng(19)
    mol=300
    dirs=rng.normal(size=(mol,3));dirs/=np.linalg.norm(dirs,axis=-1)[:,None]
    center=rng.uniform(-3000,3000,(4,mol,3))
    h0=center-0.7414*dirs[None]/2
    h1=center+0.7414*dirs[None]/2
    atoms=np.stack((h0,h1),axis=2)
    atoms6=np.vectorize(lambda z:float(f'{z:g}'))(atoms)
    a=ana.bond_precision_audit(atoms6,.7414,6)
    assert a['consistent_with_dump_precision']
    assert a['n_deviations_exceeding_rounding_bound']==0
    assert a['max_deviation_A']>.001


def test_bond_beyond_roundoff_is_rejected():
    coords=np.zeros((2,2,2,3),float)
    coords[:,:,:,0]=np.array([0.,0.8414])[None,None,:]
    a=ana.bond_precision_audit(coords,.7414,6)
    assert not a['consistent_with_dump_precision']
    assert a['n_deviations_exceeding_rounding_bound']>0


def test_precision_audit_high_precision_output():
    coords=np.zeros((2,2,2,3),float)
    coords[:,:,:,0]=np.array([10.,10.7414])[None,None,:]
    a=ana.bond_precision_audit(coords,.7414,16)
    assert a['consistent_with_dump_precision']
    assert a['max_deviation_A']<1e-10
