"""Version 0.4 tests for the paper-sourced two-site H2 configuration.
Run from the ROOT of the h2-md repository: python -m pytest -q
"""
import importlib.util
import json
from pathlib import Path
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
FILE=ROOT/'lammps'/'paper_h2'/'build_h2_paper.py'
spec=importlib.util.spec_from_file_location('build_h2_paper',FILE)
builder=importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


def test_paper_nonbonded_parameters_exact():
    assert builder.SIGMA_A == 2.7200
    assert builder.EPS_KCAL_MOL == 0.0199
    assert builder.CHARGE_E == 0.0
    assert builder.CUTOFF_A == 12.8
    assert builder.TIMESTEP_FS == 1.0


def test_generation_geometry_and_box():
    x,L = builder.build_system(4, 323.15, 10.0, seed=1)
    assert x.shape == (256,2,3)
    assert np.all(x>=0) and np.all(x<L)
    assert np.max(np.abs(np.linalg.norm(x[:,1]-x[:,0],axis=1)-builder.BOND_A))<1e-11
    assert L>2*builder.CUTOFF_A
    assert builder.build_system(4,323.15,10.0,seed=1)[1]==L


def test_pressure_box_scaling_and_guard():
    L10=builder.ideal_gas_box_angstrom(256,10,323.15)
    L20=builder.ideal_gas_box_angstrom(256,20,323.15)
    assert L20/L10 == pytest.approx(2**(-1/3))
    with pytest.raises(ValueError,match='too small'):
        builder.build_system(4,323.15,60.0)
    assert builder.build_system(5,323.15,60.0)[1]>2*builder.CUTOFF_A+2


def test_h2_lammps_input_configuration():
    script=builder.render_input(323.15,10.0,20000,50000,500)
    for token in ['units           real','pair_style      lj/cut 12.8000',
                  'pair_coeff      1 1 0.0199 2.7200 12.8000',
                  'fix             SH all shake','fix             EQ all npt',
                  'fix             PROD all nvt','special_bonds   lj 0.0 0.0 0.0',
                  'dump            TR all custom 500', 'xu yu zu',
                  'unfix           EQ','timestep        1.0']:
        assert token in script


def test_write_data_and_metadata(tmp_path):
    coords,L=builder.build_system(4,323.15,10)
    path=tmp_path/'h2_system.data'
    builder.write_data(path,coords,L)
    s=path.read_text()
    assert '512 atoms' in s and '256 bonds' in s
    assert s.count('\n1 ')>=2  # first atom and first bond appear
    assert len([line for line in s.split('Atoms # molecular')[1].split('Bonds')[0].splitlines() if line.strip()]) == 512


def test_reader_accepts_synthetic_h2_dump(tmp_path):
    # Generate physically consistent 2-frame trajectory, including PBC crossing.
    from h2md.atomistic import load_atomistic_lammps_dump
    dump=tmp_path/'h2_paper.lammpstrj'
    with dump.open('w') as f:
        for t in (0,500):
            f.write(f'ITEM: TIMESTEP\n{t}\nITEM: NUMBER OF ATOMS\n4\n'
                     'ITEM: BOX BOUNDS pp pp pp\n0 35\n0 35\n0 35\n'
                     'ITEM: ATOMS id mol type xu yu zu\n')
            f.write(f'1 1 1 {34.9+t/500:.7f} 1.0 1.0\n')
            f.write(f'2 1 1 {34.9+t/500+0.7414:.7f} 1.0 1.0\n')
            f.write(f'3 2 1 {5+t/500:.7f} 5.0 5.0\n')
            f.write(f'4 2 1 {5+t/500+0.7414:.7f} 5.0 5.0\n')
    tr=load_atomistic_lammps_dump(dump,timestep_fs=1)
    assert tr.n_frames==2 and tr.n_molecules==2
    assert tr.frame_dt==500
    np.testing.assert_allclose(tr.bond_lengths(),.7414,atol=1e-6)
    np.testing.assert_allclose(tr.com[1,:,0]-tr.com[0,:,0],1.0)
