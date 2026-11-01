"""Checks: proper atomic topology, PBC, conservation, LAMMPS mapping."""
import numpy as np
import pytest

from h2md.atomistic import (AtomicH2Trajectory, load_atomistic_lammps_dump,
                            load_atomistic_trajectory, save_atomistic_trajectory,
                            simulate_atomistic_h2)
from h2md.visualization import display_positions, write_extended_xyz


def test_atomistic_has_two_atoms_each_and_conserves_energy():
    t=simulate_atomistic_h2(equilibration_steps=100,production_steps=300,sample_every=10,seed=10)
    assert t.atoms.shape == (31,32,2,3)
    bonds=t.bond_lengths()
    assert abs(np.mean(bonds)-.3)<.05
    assert np.max(np.abs(bonds-.3))<.15
    energies=t.potential_energy+t.kinetic_energy
    assert np.ptp(energies)/abs(np.mean(energies)) < .02
    com=t.center_trajectory()
    np.testing.assert_allclose(com.positions,t.atoms.mean(axis=2))


def test_atomistic_save_and_xyz(tmp_path):
    t=simulate_atomistic_h2(equilibration_steps=20,production_steps=30,sample_every=5)
    p=tmp_path/'atoms.npz'
    save_atomistic_trajectory(t,p)
    u=load_atomistic_trajectory(p)
    np.testing.assert_array_equal(t.atoms,u.atoms)
    x=display_positions(t,0)
    np.testing.assert_allclose(np.linalg.norm(x[:,0]-x[:,1],axis=1),t.bond_lengths()[0])
    p=write_extended_xyz(t,tmp_path/'atoms.xyz',stride=2)
    txt=p.read_text()
    assert txt.startswith('64\n')
    assert txt.count('Lattice=')==4
    assert 'Properties=species:S:1:pos:R:3:molecule:I:1' in txt


def test_lammps_atomistic_parser_groups_by_mol_and_unwraps(tmp_path):
    out=[]
    for step,left,right in [(0,9.5,9.9),(10,10.0,10.4),(20,10.5,10.9)]:
        out+=['ITEM: TIMESTEP',str(step),'ITEM: NUMBER OF ATOMS','4',
              'ITEM: BOX BOUNDS pp pp pp','0 10','0 10','0 10',
              'ITEM: ATOMS id mol xu yu zu',
              # deliberately permuted atom order; mol groups correctly
              '4 2 5.4 5 5',f'2 1 {right} 2 3',
              '3 2 5.0 5 5',f'1 1 {left} 2 3']
    p=tmp_path/'atoms.lammpstrj';p.write_text('\n'.join(out)+'\n')
    t=load_atomistic_lammps_dump(p,timestep_fs=.5)
    assert t.atoms.shape==(3,2,2,3)
    assert t.frame_dt==5.0
    np.testing.assert_allclose(t.com[:,0,0],[9.7,10.2,10.7])
    np.testing.assert_allclose(t.bond_lengths(),.4)


def test_invalid_single_atom_molecule_rejected(tmp_path):
    content='\n'.join(['ITEM: TIMESTEP','0','ITEM: NUMBER OF ATOMS','4',
      'ITEM: BOX BOUNDS pp pp pp','0 10','0 10','0 10',
      'ITEM: ATOMS id mol xu yu zu',
      '1 1 1 1 1','2 2 2 2 2','3 2 3 3 3','4 2 4 4 4'])
    p=tmp_path/'bad.lammpstrj';p.write_text(content+'\n')
    with pytest.raises(ValueError,match='exactly 2'):
        load_atomistic_lammps_dump(p,timestep_fs=.5)
