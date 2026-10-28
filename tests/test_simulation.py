import numpy as np

from h2md import simulate_h2_centers
from h2md.io import load_trajectory, save_trajectory


def test_md_trajectory_has_unwrapped_positions_and_nve_energy(tmp_path):
    tr = simulate_h2_centers(cells=2, number_density=.3, equilibration_steps=30,
                             production_steps=100, sample_every=5, cutoff=2.0, seed=1)
    assert tr.positions.shape == (21, 32, 3)
    assert tr.frame_dt == .02
    assert np.isfinite(tr.positions).all()
    assert np.any((tr.positions < 0) | (tr.positions >= tr.box_length))
    energies = tr.kinetic_energy + tr.potential_energy
    assert np.ptp(energies) / np.abs(energies.mean()) < 2e-3
    p = tmp_path / 'trajectory.npz'
    save_trajectory(tr, p)
    loaded = load_trajectory(p)
    np.testing.assert_array_equal(loaded.positions, tr.positions)
    np.testing.assert_array_equal(loaded.temperature, tr.temperature)


def test_lammps_load_image_unwrapping(tmp_path):
    content = []
    for step, x, ix in [(0, 9.8, 0), (10, .3, 1), (20, .9, 1)]:
        content.extend(['ITEM: TIMESTEP', str(step), 'ITEM: NUMBER OF ATOMS', '2',
                        'ITEM: BOX BOUNDS pp pp pp', '0 10', '0 10', '0 10',
                        'ITEM: ATOMS id type x y z ix iy iz',
                        '2 1 5.0 5.0 5.0 0 0 0', f'1 1 {x} 2.0 3.0 {ix} 0 0'])
    p = tmp_path / 'example.lammpstrj'
    p.write_text('\n'.join(content) + '\n')
    from h2md.io import load_lammps_dump
    tr = load_lammps_dump(p, timestep_fs=1)
    np.testing.assert_allclose(tr.positions[:, 0, 0], [9.8, 10.3, 10.9])
    assert tr.frame_dt == 10
    assert tr.length_unit == 'angstrom'
