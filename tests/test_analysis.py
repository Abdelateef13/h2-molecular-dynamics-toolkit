import numpy as np
import pytest

from h2md import Trajectory, fit_self_diffusion, mean_squared_displacement, radial_distribution
from h2md.analysis import diffusion_in_si, mass_density_g_cm3, number_density


def test_constant_velocity_msd_is_quadratic():
    t = np.arange(50) * 0.2
    v = np.array([[0.5, 1.0, 0], [-0.7, .3, .1], [.2, .1, 1]])
    pos = t[:, None, None] * v[None, :, :]
    tr = Trajectory(pos, box_length=10, frame_dt=.2)
    times, msd = mean_squared_displacement(tr, max_lag=20)
    expected = np.mean((v * v).sum(axis=1)) * times**2
    np.testing.assert_allclose(msd, expected, rtol=1e-12, atol=1e-13)


def test_brownian_diffusion_fit_recovers_known_D():
    rng = np.random.default_rng(6)
    n_frames, n_particles, D, dt = 150, 600, 0.22, 0.05
    jumps = rng.normal(scale=np.sqrt(2*D*dt), size=(n_frames - 1, n_particles, 3))
    coords = np.concatenate([np.zeros((1, n_particles, 3)), jumps.cumsum(axis=0)], axis=0)
    t, msd = mean_squared_displacement(Trajectory(coords, 20, dt), max_lag=40)
    fitted = fit_self_diffusion(t, msd)
    assert fitted.coefficient == pytest.approx(D, rel=0.12)
    assert fitted.r_squared > 0.99


def test_uniform_ideal_gas_rdf_baseline():
    rng = np.random.default_rng(3)
    positions = rng.uniform(0, 15, size=(25, 120, 3))
    r, g = radial_distribution(Trajectory(positions, 15, 1), bins=28, frame_stride=1)
    assert r.shape == g.shape == (28,)
    assert np.mean(g[(r > 2) & (r < 6)]) == pytest.approx(1, abs=.08)


def test_real_units_and_density():
    tr = Trajectory(np.zeros((5, 4, 3)), box_length=10, frame_dt=5,
                    length_unit="angstrom", time_unit="fs")
    assert number_density(tr) == pytest.approx(.004)
    assert mass_density_g_cm3(tr) == pytest.approx(.01339, rel=.01)
    fit = fit_self_diffusion(np.arange(10.), 12 * np.arange(10.))
    assert diffusion_in_si(fit, tr) == pytest.approx(2e-5)


def test_trajectory_validation():
    with pytest.raises(ValueError):
        Trajectory(np.zeros((2, 3, 2)), 10, 1)
    with pytest.raises(ValueError):
        Trajectory(np.zeros((2, 3, 3)), -10, 1)
    with pytest.raises(ValueError):
        Trajectory(np.zeros((2, 3, 3)), 10, 0)
