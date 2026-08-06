"""Unit tests covering polarised code paths:
SimulateReflectivity and Fisher.from_sample.

These tests are designed to be fast by mocking instrument data and
RNG calls.
"""

import numpy as np
import pytest

from refnx.analysis import Parameter

from hogben.simulate import SimulateReflectivity
from hogben.utils import Fisher


def test_fisher_from_sample_with_polarised_sample(monkeypatch):
    """Fisher.from_sample should construct a Fisher object from a
    sample, call SimulateReflectivity.simulate for each model, and return
    a valid Fisher information matrix of expected shape.

    We mock SimulateReflectivity.simulate to return a minimal,
    well-formed dataset so this test is fast and deterministic.
    """

    # Prepare a minimal dataset to be returned by simulate().
    q = np.linspace(0.01, 0.1, 4)
    r = np.ones_like(q) * 0.1
    dr = np.ones_like(q) * 0.01
    counts = np.ones_like(q) * 100.0
    simulated_stack = np.vstack((q, r, dr, counts))

    # Patch the simulate method used by Fisher.from_sample to return our
    # data.
    monkeypatch.setattr(
        'hogben.utils.SimulateReflectivity.simulate',
        lambda self: simulated_stack,
    )

    class DummySample:
        def get_models(self):
            # The actual model object is not used because simulate()
            # is mocked.
            return [object()]

        def get_param_by_attribute(self, attr):
            # Return a single refnx Parameter with bounds so
            # Fisher._get_bounds can operate normally.
            return [Parameter(1.0, 'p1', (0.0, 2.0))]

    sample = DummySample()
    angle_times = [(0.7, 4, 1.0)]

    fisher = Fisher.from_sample(sample, angle_times)
    g = fisher.fisher_information

    # Expect a 1x1 Fisher matrix and a non-negative diagonal.
    assert g.shape == (1, 1)
    assert g[0, 0] >= 0


@pytest.mark.parametrize("channels", (1, 2))
def test_simulate_reflectivity_handles_channel_counts(monkeypatch, channels):
    """Ensure SimulateReflectivity.simulate handles models that return
    either a single reflectivity vector (channels=1) or multiple spin
    channels stacked (channels=2).

    The function should return an array whose number of rows equals
    4 * channels (q, r, dr, counts) per channel and columns equal the
    number of q points per angle.
    """

    # Patch instrument flux: small deterministic table
    wavelengths = np.array([1.0, 2.0])
    flux = np.array([1.0, 1.0])
    monkeypatch.setattr(
        'hogben.simulate.SimulateReflectivity._incident_flux_data',
        lambda self, polarised=False: np.column_stack((wavelengths, flux)),
    )

    # Make Poisson deterministic: return the lambda parameter as-is.
    monkeypatch.setattr('hogben.simulate.np.random.poisson', lambda lam: lam)

    points = 5  # small to keep the test fast

    # Define a model that returns either a 1D reflectivity or stacked
    # channels.
    class DummyModel:
        def __call__(self, q):
            base = np.full_like(q, 0.2)
            if channels == 1:
                return base
            else:
                # stack two spin channels (2, N)
                return np.vstack((base, base * 0.5))

    sim = SimulateReflectivity(DummyModel(), angle_times=[(0.7, points, 0.1)])

    # Call simulate with polarised=True to exercise instrument selection.
    simulation = sim.simulate(polarised=True)

    # Expect an array with rows = 4 * channels and columns = points.
    assert simulation.ndim == 2
    assert simulation.shape[0] == 4 * channels
    assert simulation.shape[1] == points

    # The first row corresponds to q; it should be finite and non-decreasing.
    q_out = simulation[0]
    assert np.all(np.isfinite(q_out))
    assert np.all(np.diff(q_out) >= 0)


def test_incident_flux_called_with_polarised_flag(monkeypatch):
    """Ensure the instrument flux loader is called with the same
    polarised flag passed to SimulateReflectivity.simulate.
    """

    called = {}

    def fake_incident(self, polarised=False):
        called['polarised'] = polarised
        # return a minimal valid array
        return np.column_stack((np.array([1.0, 2.0]), np.array([1.0, 1.0])))

    monkeypatch.setattr(
        'hogben.simulate.SimulateReflectivity._incident_flux_data',
        fake_incident,
    )

    class DummyModel:
        def __call__(self, q):
            return np.full_like(q, 0.1)

    sim = SimulateReflectivity(DummyModel(), angle_times=[(0.7, 3, 0.1)])

    # simulate with polarised=True
    sim.simulate(polarised=True)
    assert called.get('polarised') is True

    # simulate with polarised=False
    called.clear()
    sim.simulate(polarised=False)
    assert called.get('polarised') is False
