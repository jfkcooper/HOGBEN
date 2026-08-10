"""Unit tests for polarised code paths.

Covers SimulateReflectivity and Fisher.from_sample. Tests are kept in the
style used elsewhere in the test-suite: they build small refnx structures
and use monkeypatch to avoid heavy I/O.
"""

import numpy as np
import pytest

from refnx.analysis import Parameter
from refnx.reflect import SLD

from hogben.simulate import SimulateReflectivity
from hogben.utils import Fisher
from hogben.models.samples import Sample


def _simple_sample():
    """Return a small Sample instance used in tests.

    Mirrors the style used in other tests in the repository: build a
    compact layer stack and call Sample._vary_structure().
    """
    air = SLD(0, name='Air')
    layer1_thick = Parameter(100, 'Layer 1 Thickness', (50, 120))
    layer1_sld = Parameter(4, 'Layer 1 SLD', (1, 10))
    layer1_thick.optimize = True
    layer1_sld.optimize = True
    layer1 = SLD(layer1_sld, name='Layer 1')(
        thick=layer1_thick, rough=2)
    layer2 = SLD(8, name='Layer 2')(thick=150, rough=2)
    substrate = SLD(2.047, name='Substrate')(thick=0, rough=2)

    structure = air | layer1 | layer2 | substrate

    sample = Sample(structure)
    sample._vary_structure()
    return sample


def test_fisher_from_sample_with_polarised_sample(monkeypatch):
    """Fisher.from_sample should build a Fisher object from a Sample.

    We patch SimulateReflectivity.simulate to return a minimal, well-
    formed dataset so the test is fast and deterministic.
    """
    q = np.linspace(0.01, 0.1, 4)
    r = np.ones_like(q) * 0.1
    dr = np.ones_like(q) * 0.01
    counts = np.ones_like(q) * 100.0
    simulated_stack = np.vstack((q, r, dr, counts))

    monkeypatch.setattr(
        'hogben.utils.SimulateReflectivity.simulate', lambda self:
        simulated_stack,
    )

    sample = _simple_sample()
    angle_times = [(0.7, 4, 1.0)]

    fisher = Fisher.from_sample(sample, angle_times)
    g = fisher.fisher_information

    assert g.shape == (len(sample.get_param_by_attribute('vary')),
                       len(sample.get_param_by_attribute('vary')))
    assert np.all(np.diag(g) >= 0)


@pytest.mark.parametrize('channels', (1, 2))
def test_simulate_reflectivity_handles_channel_counts(monkeypatch, channels):
    """SimulateReflectivity should handle single channel models.

    The test patches the model reflectivity to return either a 1-D array
    (single channel) or a stacked 2-D array (two spin channels). The
    real implementation currently supports single-channel returns; the
    test documents that multi-channel returns raise a ValueError.
    """
    # Patch instrument flux and RNG for determinism.
    wavelengths = np.array([1.0, 2.0])
    flux = np.array([1.0, 1.0])

    monkeypatch.setattr(
        'hogben.simulate.SimulateReflectivity._incident_flux_data',
        lambda self, polarised=False: np.column_stack((wavelengths, flux)),
    )
    monkeypatch.setattr('hogben.simulate.np.random.poisson', lambda lam: lam)

    points = 5

    # Use a real model object (ReflectModel via Sample.get_models) but
    # patch SimulateReflectivity.reflectivity to control output shape.
    sample = _simple_sample()
    model = sample.get_models()[0]

    base_reflectivity = np.full(points, 0.2)

    if channels == 1:
        monkeypatch.setattr(
            'hogben.simulate.SimulateReflectivity.reflectivity',
            lambda self, q: base_reflectivity,
        )
        sim = SimulateReflectivity(model, angle_times=[(0.7, points, 0.1)])

        simulation = sim.simulate(polarised=True)
        assert simulation.ndim == 2
        assert simulation.shape[0] == 4 * channels
        assert simulation.shape[1] == points

        q_out = simulation[0]
        assert np.all(np.isfinite(q_out))
        assert np.all(np.diff(q_out) >= 0)
    else:
        # Multi-channel: reflectivity returns (2, N) stacked array.
        monkeypatch.setattr(
            'hogben.simulate.SimulateReflectivity.reflectivity',
            lambda self, q: np.vstack((base_reflectivity,
                                       base_reflectivity * 0.5)),
        )
        sim = SimulateReflectivity(model, angle_times=[(0.7, points, 0.1)])

        with pytest.raises(ValueError):
            sim.simulate(polarised=True)


def test_incident_flux_called_with_polarised_flag(monkeypatch):
    """Ensure the incident flux loader gets the polarised flag.

    This mirrors the repository test style and avoids file I/O.
    """
    called = {}

    def _fake_incident(self, polarised=False):
        called['polarised'] = polarised
        return np.column_stack((np.array([1.0, 2.0]), np.array([1.0, 1.0])))

    monkeypatch.setattr(
        'hogben.simulate.SimulateReflectivity._incident_flux_data',
        _fake_incident,
    )

    sample = _simple_sample()
    model = sample.get_models()[0]

    sim = SimulateReflectivity(model, angle_times=[(0.7, 3, 0.1)])

    sim.simulate(polarised=True)
    assert called.get('polarised') is True

    called.clear()
    sim.simulate(polarised=False)
    assert called.get('polarised') is False
