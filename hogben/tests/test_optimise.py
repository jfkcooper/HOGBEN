"""Unit tests for the methods in the optimise module"""
import pytest
import numpy as np

from hogben.optimise import Optimiser, optimise_parameters
from refnx.analysis import Parameter
from refnx.reflect import SLD
from hogben.models.samples import Sample
from hogben.models.bilayers import BilayerDMPC
from unittest.mock import MagicMock, patch


@pytest.fixture
def refnx_sample():
    """Defines a structure describing a simple sample."""
    air = SLD(0, name='Air')
    layer1_thick = Parameter(100, 'Layer 1 Thickness', (50, 120))
    layer1_sld = Parameter(4, 'Layer 1 SLD', (50, 120))
    layer1_thick.optimize = True
    layer1_sld.optimize = True
    layer1 = SLD(layer1_sld, name='Layer 1')(thick=layer1_thick, rough=2)
    layer2 = SLD(8, name='Layer 2')(thick=150, rough=2)
    substrate = SLD(2.047, name='Substrate')(thick=0, rough=2)
    structure = air | layer1 | layer2 | substrate
    sample = Sample(structure)
    sample._vary_structure()
    return sample


class MockSample:
    """
    Mocks a HOGBEN sample, so that class methods can be mocked
    and tracked
    """
    def __init__(self):
        """Iniitalize the mocked sample"""
        self.params = [MagicMock(name='param1'), MagicMock(name='param2')]
        self.simulate_reflectivity = MagicMock()
        self.sld_profile = MagicMock()
        self.get_models = MagicMock()
        self.is_magnetic = MagicMock()
        self.polarised = False

    def get_param_by_attribute(self, _attribute):
        """
        Return mocked parameters when requesting parameters from attribute
        """
        return self.params


class MockFisher:
    """Mocks the Fisher class, to return a simple minimum eigenvalue of 1.0 """

    def __init__(self, _sample, _angle_times, inst_or_path='OFFSPEC'):
        """Initialize the mocked Fisher class and set eigenval to 1"""
        self.min_eigenval = 1.0

    @staticmethod
    def from_sample(sample, angle_times, inst_or_path='OFFSPEC'):
        """Create a Mocked Fisher object from sample"""
        return MockFisher(sample, angle_times, inst_or_path)


@patch('hogben.optimise.Optimiser._Optimiser__optimise')
def test_optimise_angle_times_length(mock_optimise, refnx_sample):
    """
    Tests that the optimise_angle_times method outputs the correct amount of
    angles and counting times.
    """
    num_angles = 2
    optimiser = Optimiser(refnx_sample)

    # Mock values retrieved from previous run
    mock_optimise.return_value = np.array([0.8847156, 0.88834418,
                                           0.00139696,
                                           0.99860304]), -0.7573710562837207
    angles, splits, _ = optimiser.optimise_angle_times(num_angles,
                                                       angle_bounds=(0.2, 2.3),
                                                       verbose=False)
    assert len(angles) == num_angles and len(splits) == num_angles


@patch('hogben.optimise.Optimiser._Optimiser__optimise')
def test_optimise_contrasts(mock_optimise):
    """
    Tests that the optimise_contrasts method outputs the correct amount of
    contrasts and counting times.
    """
    optimiser = Optimiser(BilayerDMPC())
    num_contrasts = 3
    angle_times = [(0.7, 100, 10), (2.3, 100, 40)]

    # Mock values retreived from previous run
    mock_optimise.return_value = (
        np.array([-0.56, 2.15, 6.36, 0.17, 0.28, 0.56]), -0.18
    )
    contrasts, splits, _ = optimiser.optimise_contrasts(num_contrasts,
                                                        angle_times,
                                                        workers=-1,
                                                        verbose=False)
    assert len(contrasts) == num_contrasts and len(splits) == num_contrasts


@patch('hogben.optimise.Optimiser._Optimiser__optimise')
def test_optimise_underlayers(mock_optimise):
    """
    Tests that the optimise_contrasts method outputs the correct amount of
    contrasts and counting times.
    """
    optimiser = Optimiser(BilayerDMPC())

    num_underlayers = 3
    angle_times = [(0.7, 100, 10), (2.3, 100, 40)]
    contrasts = [-0.56, 6.36]
    thick_bounds = (0, 500)
    sld_bounds = (1, 9)

    # Mock values retreived from previous run
    mock_optimise.return_value = (
        np.array([-0.56, 2.15, 6.36, 0.17, 0.28, 0.56]), -0.18
    )

    contrasts, splits, _ = optimiser.optimise_underlayers(num_underlayers,
                                                          angle_times,
                                                          contrasts,
                                                          thick_bounds,
                                                          sld_bounds,
                                                          verbose=False)
    assert len(contrasts) == num_underlayers and len(splits) == num_underlayers


def test_angle_times_func_result(refnx_sample):
    """Checks that the angle_times_func method gives the correct result"""
    angle_time_split = [0.3, 1.3, 0.8, 0.2]  # [angle, angle, time, time]
    num_angles = 2
    contrasts = [3, 14, -2]
    points = 100
    total_time = 10000

    optimiser = Optimiser(refnx_sample)
    result = optimiser._angle_times_func(angle_time_split, num_angles,
                                         contrasts, points, total_time)

    expected_result = -1.7716709840530174
    np.testing.assert_allclose(result, expected_result, rtol=1e-06)


def test_contrasts_func_result():
    """Checks that the _contrasts_func method gives the correct result"""
    contrasts_time = [0.3, 9.3, 0.8, 0.2]  # [SLD, SLD, time, time]
    num_contrasts = 2
    angle_splits = [(0.7, 100, 0.6), (2.3, 100, 0.4)]
    total_time = 100000

    optimiser = Optimiser(BilayerDMPC())
    result = optimiser._contrasts_func(contrasts_time, num_contrasts,
                                       angle_splits, total_time)

    expected_result = -0.199884
    np.testing.assert_allclose(result, expected_result, rtol=1e-06)


def test_underlayers_func():
    """Checks that the _underlayers_func method gives the correct result"""
    thickness_SLD = [50, 20, -5, 10]  # [thickness, thickness, SLD, SLD]

    bilayer = BilayerDMPC()
    optimiser = Optimiser(bilayer)
    num_underlayers = 2
    contrasts = [-0.56, 6.36]

    angle_times = [(0.7, 100, 10000), (2.3, 100, 10000)]
    result = optimiser._underlayers_func(thickness_SLD, num_underlayers,
                                         angle_times, contrasts)

    expected_result = -1.500100628
    np.testing.assert_allclose(result, expected_result, rtol=1e-06)


@patch('hogben.optimise.Optimiser.optimise_parameters')
@patch('hogben.optimise.scan_parameters')
@patch('hogben.optimise.Fisher.from_sample', new=MockFisher)
def test_optimise_parameters(mock_scan_parameters, mock_optimise_parameters):
    """
    Runs the general `optimise_parameters` workflow, and makes sure the output
    is handled correctly, and data visualization is called.
    """
    mock_optimise_parameters.return_value = ([0.5, 0.5], 2.0)
    mock_scan_parameters.return_value = None

    sample = MockSample()
    angle_times = [1, 2, 3]
    result = optimise_parameters(sample, angle_times, visualise=True)

    # Check if the sample parameters were updated correctly
    assert result.params[0].value == 0.5
    assert result.params[1].value == 0.5

    # Check whether the visualisation functions are called
    sample.sld_profile.assert_called_once()
    sample.simulate_reflectivity.assert_called_once()
    mock_scan_parameters.assert_called_once()


def test_optimise_params_length(refnx_sample):
    """
    Runs the optimiser workflow and tests that the optimise_parameters method
    outputs the correct amount of counting parameters.
    """
    angle_times = [(0.7, 100, 10000), (2.3, 100, 10000)]
    params = refnx_sample.get_param_by_attribute('optimize')
    optimiser = Optimiser(refnx_sample)

    # Mock values retrieved from previous run
    values, _ = optimiser.optimise_parameters(angle_times, verbose=False)
    assert len(values) == len(params)


def test_parameter_func(refnx_sample):
    """Checks that the _underlayers_func method gives the correct result"""
    optimiser = Optimiser(refnx_sample)
    params = refnx_sample.get_param_by_attribute('optimize')
    params.sort(key=lambda x: x.name)
    values = [3.3, 109]
    angle_times = [(0.7, 100, 10000), (2.3, 100, 10000)]
    result = optimiser._parameter_func(values, params, angle_times)
    expected_result = -15.438106789796723
    np.testing.assert_allclose(result, expected_result, rtol=1e-06)
