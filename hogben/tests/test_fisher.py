"""Unit tests for the methods that relate directly to the Fisher information"""

import copy

import numpy as np
import pytest
import refnx

from hogben.models.samples import Sample
from hogben.utils import Fisher
from refnx.reflect import SLD
from unittest.mock import Mock, patch


QS = [np.array([0.1, 0.2, 0.4, 0.6, 0.8])]
COUNTS = [np.ones(len(QS[0])) * 100]


@pytest.fixture
def model():
    """Define a bilayer sample, and return the associated refnx model"""
    # Define sample
    air = SLD(0, name='Air')
    layer1 = SLD(4, name='Layer 1')(thick=60, rough=8)
    layer2 = SLD(8, name='Layer 2')(thick=150, rough=2)
    substrate = SLD(2.047, name='Substrate')(thick=0, rough=2)
    layer1.thick.bounds = (50, 70)
    layer2.thick.bounds = (140, 160)
    layer1.rough.bounds = (7, 9)
    layer2.rough.bounds = (1, 3)
    sample = air | layer1 | layer2 | substrate
    model = refnx.reflect.ReflectModel(sample, scale=1, bkg=5e-6, dq=2)
    # Define model
    model.xi = [layer1.rough, layer2.rough, layer1.thick, layer2.thick]
    return model


@pytest.fixture
def mock_model():
    """
    Create a mock of the refnx model with a given set of parameters and their
    bounds
    """
    # Parameters described as tuples: (value, lower bound, upper bound)
    parameter_values = [
        (20, 15, 25),
        (50, 45, 55),
        (10, 7.5, 8.5),
        (2, 1.5, 2.5),
    ]

    # Fill parameter values and bounds
    parameters = [
        Mock(
            spec=refnx.analysis.Parameter,
            value=value,
            bounds=Mock(lb=lb, ub=ub),
        )
        for value, lb, ub in parameter_values
    ]
    model = Mock(spec=refnx.reflect.ReflectModel, xi=parameters)
    model.xi = parameters
    return model


@pytest.fixture
def refnx_two_solvents():
    """Defines a structure describing a simple sample with two solvents"""
    H2O = SLD(-0.52, name='H2O')
    D2O = SLD(6.19, name='D2O')
    layer1 = SLD(4, name='Layer 1')(thick=60, rough=8)
    layer2 = SLD(8, name='Layer 2')(thick=150, rough=2)
    substrate = SLD(2.047, name='Substrate')(thick=0, rough=2)
    structure_H2O = H2O | layer1 | layer2 | substrate
    structure_D2O = D2O | layer1 | layer2 | substrate
    return [structure_H2O, structure_D2O]


def generate_reflectivity_data():
    """
    Generates predefined reflectivity.The reflectivity values are yielded
    alternatingly between two predefined lists of reflectivity values,
    simulating a change in reflectivity between two data points
    """
    r = [[1.0, 0.5, 0.4, 0.2, 0.1], [0.95, 0.45, 0.35, 0.15, 0.05]]
    while True:
        yield r[0]
        yield r[1]


def test_fisher_workflow(model):
    """
    Runs the entire fisher workflow for the refnx model, and checks that the
    corresponding results are consistent with the expected values
    """
    g = Fisher(QS, model.xi, COUNTS,
               [model]).fisher_information
    expected_fisher = [
        [5.17704306e-06, 2.24179068e-06, -5.02221954e-07, -7.91886209e-07],
        [2.24179068e-06, 1.00559528e-06, -2.09433754e-07, -3.18583142e-07],
        [-5.02221954e-07, -2.09433754e-07, 5.75647233e-08, 1.03142100e-07],
        [-7.91886209e-07, -3.18583142e-07, 1.03142100e-07, 1.99470835e-07],
    ]
    np.testing.assert_allclose(g, expected_fisher, rtol=1e-08)


def test_fisher_multiple_backgrounds(refnx_two_solvents):
    """
    Tests whether the Fisher information still behaves as expected when using
    two different backgrounds in a sample.
    """
    sample = Sample(refnx_two_solvents, bkg=[1e-6, 5e-6])
    sample._vary_structure()
    angle_times = [(0.7, 100, 100000), (2.0, 100, 100000)]
    fisher = Fisher.from_sample(sample, angle_times)
    eigenval_ref = 257.9191161169581
    np.testing.assert_allclose(eigenval_ref,
                               fisher.min_eigenval, rtol=1e-08)


@patch('hogben.utils.SimulateReflectivity.reflectivity')
def test_fisher_analytical_values(mock_reflectivity, mock_model):
    """
    Tests that the values of the calculated Fisher information matrix (FIM)
    are calculated correctly when no importance scaling is given.

    The FIM is calculated using matrix multiplication given:
    g = J.T x M x J


    Where J describes the Jacobian of the reflectance with respect to the
    parameter value, and M describe the diagonal matrix of the incident count
    divided by the model reflectances. J.T describes the transpose of J.

    For this unit test the values for J and M are known:
    J = [-0.25, -0.1 , -0.5 ],
        [-0.25, -0.1 , -0.5 ],
        [-0.25, -0.1 , -0.5 ],
        [-0.25, -0.1 , -0.5 ],
        [-0.25, -0.1 , -0.5 ]

    M is given by by:
    M =  [ 100.,    0.,    0.,    0.,    0.],
         [   0.,  200.,    0.,    0.,    0.],
         [   0.,    0.,  250.,    0.,    0.],
         [   0.,    0.,    0.,  500.,    0.],
         [   0.,    0.,    0.,    0., 1000.]

    Resulting in g = J.T x M x J:
    g =  [128.125,  51.25 , 256.25 ],
       [ 51.25 ,  20.5  , 102.5  ],
       [256.25 , 102.5  , 512.5  ]

    After this, the elements are scaled to the bounds of each unit. Using:
    g_scaled = H.T * g *H
    Where H is a diagonal matrix where the diagonal elements for each
    parameter are given by H_ij = 1/(upper_bound - lower_bound), resulting in:
    H = [0.1, 0. , 0. ],
       [0. , 0.1, 0. ],
       [0. , 0. , 1. ]

    Which should finally result in the FIM matrix equal to:
    g = [1.28125, 0.5125, 25.625],
        [0.5125, 0.205, 10.25],
        [25.625, 10.25, 512.5]
    """
    xi = mock_model.xi[:3]
    mock_reflectivity.side_effect = generate_reflectivity_data()
    g_correct = [
        [1.28125, 0.5125, 25.625],
        [0.5125, 0.205, 10.25],
        [25.625, 10.25, 512.5],
    ]

    g_reference = Fisher(QS, xi, COUNTS, [mock_model]).fisher_information
    np.testing.assert_allclose(g_reference, g_correct, rtol=1e-08)


@patch('hogben.utils.SimulateReflectivity.reflectivity')
def test_fisher_importance_scaling(mock_reflectivity, mock_model):
    """
    Tests that the values of the calculated Fisher information matrix
    are calculated correctly when an importance scaling is applied.

    The importance scaling is applied by scaling each parameter of the FIM
    to a given importance value using g_scaled = g * importance
    Where g is the unscaled FIM, and importance is a diagonal matrix with
    the importance scaling of each parameter on the diagonals. For this unit
    test the importance matrix is equal to:
    importance = [1, 0 , 0]
                 [0, 2, 0]
                [0, 0, 3]
    Yielding a FIM where every column should be scaled by the corresponding
    diagonal in the importance matrix:
    g = [1.28125, 1.025, 76.875],
        [0.5125, 0.41, 30.75],
        [25.625, 20.5, 1537.5]
    """
    xi = mock_model.xi[:3]
    for index, param in enumerate(xi):
        param.importance = index + 1
    mock_reflectivity.side_effect = generate_reflectivity_data()
    g_correct = [
        [1.28125, 1.025, 76.875],
        [0.5125, 0.41, 30.75],
        [25.625, 20.5, 1537.5],
    ]
    g_reference = Fisher(QS, xi, COUNTS, [mock_model]).fisher_information
    np.testing.assert_allclose(g_reference, g_correct, rtol=1e-08)


@pytest.mark.parametrize('step', (0.01, 0.0075, 0.0025, 0.001, 0.0001))
def test_fisher_consistent_steps(step, model):
    """
    Tests whether the Fisher information remains mostly consistent when
    changing step size using the refnx model
    """
    model.xi = model.xi[:1]
    g_reference = Fisher(QS, model.xi, COUNTS, [model],
                         step=0.005).fisher_information
    g_compare = Fisher(QS, model.xi, COUNTS, [model],
                       step=step).fisher_information
    np.testing.assert_allclose(g_reference, g_compare, rtol=1e-02)


@patch('hogben.utils.SimulateReflectivity.reflectivity')
@pytest.mark.parametrize('model_params', (1, 2, 3, 4))
def test_fisher_shape(mock_reflectivity, model_params, mock_model):
    """
    Tests whether the shape of the Fisher information matrix remains
     correct when changing the amount of parameters
    """
    xi = mock_model.xi[:model_params]

    mock_reflectivity.side_effect = generate_reflectivity_data()

    expected_shape = (model_params, model_params)

    g = Fisher(QS, xi, COUNTS, [mock_model]).fisher_information
    np.testing.assert_array_equal(g.shape, expected_shape)


@patch('hogben.utils.SimulateReflectivity.reflectivity')
@pytest.mark.parametrize(
    'qs',
    (
        np.arange(0.001, 1.0, 0.25),
        np.arange(0.001, 1.0, 0.10),
        np.arange(0.001, 1.0, 0.05),
        np.arange(0.001, 1.0, 0.01),
    ),
)
def test_fisher_diagonal_non_negative(mock_reflectivity,
                                      qs, mock_model):
    """Tests whether the diagonal values in the Fisher information matrix
    are all zero or greater"""
    mock_reflectivity.side_effect = (np.random.rand(len(qs)) for _ in range(9))
    counts = [np.ones(len(qs)) * 100]
    g = Fisher([qs], mock_model.xi, counts, [mock_model]).fisher_information
    assert np.all(np.diag(g)) >= 0


@pytest.mark.parametrize('model_params', (1, 2, 3, 4))
def test_fisher_no_data(model_params, mock_model):
    """Tests whether a model with zero data points properly returns an empty
    matrix of the correct shape"""
    xi = mock_model.xi[:model_params]
    g = Fisher([], xi, COUNTS, [mock_model]).fisher_information
    np.testing.assert_equal(g, np.zeros((len(xi), len(xi))))


@patch('hogben.utils.SimulateReflectivity.reflectivity')
def test_fisher_no_parameters(mock_reflectivity, mock_model):
    """Tests whether a model without any parameters properly returns a
    zero array"""
    mock_reflectivity.side_effect = generate_reflectivity_data()
    g = Fisher(QS, [], COUNTS, [model]).fisher_information
    np.testing.assert_equal(g.shape, (0, 0))


def test_fisher_doubling_with_two_identical_models(model):
    """
    Tests that using two identical models with the same q-points and counts
    correctly doubles the values on the elements in the Fisher information
    matrix
    """
    g_single = Fisher(QS, model.xi, COUNTS, [model],
                      0.005).fisher_information
    counts = [COUNTS[0], COUNTS[0]]
    qs = [QS[0], QS[0]]
    g_double = Fisher(qs, model.xi, counts, [model, model],
                      0.005).fisher_information
    np.testing.assert_allclose(g_double, g_single * 2, rtol=1e-08)


def test_multiple_models_shape(model):
    """
    Tests that shape of the Fisher information matrix is equal to the total
    sum of parameters over all models.
    """
    model_2 = copy.deepcopy(model)
    model_2.xi = model_2.xi[:-1]
    xi = model.xi + model_2.xi
    xi_length = len(model.xi) + len(model_2.xi)
    counts = [COUNTS[0], COUNTS[0]]
    qs = [QS[0], QS[0]]
    g_double = Fisher(qs, xi, counts, [model, model_2],
                      0.005).fisher_information
    np.testing.assert_equal(g_double.shape, (xi_length, xi_length))
