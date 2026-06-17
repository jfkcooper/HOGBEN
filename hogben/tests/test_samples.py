import os
import tempfile

import numpy as np
import pytest
import matplotlib
import hogben.models.samples as samples

from hogben.models.base import MagneticSLD
from hogben.models.samples import Sample
from hogben.simulate import SimulateReflectivity
from hogben.utils import Fisher
from refnx.reflect import SLD
from unittest.mock import Mock, patch


@pytest.fixture
def refnx_sample():
    """Defines a structure describing a simple sample."""
    air = SLD(0, name='Air')
    layer1 = SLD(4, name='Layer 1')(thick=60, rough=8)
    layer2 = SLD(8, name='Layer 2')(thick=150, rough=2)
    substrate = SLD(2.047, name='Substrate')(thick=0, rough=2)
    structure = air | layer1 | layer2 | substrate
    return Sample(structure)


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


@pytest.fixture
def refnx_three_solvents():
    """Defines a structure describing a simple sample with three solvents"""
    H2O = SLD(-0.52, name='H2O')
    D2O = SLD(6.19, name='D2O')
    SMW = SLD(2.07, name='SMW')
    layer1 = SLD(4, name='Layer 1')(thick=60, rough=8)
    layer2 = SLD(8, name='Layer 2')(thick=150, rough=2)
    substrate = SLD(2.047, name='Substrate')(thick=0, rough=2)
    structure_H2O = H2O | layer1 | layer2 | substrate
    structure_D2O = D2O | layer1 | layer2 | substrate
    structure_SMW = SMW | layer1 | layer2 | substrate
    return [structure_H2O, structure_D2O, structure_SMW]


@pytest.fixture
def refnx_magnetic_structure():
    """Defines a structure describing a sample with one magnetic layer."""
    air = SLD(0, name='Air')
    layer1 = SLD(3, name='Layer 1')(thick=60, rough=8)
    layer2 = SLD(8, name='Layer 2')(thick=150, rough=2)
    mag_layer = MagneticSLD(6, 2, name='Magnetic layer')(thick=75, rough=2)
    substrate = SLD(2.047, name='Substrate')(thick=0, rough=2)
    structure = air | layer1 | layer2 | mag_layer | substrate
    return structure


@pytest.fixture
def refnx_magnetic_structure_multiple_layers():
    """Defines a structure describing a sample with one magnetic layer."""
    air = SLD(0, name='Air')
    D2O = SLD(6.5, name='Air')
    layer1 = SLD(3, name='Layer 1')(thick=60, rough=8)
    layer2 = MagneticSLD(8, 2, name='Layer 2')(thick=150, rough=2)
    mag_layer = MagneticSLD(6, 2, name='Magnetic layer')(thick=75, rough=2)
    substrate = SLD(2.047, name='Substrate')(thick=0, rough=2)
    structure_D2O = D2O | layer1 | layer2 | mag_layer | substrate
    structure = air | layer1 | layer2 | mag_layer | substrate
    return [structure, structure_D2O]


def mock_save_plot(fig: matplotlib.figure.Figure,
                   save_path: str,
                   filename: str) -> None:
    """
    A mocked version of the hogben.utils.save_plot method, where a lower
    dpi is used when saving a figure

    Args:
        fig: The matplotlib figure to be plotted
        save_path: The path where the figure will be saved
        filename: The file name of the figure without the png extension
    """
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    file_path = os.path.join(save_path, filename + '.png')
    fig.savefig(file_path, dpi=40)


def test_magnetic_sample_SLD(refnx_magnetic_structure):
    """
    Tests whether the spin-up and spin-down directions of the magnetic
    structure are as expected
    """
    sample = Sample(refnx_magnetic_structure)
    assert sample.structures[0][3].sld.real.value == 8
    assert sample.structures[1][3].sld.real.value == 4


def test_magnetic_sample_SLD_unpolarised(refnx_magnetic_structure):
    """
    Tests whether the unpolarised state of the magnetic
    structure are as expected
    """
    sample = Sample(refnx_magnetic_structure, polarised=False)
    assert sample.structures[0][3].sld.real.value == 6


def test_magnetic_sample_length(refnx_magnetic_structure,
                                refnx_magnetic_structure_multiple_layers):
    """
    Tests whether the amount of structures for a magnetic layer is twice as
    much as the defined amount of structures (one per spin-state).
    """
    sample_single = Sample(refnx_magnetic_structure)
    sample_double = Sample(refnx_magnetic_structure_multiple_layers)

    assert len(sample_single.structures) == 2 * len(sample_single._structures)
    assert len(sample_double.structures) == 2 * len(sample_double._structures)


def test_sample_with_multiple_bkg_order(refnx_three_solvents):
    """
    Tests whether the order of the sample backgrounds is still as expected
    when using multiple backgrounds
    """
    bkg = [2e-6, 5e-6, 1e-5]
    sample = Sample(refnx_three_solvents, bkg=bkg)
    assert bkg == sample.bkg
    assert refnx_three_solvents == sample.structures


def test_sample_with_multiple_scale_order(refnx_three_solvents):
    """
    Tests whether the order of the sample scales is still as expected
    when using multiple scales
    """
    scale = [2, 5, 1]
    sample = Sample(refnx_three_solvents, scale=scale)
    assert scale == sample.scale
    assert refnx_three_solvents == sample.structures


def test_sample_with_multiple_dq_order(refnx_three_solvents):
    """
    Tests whether the order of the sample backgrounds is still as expected
    when using multiple backgrounds
    """
    dq = [2, 5, 1]
    sample = Sample(refnx_three_solvents, dq=dq)
    assert dq == sample.dq
    assert refnx_three_solvents == sample.structures


@pytest.mark.parametrize('bkg', ([1e-6],
                                 [1e-6, 5e-6, 2e-6],
                                 [1e-6, 5e-6, 2e-6, 4e-6])
                         )
def test_sample_with_multiple_bkg_length(refnx_two_solvents, bkg):
    """
    Tests whether a ValueError is properly raised when a list of backgrounds
    is given to a sample that does not equal the amount of structures
    """
    with pytest.raises(ValueError):
        Sample(refnx_two_solvents, bkg=bkg)


@pytest.mark.parametrize('scale', ([1],
                                   [1, 5, 2],
                                   [1, 5, 2, 4])
                         )
def test_sample_with_multiple_scales_length(refnx_two_solvents, scale):
    """
    Tests whether a ValueError is properly raised when a list of scales
    is given to a sample that does not equal the amount of structures
    """
    with pytest.raises(ValueError):
        Sample(refnx_two_solvents, scale=scale)


@pytest.mark.parametrize('dq', ([1],
                                [1, 5, 2],
                                [1, 5, 2, 4])
                         )
def test_sample_with_multiple_dq_length(refnx_two_solvents, dq):
    """
    Tests whether a ValueError is properly raised when a list of dq's
    is given to a sample that does not equal the amount of structures
    """
    with pytest.raises(ValueError):
        Sample(refnx_two_solvents, dq=dq)


@pytest.mark.parametrize('label', (['1'],
                                   ['1', '2', '3'],
                                   ['1', '2', '3', '4'])
                         )
def test_sample_with_labels_length(refnx_two_solvents, label):
    """
    Tests whether a ValueError is properly raised when a list of labels
    is given to a sample that does not equal the amount of structures
    """
    with pytest.raises(ValueError):
        Sample(refnx_two_solvents, labels=label)


def test_sample_with_labels_type(refnx_two_solvents):
    """
    Tests whether a TypeError is properly raised when the labels are not given
    as a list
    """
    label = 'Structure 1'
    with pytest.raises(TypeError):
        Sample(refnx_two_solvents, labels=label)


@pytest.mark.parametrize('label', (['1', 2],
                                   [1, '2'])
                         )
def test_sample_with_labels_string_type(refnx_two_solvents, label):
    """
    Tests whether a TypeError is properly raised when any of the labels are
    not given as a string
    """
    with pytest.raises(TypeError):
        Sample(refnx_two_solvents, labels=label)


def test_angle_info(refnx_sample):
    """
    Tests whether the angle_info function correctly calculates the Fisher
    information, and outputs the same values as if the functions were called
    manually.
    """
    # Get Fisher information from tested unit
    angle_times = [(0.7, 100, 100000), (2.0, 100, 100000)]

    # Get Fisher information directly
    model = refnx_sample.get_models()[0]
    sim = SimulateReflectivity(model, angle_times)
    data = sim.simulate()
    qs, counts, models = [data[0]], [data[3]], [model]
    g = Fisher(qs, refnx_sample.params, counts, models).fisher_information
    angle_info = refnx_sample.angle_info(angle_times).fisher_information

    np.testing.assert_allclose(g, angle_info, rtol=1e-08)


@patch('hogben.models.samples.Sample._get_sld_profile')
@patch('hogben.models.samples.save_plot', side_effect=mock_save_plot)
def test_sld_profile_valid_figure(_mock_save_plot,
                                  mock_sld_profile, refnx_sample):
    """
    Tests whether the sld_profile function succesfully outputs a figure
    """
    mock_sld_profile.return_value = [([0, 10, 60, 110, 160, 210],
                                      [4, 9, -2, 9, -2, 9])]

    # Use temporary directory, so it doesn't leave any files after testing
    with tempfile.TemporaryDirectory() as temp_dir:
        refnx_sample.sld_profile(temp_dir)
        img_test = os.path.join(temp_dir, refnx_sample.name, 'sld_profile.png')
        assert os.path.isfile(img_test)


@patch('hogben.models.samples.Sample._get_reflectivity_profile')
@patch('hogben.models.samples.save_plot', side_effect=mock_save_plot)
def test_reflectivity_profile_valid_figure(_mock_save_plot,
                                           _mock_reflectivity_profile,
                                           refnx_sample):
    """
    Tests whether the reflectivity_profile function succesfully outputs a
    figure
    """
    _mock_reflectivity_profile.return_value = [([0, 0.05, 0.1, 0.15, 0.2],
                                                [1, 0.9, 0.8, 0.75, 0.8])]
    # Use temporary directory, so it doesn't leave any files after testing
    with tempfile.TemporaryDirectory() as temp_dir:
        refnx_sample.reflectivity_profile(temp_dir)
        img_test = os.path.join(temp_dir, refnx_sample.name,
                                'reflectivity_profile.png')
        assert os.path.isfile(img_test)


@patch('hogben.models.samples.save_plot', side_effect=mock_save_plot)
def test_sld_profile_length(_mock_save_plot, refnx_sample):
    """
    Tests whether _get_sld_profile() succesfully retrieves two arrays with
    equal lengths, representing an SLD profile that can be plotted in a figure
    """
    z, slds = refnx_sample._get_sld_profile()[0]
    assert len(z) == len(slds)
    assert len(z) > 0  # Make sure arrays are not empty


def test_reflectivity_profile_positive(refnx_sample):
    """
    Tests whether _get_reflectivity_profile() succesfully obtains reflectivity
    values that are all positively valued
    """
    q, r = refnx_sample._get_reflectivity_profile(0.005,
                                                  0.4,
                                                  500,
                                                  1,
                                                  1e-7,
                                                  2)[0]
    assert np.all(np.greater(r, 0.0))


def test_reflectivity_invalid_structure():
    """
    Test whether a RunTimeError is correctly given when an invalid sample
    structure is used in get_reflectivity_profile
    """
    sample = Mock(spec=None)
    with pytest.raises(TypeError):
        Sample._get_reflectivity_profile(sample, 0.005, 0.4, 500, 1, 1e-7, 2)


def test_sld_invalid_structure():
    """
    Test whether a RunTimeError is correctly given when an invalid sample
    structure is used in get_sld_profile
    """
    sample = Mock(spec=None)
    with pytest.raises(TypeError):
        Sample._get_sld_profile(sample)[0]


def test_vary_structure_invalid_structure():
    """
    Test whether a RunTimeError is correctly given when an invalid sample
    structure is used in _vary_structure
    """
    structure = Mock(spec=None)
    with pytest.raises(TypeError):
        Sample._vary_structure(structure)


def test_reflectivity_profile_length(refnx_sample):
    """
    Tests whether _get_reflectivity_profile() succesfully retrieves two arrays
    with equal lengths, representing a reflectivity profile that can be
    plotted in a figure.
    """
    q, r = refnx_sample._get_reflectivity_profile(0.005,
                                                  0.4,
                                                  500,
                                                  1,
                                                  1e-7,
                                                  2)[0]
    assert len(q) == len(r)
    assert len(q) > 0  # Make sure array is not empty


@patch('hogben.models.samples.save_plot', side_effect=mock_save_plot)
def test_main_function(_mock_save_plot):
    """
    Tests whether the main function runs properly and creates a figure for
    all defined model types.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        save_path = os.path.join(temp_dir, 'results')
        samples.run_main(save_path)

        for subfolder in os.listdir(save_path):
            reflectivity_profile = os.path.join(save_path, subfolder,
                                                'reflectivity_profile.png')
            sld_profile = os.path.join(save_path, subfolder, 'sld_profile.png')
            assert os.path.isfile(reflectivity_profile)
            assert os.path.isfile(sld_profile)
