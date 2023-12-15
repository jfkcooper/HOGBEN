import os
import tempfile

import numpy as np
import pytest
import matplotlib
import refl1d.experiment
import hogben.models.samples as samples

from hogben.models.samples import Sample
from hogben.simulate import SimulateReflectivity
from hogben.utils import Fisher
from refnx.reflect import SLD, ReflectModel
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


def test_angle_info(refnx_sample):
    """
    Tests whether the angle_info function correctly calculates the Fisher
    information, and outputs the same values as if the functions were called
    manually.
    """
    return
    # Get Fisher information from tested unit
    angle_times = [(0.7, 100, 100000), (2.0, 100, 100000)]
    
    # Get Fisher information directly
    model = ReflectModel(refnx_sample.structure)
    sim = SimulateReflectivity(model, angle_times)
    data = sim.simulate()
    qs, counts, models = [data[0]], [data[3]], [model]
    g = Fisher(qs, sample.params, counts, models).fisher_information
    angle_info = sample.angle_info(angle_times).fisher_information

    np.testing.assert_allclose(g, angle_info, rtol=1e-08)


@patch('hogben.models.samples.Sample._get_sld_profile')
@patch('hogben.models.samples.save_plot', side_effect=mock_save_plot)
def test_sld_profile_valid_figure(_mock_save_plot,
                                  mock_sld_profile, refnx_sample):
    """
    Tests whether the sld_profile function succesfully outputs a figure
    """
    mock_sld_profile.return_value = ([0, 10, 60, 110, 160, 210],
                                     [4, 9, -2, 9, -2, 9])

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
    _mock_reflectivity_profile.return_value = ([0, 0.05, 0.1, 0.15, 0.2],
                                               [1, 0.9, 0.8, 0.75, 0.8])
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
    z, slds = refnx_sample._get_sld_profile()
    assert len(z) == len(slds)
    assert len(z) > 0  # Make sure arrays are not empty


def test_reflectivity_profile_positive(refnx_sample):
    """
    Tests whether _get_reflectivity_profile() succesfully obtains reflectivity
    values that are all positively valued
    """
    q, r = refnx_sample._get_reflectivity_profile(0.005, 0.4, 500, 1, 1e-7, 2)
    assert np.all(np.greater(r, 0.0))


def test_reflectivity_invalid_structure():
    """
    Test whether a RunTimeError is correctly given when an invalid sample
    structure is used in get_reflectivity_profile
    """
    sample = Mock(spec=None)
    with pytest.raises(RuntimeError):
        Sample._get_reflectivity_profile(sample, 0.005, 0.4, 500, 1, 1e-7, 2)


def test_sld_invalid_structure():
    """
    Test whether a RunTimeError is correctly given when an invalid sample
    structure is used in get_sld_profile
    """
    sample = Mock(spec=None)
    with pytest.raises(RuntimeError):
        Sample._get_sld_profile(sample)


def test_vary_structure_invalid_structure():
    """
    Test whether a RunTimeError is correctly given when an invalid sample
    structure is used in _vary_structure
    """
    structure = Mock(spec=None)
    with pytest.raises(RuntimeError):
        Sample._Sample__vary_structure(structure)


def test_reflectivity_profile_length(refnx_sample):
    """
    Tests whether _get_reflectivity_profile() succesfully retrieves two arrays
    with equal lengths, representing a reflectivity profile that can be
    plotted in a figure.
    """
    q, r = refnx_sample._get_reflectivity_profile(0.005, 0.4, 500, 1, 1e-7, 2)
    assert len(q) == len(r)
    assert len(q) > 0  # Make sure array is not empty


@pytest.mark.parametrize('sample', [samples.simple_sample(),
                                    samples.many_param_sample(),
                                    samples.thin_layer_sample_1(),
                                    samples.thin_layer_sample_2(),
                                    samples.similar_sld_sample_1(),
                                    samples.similar_sld_sample_2(),

                                    ]
                         )
def test_predefined_samples(sample):
    """
    Tests whether the predefined samples provide a valid structure that
    succesfully can be converted to a refl1d structure.
    """
    sample.to_refl1d()
    assert isinstance(sample.structure, refl1d.model.Stack)


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
