"""Tests for the simulation module"""

import pytest

import numpy as np
from refnx.reflect import SLD, ReflectModel

from hogben.simulate import SimulateReflectivity


@pytest.fixture(scope='module')
def refnx_structure():
    """Defines a refnx structure describing a simple sample."""
    air = SLD(0, name='Air')
    layer1 = SLD(4, name='Layer 1')(thick=100, rough=2)
    layer2 = SLD(8, name='Layer 2')(thick=150, rough=2)
    substrate = SLD(2.047, name='Substrate')(thick=0, rough=2)

    sample_1 = air | layer1 | layer2 | substrate
    return sample_1


@pytest.fixture(scope='module')
def refnx_model(refnx_structure):
    """A refnx model for use in testing"""
    model = ReflectModel(refnx_structure)
    model.bkg = 1e-6
    model.dq = 2
    model.scale = 1.0
    return model


class TestSimulate:
    """A class to test the SimulateReflectivity methods"""
    angle_times = [(0.3, 100, 1000)]  # (Angle, Points, Time)
    instrument = 'OFFSPEC'

    def test_data_streaming(self, refnx_model):
        """
        Tests that without an input for the datafile,
        the correct one is picked up
        """
        sim = SimulateReflectivity(refnx_model,
                                   self.angle_times,
                                   self.instrument)
        simulated_datapoints = sim.simulate()
        zeros = np.zeros_like(simulated_datapoints)
        assert np.all(zeros <= simulated_datapoints[3])  # Counts

        # Check that the default instrument also works
        sim = SimulateReflectivity(refnx_model,
                                   self.angle_times)
        simulated_datapoints = sim.simulate()
        zeros = np.zeros_like(simulated_datapoints)
        assert np.all(zeros <= simulated_datapoints[3])  # Counts

    def test_incident_flux_data(self):
        """
        Tests that the `_incident_flux_data` function correctly returns
        count data when instruments are loaded, and raises errors when
        given an incorrect path
        """
        # Test that non-polarised instruments work and have counts in
        for instrument in ['OFFSPEC', 'SURF', 'INTER', 'POLREF']:
            sim_no_pol = SimulateReflectivity(None,
                                              angle_times=self.angle_times,
                                              inst_or_path=instrument)

            assert (sim_no_pol._incident_flux_data(polarised=False).shape[1]
                    == 2)
            counts = sim_no_pol._incident_flux_data(polarised=False)[:, 1]
            assert (np.sum(counts) > 100)

        # Test that polarised instruments work and have counts in
        for instrument in ['OFFSPEC', 'POLREF']:
            sim_pol = SimulateReflectivity(None,
                                           angle_times=self.angle_times,
                                           inst_or_path=instrument)

            assert (sim_pol._incident_flux_data(polarised=True).shape[1]
                    == 2)
            assert (np.sum(sim_pol._incident_flux_data(polarised=True)[:, 1])
                    > 100)

        # Test that a non-existing path raises an error
        with pytest.raises(FileNotFoundError):
            sim_wrong_path = SimulateReflectivity(None,
                                                  angle_times=self.angle_times,
                                                  inst_or_path='doesnt_exist')

            sim_wrong_path._incident_flux_data()

        # Test that a blank instrument raises an error
        with pytest.raises(FileNotFoundError):
            sim_no_path = SimulateReflectivity(None,
                                               angle_times=self.angle_times,
                                               inst_or_path='')

            sim_no_path._incident_flux_data()

        # Test that a valid non-polarised instrument can't be used
        # for a polarised simulation
        with pytest.raises(FileNotFoundError):
            sim_not_in_pol = SimulateReflectivity(None,
                                                  angle_times=self.angle_times,
                                                  inst_or_path='SURF')

            sim_not_in_pol._incident_flux_data(polarised=True)

    def test_number_of_points(self, refnx_model):
        """
        Tests that the number of points generated in the
        simulation is what we define (even when many are zero)
        """
        NPOINTS = 100
        angle_times = [(0.3, NPOINTS, 10000)]
        sim = SimulateReflectivity(refnx_model,
                                   angle_times=angle_times).simulate()
        assert len(sim[0]) == NPOINTS
        angle_times = [(0.3, NPOINTS, 5)]
        sim = SimulateReflectivity(refnx_model,
                                   angle_times=angle_times).simulate()
        assert len(sim[0]) == NPOINTS
        angle_times = [(0.3, NPOINTS, 5), (2.0, NPOINTS, 50)]
        sim = SimulateReflectivity(refnx_model,
                                   angle_times=angle_times).simulate()
        assert len(sim[0]) == NPOINTS * 2

    def test_reflectivity(self, refnx_model):
        """
        Checks that a refnx model reflectivity generated through
        `hogben.reflectivity` is always greater than zero.
        """
        sim = SimulateReflectivity(refnx_model)
        ideal_reflectivity = sim.reflectivity(np.linspace(0.001, 0.3, 200))

        np.testing.assert_array_less(np.zeros(len(ideal_reflectivity)),
                                     ideal_reflectivity)

    @pytest.mark.parametrize('polarised', (True, False))
    def test_run_experiment(self, refnx_model, polarised):
        """
        Checks the output of _run_experiment gives the right
        length of outputs
        """
        sim = SimulateReflectivity(refnx_model,
                                   self.angle_times,
                                   self.instrument)

        q_binned, r_noisy, r_error, counts_incident = (
            sim._run_experiment(*self.angle_times[0], polarised=polarised))

        for item in [q_binned, r_noisy, r_error, counts_incident]:
            assert len(item) == self.angle_times[0][1]

        np.testing.assert_array_less(np.zeros_like(q_binned), q_binned)

    @pytest.mark.parametrize('polarised', (True, False))
    def test_simulate_multiple_angles(self, refnx_model, polarised):
        """
        Checks that simulated reflectivity data points and simulated neutron
        counts generated through `hogben.simulate` are always greater than
        zero (given a long count time).
        """
        angle_times = [(0.3, 100, 1000), (2.3, 150, 1000)]

        sim = SimulateReflectivity(refnx_model, angle_times, self.instrument)
        q_binned, r_noisy, r_error, counts_incident = (
            sim.simulate(polarised=polarised))

        for item in [q_binned, r_noisy, r_error, counts_incident]:
            assert len(item) == sum(condition[1] for condition in angle_times)

        np.testing.assert_array_less(np.zeros_like(q_binned), q_binned)
