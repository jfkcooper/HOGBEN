"""Contains the class to define magnetic samples"""

import os

import matplotlib.pyplot as plt

import numpy as np

import refl1d.material
import refl1d.model
import refl1d.probe
import refl1d.experiment
import refl1d.magnetism

import bumps.parameter
import bumps.fitproblem

from hogben.models.base import BaseSample, VariableUnderlayer
from hogben.utils import Fisher, Sampler, save_plot

plt.rcParams['figure.figsize'] = (9, 7)
plt.rcParams['figure.dpi'] = 600


class SampleYIG(BaseSample, VariableUnderlayer):
    """Defines a magnetic model describing yttrium iron garnet (YIG) film
       grown on a yttrium aluminium garnet (YAG) substrate.

    Attributes:
        name (str): name of the magnetic sample.
        data_path (str): path to directory containing measured data.
        labels (list): label for each measured data set.
        mag_angle (float): magnetic angle in degrees,
        scale (float): experimental scale factor
        bkg (float): level of instrument background noise.
        dq (float): instrument resolution.
        pt_sld (bumps.parameter.Parameter): platinum layer SLD.
        pt_thick (bumps.parameter.Parameter): platinum layer thickness.
        pt_rough (bumps.parameter.Parameter): air/platinum roughness.
        pt_mag (bumps.parameter.Parameter): platinum layer magnetic SLD.
        intermediary_sld (bumps.parameter.Parameter): intermediary layer SLD.
        intermediary_thick (bumps.parameter.Parameter): intermediary layer
                                                        thickness.
        intermediary_rough (bumps.parameter.Parameter): platinum/intermediary
                                                        roughness.
        yig_sld (bumps.parameter.Parameter): YIG layer SLD.
        yig_thick (bumps.parameter.Parameter): YIG layer thickness.
        yig_rough (bumps.parameter.Parameter): intermediary/YIG roughness.
        yig_mag (bumps.parameter.Parameter): yig layer magnetic SLD.
        yag_sld (bumps.parameter.Parameter): YAG substrate SLD.
        yag_rough (bumps.parameter.Parameter): YIG/YAG roughness.
        params (list): parameters of the model.
        structure (refl1d.model.Stack): Refl1D representation of sample.
        experiment (refl1d.experiment.Experiment): fittable experiment for
                                                   sample.

    """

    def __init__(self):
        """Initialise YIG sample and set the sample properties"""
        self.name = 'YIG_sample'
        self.data_path = os.path.join(
            os.path.dirname(__file__), '..', 'data', 'YIG_sample'
        )

        self.labels = ['Up', 'Down']
        self.mag_angle = 90
        self.scale = 1.025
        self.bkg = 4e-7
        self.dq = 2.8

        # Define the parameters of the model.
        self.pt_sld = bumps.parameter.Parameter(5.646, name='Pt SLD')
        self.pt_thick = bumps.parameter.Parameter(21.08, name='Pt Thickness')
        self.pt_rough = bumps.parameter.Parameter(
            8.211, name='Air/Pt Roughness'
        )
        self.pt_mag = bumps.parameter.Parameter(0, name='Pt Magnetic SLD')

        self.intermediary_sld = bumps.parameter.Parameter(
            4.678, name='Intermediary SLD'
        )
        self.intermediary_thick = bumps.parameter.Parameter(
            19.67, name='Intermediary Thickness'
        )
        self.intermediary_rough = bumps.parameter.Parameter(
            2, name='Pt/Intermediary Roughness'
        )

        self.yig_sld = bumps.parameter.Parameter(5.836, name='YIG SLD')
        self.yig_thick = bumps.parameter.Parameter(713.8, name='YIG Thickness')
        self.yig_rough = bumps.parameter.Parameter(
            13.55, name='Intermediary/YIG Roughness'
        )
        self.yig_mag = bumps.parameter.Parameter(
            0.349, name='YIG Magnetic SLD'
        )

        self.yag_sld = bumps.parameter.Parameter(5.304, name='YAG SLD')
        self.yag_rough = bumps.parameter.Parameter(
            30, name='YIG/YAG Roughness'
        )

        """
        self.params = [self.pt_sld,
                       self.pt_thick,
                       self.pt_rough,
                       self.pt_mag,
                       self.intermediary_sld,
                       self.intermediary_thick,
                       self.intermediary_rough,
                       self.yig_sld,
                       self.yig_thick,
                       self.yig_rough,
                       self.yig_mag,
                       self.yag_sld,
                       self.yag_rough]
        """

        self.params = [self.pt_mag]

        self.pt_sld.range(5, 6)
        self.pt_thick.range(2, 30)
        self.pt_rough.range(0, 9)
        self.pt_mag.range(0, 0.1)

        self.intermediary_sld.range(4.5, 5.5)
        self.intermediary_thick.range(0, 25)
        self.intermediary_rough.range(2, 10)

        self.yig_sld.range(5, 6)
        self.yig_thick.range(100, 900)
        self.yig_rough.range(0, 70)
        self.yig_mag.range(0.2, 0.5)

        self.yag_sld.range(4.5, 5.5)
        self.yag_rough.range(20, 30)

        # Load the experimentally-measured data for the sample.
        self.__create_experiment()

    def __create_experiment(self):
        """Creates an experiment corresponding to the measured data."""
        # Load the "up" and "down" spin state data sets.
        file_path_up = os.path.join(self.data_path, 'YIG_up.dat')
        file_path_down = os.path.join(self.data_path, 'YIG_down.dat')

        pp = refl1d.probe.load4(
            file_path_up, sep='\t', intensity=self.scale, background=self.bkg
        )

        mm = refl1d.probe.load4(
            file_path_down, sep='\t', intensity=self.scale, background=self.bkg
        )

        # Set the resolution to be constant dQ/q.
        self.__set_dq(pp)
        self.__set_dq(mm)

        # Combine the up and down probes into a single polarised Q probe.
        probe = refl1d.probe.PolarizedQProbe(xs=(pp, None, None, mm), name='')

        # Define the experiment using the probe and sample structure.
        self.structure = self.using_conditions()
        self.experiment = refl1d.experiment.Experiment(
            sample=self.structure, probe=probe
        )

    def using_conditions(self, yig_thick=None, pt_thick=None):
        """Creates a structure representing the YIG sample measured using
           given measurement conditions.

        Args:
            yig_thick (float): thickness of YIG layer to use.
            pt_thick (float): thickness of platinum layer to use.

        Returns:
            refl1d.model.Stack: structure defined using given conditions.

        """
        # If not given a YIG thickness, use the fitted value.
        if yig_thick is None:
            yig_thick = self.yig_thick

        # If not given a platinum thickness, use the fitted value.
        if pt_thick is None:
            pt_thick = 0
        else:
            # Subtract existing thickness and check result is non-negative.
            pt_thick -= self.pt_thick.value
            assert pt_thick >= 0

        # 1 uB/atom = 1.638
        pt_magnetism = refl1d.magnetism.Magnetism(
            rhoM=self.pt_mag, thetaM=self.mag_angle
        )
        yig_magnetism = refl1d.magnetism.Magnetism(
            rhoM=self.yig_mag, thetaM=self.mag_angle
        )

        # Define the sample structure.
        air = refl1d.material.SLD(rho=0, name='Air')
        pt = refl1d.material.SLD(rho=self.pt_sld, name='Pt')(
            self.pt_thick, self.pt_rough, magnetism=pt_magnetism
        )
        intermediary = refl1d.material.SLD(
            rho=self.intermediary_sld, name='Intermediary'
        )(self.intermediary_thick, self.intermediary_rough)
        yig = refl1d.material.SLD(rho=self.yig_sld, name='YIG')(
            yig_thick, self.yig_rough, magnetism=yig_magnetism
        )
        yag = refl1d.material.SLD(rho=self.yag_sld, name='Substrate')(
            0, self.yag_rough
        )

        # Define the added platinum thickness, if applicable.
        pt_extra = refl1d.material.SLD(rho=self.pt_sld, name='Pt Extra')(
            pt_thick, 0
        )

        # Add the extra platinum layer if requested.
        if pt_thick == 0:
            return yag | yig | intermediary | pt | air
        else:
            return yag | yig | intermediary | pt_extra | pt | air

    def angle_info(self, angle_times, contrasts=None):
        """Calculates the Fisher information matrix for the YIG sample
           measured over a number of angles.

        Args:
            angle_times (list): points and times for each angle to simulate.
            contrasts (list): not applicable.

        Returns:
            numpy.ndarray: Fisher information matrix.

        """
        # Simulate the polarised experiment.
        models, datasets = [], []

        """simulate_magnetic(self.structure, angle_times,
                                             scale=1, bkg=5e-7, dq=2,
                                             pp=True, pm=False,
                                             mp=False, mm=True)
        """
        # Calculate the Fisher information matrix.
        qs = [data[:, 0] for data in datasets]
        counts = [data[:, 3] for data in datasets]
        return Fisher(qs, self.params, counts, models)

    def underlayer_info(self, angle_times, yig_thick, pt_thick):
        """Calculates the Fisher information matrix for the YIG sample
           with given YIG and platinum layer thicknesses.

        Args:
            angle_times (list): points and times for each angle to simulate.
            yig_thick (float): YIG layer thickness to use.
            pt_thick (float): platinum layer thickness to use.

        Returns:
            Fisher: Fisher information object.

        """
        # Create a structure with the given YIG and Pt thicknesses.
        # structure = self.using_conditions(yig_thick, pt_thick)

        # Simulate a polarised measurement of the structure.
        models, datasets = [], []
        """simulate_magnetic(structure, angle_times,
                                             scale=1, bkg=5e-7, dq=2,
                                             pp=True, pm=False,
                                             mp=False, mm=True)
        """
        # Calculate the Fisher information matrix.
        qs = [data[:, 0] for data in datasets]
        counts = [data[:, 3] for data in datasets]
        return Fisher(qs, self.params, counts, models)

    def __set_dq(self, probe):
        """Sets the resolution of a given `probe` to be constant dQ/Q.

        Args:
            probe (refl1d.probe.QProbe): probe to set the resolution for.

        """
        # Transform the resolution from refnx to Refl1D format.
        dq = self.dq / (100 * np.sqrt(8 * np.log(2)))

        q_array = probe.Q

        # Calculate the dQ array and update the QProbe.
        dq_array = probe.Q * dq
        probe.dQ = dq_array

        # Adjust probe calculation for constant resolution.
        argmin, argmax = np.argmin(q_array), np.argmax(q_array)
        probe.calc_Qo = np.linspace(q_array[argmin] - 3.5 * dq_array[argmin],
                                    q_array[argmax] + 3.5 * dq_array[argmax],
                                    21 * len(q_array))

    def sld_profile(self, save_path):
        """Plots the SLD profile of the sample.

        Args:
            save_path (str): path to directory to save SLD profile to.

        """
        # Get the SLD profile values.
        z, slds, _, slds_mag, _ = self.experiment.magnetic_smooth_profile()

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)

        # Plot the SLD profile.
        ax.plot(z, slds, color='black', label='SLD')
        ax.plot(z, slds_mag, color='red', label='Magnetic SLD')

        x_label = r'$\mathregular{Distance\ (\AA)}$'
        y_label = r'$\mathregular{SLD\ (10^{-6} \AA^{-2})}$'

        ax.set_xlabel(x_label, fontsize=11, weight='bold')
        ax.set_ylabel(y_label, fontsize=11, weight='bold')
        ax.legend()

        # Save the plot.
        save_path = os.path.join(save_path, self.name)
        save_plot(fig, save_path, 'sld_profile')

    def reflectivity_profile(self, save_path):
        """Plots the reflectivity profile of the YIG sample.

        Args:
            save_path (str): path to directory to save profile to.

        """
        fig = plt.figure(figsize=(6, 7))
        ax = fig.add_subplot(111)

        colours = ['b', 'g']
        count = 0
        probes = self.experiment.probe.xs
        reflectance = self.experiment.reflectivity()

        # Plot the data for the "up" and "down" spin states.
        for probe, qr in zip(probes, reflectance):
            if qr is not None:
                # Plot the experimentally-measured data.
                ax.errorbar(probe.Q, probe.R, probe.dR, marker='o', ms=2,
                            lw=0, elinewidth=0.5, capsize=0.5,
                            color=colours[count],
                            label=self.labels[count] + ' Data')

                # Plot the model reflectivity.
                ax.plot(probe.Q, qr[1], zorder=20,
                        color=colours[count],
                        label=self.labels[count] + ' Fitted')

                count += 1

        x_label = r'$\mathregular{Q\ (Å^{-1})}$'
        y_label = 'Reflectivity (arb.)'

        ax.set_xlabel(x_label, fontsize=11, weight='bold')
        ax.set_ylabel(y_label, fontsize=11, weight='bold')
        ax.set_yscale('log')
        ax.legend(loc='lower left')

        # Save the plot.
        save_path = os.path.join(save_path, self.name)
        save_plot(fig, save_path, 'reflectivity_profile')

    def nested_sampling(self, angle_times, save_path, filename, dynamic=False):
        """Runs nested sampling on simulated data of the YIG sample.

        Args:
            angle_times (list): points and times for each angle to simulate.
            save_path (str): path to directory to save corner plot to.
            filename (str): file name to use when saving corner plot.
            dynamic (bool): whether to use static or dynamic nested sampling.

        """
        # Define the object to sample.
        objective = bumps.fitproblem.FitProblem(self.experiment)

        # Sample the objective using nested sampling.
        sampler = Sampler(objective)
        fig = sampler.sample(dynamic=dynamic)

        # Save the sampling corner plot.
        save_path = os.path.join(save_path, self.name)
        save_plot(fig, save_path, 'nested_sampling_' + filename)


if __name__ == '__main__':
    print('Magnetism as it was coded currently does not work in this version')
    """save_path = '../results'

    # Save the SLD and reflectivity profiles of the YIG sample.
    yig_sample = SampleYIG()
    yig_sample.sld_profile(save_path)
    yig_sample.reflectivity_profile(save_path)

    # Close the plots.
    plt.close('all')"""
