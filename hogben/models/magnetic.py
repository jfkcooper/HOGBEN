"""Contains the class to define a native refnx YIG magnetic sample."""

import os

import matplotlib.pyplot as plt

import refnx.analysis
import refnx.dataset
import refnx.reflect
from refnx.analysis import Parameter
from refnx.reflect import ReflectModel, SLD
from refnx.reflect.structure import MagneticSlab, Slab

from hogben.models.base import BaseSample, VariableUnderlayer
from hogben.utils import Fisher, Sampler, save_plot

plt.rcParams['figure.figsize'] = (9, 7)
plt.rcParams['figure.dpi'] = 600


class SampleYIG(BaseSample, VariableUnderlayer):
    """Defines a magnetic YIG sample composed of a platinum layer,
    an intermediary layer, a magnetic YIG layer, and a YAG substrate.

    Attributes:
        name (str): Name of the sample.
        data_path (str): Path to the measured YIG data files.
        labels (list): Spin-state labels for the measured data.
        mag_angle (float): Magnetic field angle in degrees.
        scale (list): Scale factors applied to the up/down measurements.
        bkg (list): Background values for each spin-state dataset.
        dq (list): Resolution values used for the reflectivity models.
        pt_sld (Parameter): SLD of the Pt layer.
        pt_thick (Parameter): Thickness of the Pt layer.
        pt_rough (Parameter): Air/Pt roughness.
        pt_mag (Parameter): Magnetic SLD of the Pt layer.
        intermediary_sld (Parameter): SLD of the intermediary layer.
        intermediary_thick (Parameter): Thickness of the intermediary layer.
        intermediary_rough (Parameter): Pt/intermediary roughness.
        yig_sld (Parameter): Non-magnetic SLD of the YIG layer.
        yig_thick (Parameter): Thickness of the YIG layer.
        yig_rough (Parameter): Intermediary/YIG roughness.
        yig_mag (Parameter): Magnetic SLD of the YIG layer.
        yag_sld (Parameter): SLD of the YAG substrate.
        yag_rough (Parameter): YIG/YAG roughness.
        params (list): Parameters that are varied during fitting.
        structures (list): Native `refnx` structures for the sample.
    """

    def __init__(self):
        """Initialise the YIG sample and define its parameters."""
        super().__init__()

        self.name = 'YIG_sample'
        self.data_path = os.path.join(
            os.path.dirname(__file__), '..', 'data', 'YIG_sample'
        )
        self.labels = ['Up', 'Down']
        self.mag_angle = 90.0
        self.scale = [1.025, 1.025]
        self.bkg = [4e-7, 4e-7]
        self.dq = [2.8, 2.8]

        self.pt_sld = Parameter(5.646, name='Pt SLD', bounds=(5, 6))
        self.pt_thick = Parameter(21.08, name='Pt Thickness', bounds=(2, 30))
        self.pt_rough = Parameter(8.211, name='Air/Pt Roughness',
                                  bounds=(0, 9))
        self.pt_mag = Parameter(0.0, name='Pt Magnetic SLD',
                                bounds=(0, 0.1), vary=True)

        self.intermediary_sld = Parameter(4.678, name='Intermediary SLD',
                                          bounds=(4.5, 5.5))
        self.intermediary_thick = Parameter(19.67,
                                            name='Intermediary Thickness',
                                            bounds=(0, 25))
        self.intermediary_rough = Parameter(
            2, name='Pt/Intermediary Roughness', bounds=(2, 10)
        )

        self.yig_sld = Parameter(5.836, name='YIG SLD', bounds=(5, 6))
        self.yig_thick = Parameter(713.8, name='YIG Thickness',
                                   bounds=(100, 900))
        self.yig_rough = Parameter(13.55, name='Intermediary/YIG Roughness',
                                   bounds=(0, 70))
        self.yig_mag = Parameter(0.349, name='YIG Magnetic SLD',
                                 bounds=(0.2, 0.5), vary=True)

        self.yag_sld = Parameter(5.304, name='YAG SLD', bounds=(4.5, 5.5))
        self.yag_rough = Parameter(30, name='YIG/YAG Roughness',
                                   bounds=(20, 30))

        self.params = [self.pt_mag, self.yig_mag]
        self._structures = [self.using_conditions()]
        self.polarised = True

        self._create_objectives()

    def _create_objectives(self):
        """Create refnx objectives for the measured YIG spin states."""
        file_paths = [
            os.path.join(self.data_path, 'YIG_up.dat'),
            os.path.join(self.data_path, 'YIG_down.dat'),
        ]

        self.objectives = []
        for structure, file_path, scale, bkg in zip(
            self.structures,
            file_paths,
            self.scale,
            self.bkg,
        ):
            model = ReflectModel(
                structure, scale=scale, bkg=bkg, dq=self.dq[0]
            )
            data = refnx.dataset.ReflectDataset(file_path)
            self.objectives.append(refnx.analysis.Objective(model, data))

    def using_conditions(self, yig_thick=None, pt_thick=None):
        """Create a native refnx structure for the YIG sample.

        Args:
            yig_thick (float | None): YIG thickness to use.
            pt_thick (float | None): total Pt thickness to use.

        Returns:
            refnx.reflect.Structure: YIG sample structure.
        """
        yig_thick = self.yig_thick.value if yig_thick is None else yig_thick
        pt_total = self.pt_thick.value if pt_thick is None else float(pt_thick)

        if pt_total < self.pt_thick.value:
            raise ValueError(
                'Pt thickness cannot be smaller than the base Pt layer'
            )

        pt_extra = pt_total - self.pt_thick.value

        air = SLD(0, name='Air')
        pt = MagneticSlab(
            thick=self.pt_thick.value,
            sld=self.pt_sld,
            rough=self.pt_rough,
            rhoM=self.pt_mag,
            thetaM=self.mag_angle,
            name='Pt',
        )
        intermediary = Slab(
            self.intermediary_thick,
            self.intermediary_sld,
            self.intermediary_rough,
            name='Intermediary',
        )
        yig = MagneticSlab(
            thick=yig_thick,
            sld=self.yig_sld,
            rough=self.yig_rough,
            rhoM=self.yig_mag,
            thetaM=self.mag_angle,
            name='YIG',
        )
        yag = Slab(0, self.yag_sld, self.yag_rough, name='Substrate')

        structure = air | pt | intermediary | yig | yag
        if pt_extra > 0:
            pt_extra_layer = Slab(pt_extra, self.pt_sld, 0, name='Pt Extra')
            structure = air | pt | pt_extra_layer | intermediary | yig | yag

        return structure

    def angle_info(self, angle_times, contrasts=None, inst_or_path='OFFSPEC'):
        """Calculate the Fisher information for the measured YIG geometry."""
        return Fisher.from_sample(self, angle_times, inst_or_path)

    def underlayer_info(self, angle_times, yig_thick, pt_thick):
        """Calculate the Fisher information for a given YIG/Pt geometry."""
        temp_structure = self.using_conditions(
            yig_thick=yig_thick,
            pt_thick=pt_thick,
        )

        from hogben.models.samples import Sample

        temp_sample = Sample(
            temp_structure,
            scale=1.025,
            bkg=4e-7,
            dq=2.8,
            polarised=True,
        )
        return Fisher.from_sample(
            temp_sample, angle_times, inst_or_path='OFFSPEC'
        )

    def sld_profile(self, save_path):
        """Plot the SLD profile for the YIG sample."""
        z_up, slds_up = self.structures[0].sld_profile()
        z_down, slds_down = self.structures[1].sld_profile()

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)

        ax.plot(z_up, slds_up, color='black', label='Spin-up SLD')
        ax.plot(z_down, slds_down, color='red', label='Spin-down SLD')

        ax.set_xlabel(r'$\mathregular{Distance\ (\AA)}$', fontsize=11,
                      weight='bold')
        ax.set_ylabel(r'$\mathregular{SLD\ (10^{-6} \AA^{-2})}$',
                      fontsize=11, weight='bold')
        ax.legend()

        save_path = os.path.join(save_path, self.name)
        save_plot(fig, save_path, 'sld_profile')

    def reflectivity_profile(self, save_path):
        """Plot the reflectivity profile for the measured YIG sample."""
        fig = plt.figure(figsize=(6, 7))
        ax = fig.add_subplot(111)

        colours = ['b', 'g']
        for i, objective in enumerate(self.objectives):
            q = objective.data.x
            r = objective.data.y
            model = objective.model(q)

            ax.errorbar(q, r, yerr=objective.data.yerr,
                        marker='o', ms=2, lw=0, elinewidth=0.5,
                        capsize=0.5, color=colours[i],
                        label=self.labels[i] + ' Data')
            ax.plot(q, model, zorder=20, color=colours[i],
                    label=self.labels[i] + ' Fitted')

        ax.set_xlabel(r'$\mathregular{Q\ (Å^{-1})}$', fontsize=11,
                      weight='bold')
        ax.set_ylabel('Reflectivity (arb.)', fontsize=11, weight='bold')
        ax.set_yscale('log')
        ax.legend(loc='lower left')

        save_path = os.path.join(save_path, self.name)
        save_plot(fig, save_path, 'reflectivity_profile')

    def nested_sampling(self, angle_times, save_path, filename, dynamic=False):
        """Run nested sampling on the measured YIG sample."""
        objectives = []
        for structure in self.structures:
            model = ReflectModel(structure)
            data = refnx.dataset.ReflectDataset(
                os.path.join(
                    self.data_path,
                    f'YIG_{self.labels[0].lower()}.dat'
                )
            )
            objectives.append(refnx.analysis.Objective(model, data))

        global_objective = refnx.analysis.GlobalObjective(objectives)
        global_objective.varying_parameters = lambda: self.params

        sampler = Sampler(global_objective)
        fig = sampler.sample(dynamic=dynamic)

        save_path = os.path.join(save_path, self.name)
        save_plot(fig, save_path, 'nested_sampling_' + filename)