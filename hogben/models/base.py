"""Base classes for different sample types """

import os
from abc import ABC, abstractmethod
from typing import Optional

import matplotlib.pyplot as plt

import numpy as np

import refnx.dataset
import refnx.reflect
import refnx.analysis
from refnx.reflect import ReflectModel, PolarisedReflectModel
from refnx.reflect.structure import Slab, MagneticSlab
from refnx._lib import flatten

from hogben.simulate import SimulateReflectivity
from hogben.utils import Fisher, Sampler, save_plot
from itertools import repeat


plt.rcParams['figure.figsize'] = (9, 7)
plt.rcParams['figure.dpi'] = 600


class VariableAngle(ABC):
    """Abstract class representing whether the measurement angle of a sample
       can be varied."""

    @abstractmethod
    def angle_info(self):
        """Calculates the Fisher information matrix for a sample measured
        over a number of angles."""
        pass


class VariableContrast(ABC):
    """Abstract class representing whether the contrast of a sample
       dan be varied."""

    @abstractmethod
    def contrast_info(self):
        """Calculates the Fisher information matrix for a sample with contrasts
           measured over a number of angles."""
        pass


class VariableUnderlayer(ABC):
    """Abstract class representing whether the underlayer(s) of a sample
       can be varied."""

    @abstractmethod
    def underlayer_info(self):
        """Calculates the Fisher information matrix for a sample with
        underlayers, and contrasts measured over a number of angles."""
        pass


class BaseSample(VariableAngle):
    """Abstract class representing a "standard" neutron reflectometry sample
    defined by a series of contiguous layers."""

    def __init__(self):
        """Initialise the sample and define class attributes"""
        self._structures = []
        self._bkg = None
        self._dq = None
        self._scale = None
        self.polarised = self.is_magnetic()

    def get_structures(self) -> list:
        """
        Get a list of the possible sample structures.
        """
        spin_structures = []
        if self.polarised:
            for structure in self._structures:
                up_structure = structure.copy()
                down_structure = structure.copy()
                for i, layer in enumerate(structure):
                    if isinstance(layer, MagneticSlab):
                        up_structure[i] = self._spin_structure(layer, 'up')
                        down_structure[i] = self._spin_structure(layer, 'down')
                if self.is_magnetic():
                    spin_structures.extend([up_structure, down_structure])
                else:
                    spin_structures.extend([structure.copy()])
            return spin_structures
        return self._structures

    @property
    def structures(self):
        """Return the structures that belong to the sample"""
        return self.get_structures()

    @structures.setter
    def structures(self, structures: list):
        self._structures = structures

    def is_magnetic(self) -> bool:
        """Checks whether the sample contains at least one magnetic layer"""
        for structure in self._structures:
            if getattr(structure, 'is_magnetic', False):
                return True
        return False

    def get_param_by_attribute(self, attr: str) -> list:
        """
        Get all parameters defined in the sample model that have a given
        attribute. Returns a list with all parameters with this attribute,
        e.g. `attr='vary'` returns all varying parameters.

        Args:
            attr (str): The attribute to filter for

        Returns:
            list: A list of all parameters with the given attribute
        """
        params = []
        for model in self.get_models():
            for p in flatten(model.parameters):
                if hasattr(p, attr) and getattr(p, attr):
                    params.append(p)
                    continue
                # Get parameters that are coupled to model attributes as
                # dependencies:
                if p._deps:
                    params.extend([_p for _p in p.dependencies() if
                                   hasattr(_p, attr) and getattr(_p, attr)])
        return list(set(params))

    def get_models(self) -> list:
        """
        Generates a refnx `ReflectModel` for each structure associated with the
        all structures of the Sample, and returns these in a list.
        """
        dq_iter = (
            self.dq
            if isinstance(self.dq, (list, tuple))
            else repeat(self.dq)
        )

        return [
            (
                PolarisedReflectModel(
                    structure,
                    scales=scale,
                    bkgs=bkg,
                    dq=dq,
                )
                if self.is_magnetic()
                else ReflectModel(
                    structure,
                    scale=scale,
                    bkg=bkg,
                    dq=dq,
                )
            )
            for structure, scale, bkg, dq in zip(
                self.structures, self.scale, self.bkg, dq_iter
            )
        ]

    def _spin_structure(self, slab, spin):
        """Convert a refnx MagneticSlab into a spin-dependent Slab."""
        sld_param = slab.sld.parameters[0]

        if spin == 'up':
            sld_value = sld_param + slab.rhoM

        elif spin == 'down':
            sld_value = sld_param - slab.rhoM

        else:
            raise ValueError("spin must be 'up' or 'down'")

        return Slab(
            thick=slab.thick,
            sld=sld_value,
            rough=slab.rough,
            vfsolv=0,
            name=f'{slab.name} {spin}',
            interface=slab.interfaces,
        )

    def simulate_reflectivity(self, angle_times,
                              inst_or_path='OFFSPEC') -> None:
        """
        Plot a simulated reflectivity curve given a set of `angle_times` and
        the neutron instrument.

        Args:
            angle_times (list): points and times for each angle.
            inst_or_path (str): either the name of an instrument already in
                                HOGBEN or the path to a direct beam file,
                                defaults to 'OFFSPEC'
        """
        if not isinstance(angle_times[0], list):
            angle_times = [angle_times for _ in self.get_models()]

        # Plot the model and simulated reflectivity against Q.
        fig = plt.figure()
        ax = fig.add_subplot(111)
        current_xmax = 0

        for i, model in enumerate(self.get_models()):
            data = SimulateReflectivity(model, angle_times[i],
                                        inst_or_path).simulate()
            # Extract each column of the simulated `data`.
            q, r, dr, _ = data[0], data[1], data[2], data[3]

            # Calculate the model reflectivity.
            r_model = SimulateReflectivity(model, angle_times[i],
                                           inst_or_path).reflectivity(q)

            label = f', {self.labels[i]}' if len(self.structures) > 1 else ''

            # Model reflectivity.
            ax.plot(q, r_model, zorder=20, label=f'Model Reflectivity{label}')

            # Simulated reflectivity
            ax.errorbar(q, r, dr, marker='o', ms=3, lw=0,
                        elinewidth=1, capsize=1.5,
                        label=f'Simulated Data{label}',
                        )
            if max(q) > current_xmax:
                current_xmax = max(q)

    def plot_sensitivity_profile(self, q=None, counts=None, show=True,
                                 ax=None):
        """
        Plot the average squared sensitivity divided by reflectivity.

        For each Q value the sensitivity is calculated as the mean of the
        squared parameter gradients, where each parameter is normalized by its
        importance value (defaulting to 1 if not set) and then divided by the
        reflectivity at that Q point.

        Args:
            q (array-like): Q values to evaluate. Defaults to a logarithmic
                span from 0.001 to 0.3 Å^-1.
            counts (array-like): incident count values corresponding to each Q
                point. Currently unused.
            show (bool): whether to display the plot immediately.
            ax (matplotlib.axes.Axes): optional axes object to draw on.

        Returns:
            matplotlib.axes.Axes: the axes containing the sensitivity profile.
        """
        if q is None:
            q = np.geomspace(0.001, 0.3, 1000)

        q = np.asarray(q, dtype=float)

        if q.ndim != 1:
            raise ValueError('q must be a one-dimensional array')

        xi = self.get_param_by_attribute('vary')

        if len(xi) == 0:
            xi = self.get_param_by_attribute('optimize')

        if len(xi) == 0:
            raise ValueError(
                'No varying parameters are available for plotting'
            )

        models = self.get_models()

        if not models:
            raise ValueError('No models are available for plotting')

        if ax is None:
            _, ax = plt.subplots()

        ax_reflectivity = ax.twinx()

        for index, model in enumerate(models):

            q_values = np.asarray(q.copy(), dtype=float)

            fisher = Fisher(
                qs=[q_values],
                xi=xi,
                counts=None,
                models=[model],
            )

            J = fisher._get_gradient_matrix()

            reflectivity = SimulateReflectivity(
                model
            ).reflectivity(q_values)

            importance = np.array(
                [
                    param.importance
                    if hasattr(param, 'importance')
                    else 1.0
                    for param in xi
                ],
                dtype=float,
            )

            sensitivity = np.mean(
                (J ** 2) / importance[np.newaxis, :],
                axis=1,
            )

            sensitivity /= np.clip(
                np.abs(reflectivity),
                1e-30,
                None,
            )

            if len(self.labels) > index:
                base_label = self.labels[index]
            else:
                base_label = f'Model {index + 1}'

            sens_line, = ax.plot(
                q_values,
                sensitivity,
                label=f'{base_label} Sensitivity',
            )

            ax_reflectivity.plot(
                q_values,
                np.abs(reflectivity),
                linestyle='--',
                color=sens_line.get_color(),
                alpha=0.8,
                label=f'{base_label} Reflectivity',
            )

        ax.set_xlabel(r'$\mathregular{Q\ (Å^{-1})}$')
        ax.set_ylabel('Sensitivity (Arb. Units)')
        ax_reflectivity.set_ylabel('Reflectivity')

        ax.set_title('Sensitivity Profile')

        ax.set_yscale('log')
        ax_reflectivity.set_yscale('log')

        ax.set_xlim(q.min(), q.max())

        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax_reflectivity.get_legend_handles_labels()

        ax.legend(
            lines1 + lines2,
            labels1 + labels2,
            loc='best',
        )

        if show:
            plt.show()

        return ax


    @abstractmethod
    def nested_sampling(self):
        """Runs nested sampling on measured or simulated data of the sample."""
        pass


class BaseLipid(BaseSample, VariableContrast, VariableUnderlayer):
    """Abstract class representing the base class for a lipid model."""

    def __init__(self):
        """
        Initialize a BaseLipid object sample, and loads the
        experimentally measured data
        """
        super().__init__()
        self._create_objectives()  # Load experimentally-measured data.

    @abstractmethod
    def _create_objectives(self):
        """Loads the measured data for the lipid sample."""
        pass

    def angle_info(self, angle_times, contrasts=None,
                   inst_or_path='OFFSPEC'):
        """Calculates the Fisher information matrix for the lipid sample
           measured over a number of angles.

        Args:
            angle_times (list): points and times for each angle to simulate.
            contrasts (list): SLDs of contrasts to simulate.

        Returns:
            numpy.ndarray: Fisher information matrix.

        """
        return self.__conditions_info(angle_times, contrasts, None,
                                      inst_or_path)

    def contrast_info(self, angle_times, contrasts, inst_or_path='OFFSPEC'):
        """Calculates the Fisher information matrix for the lipid sample
           with contrasts measured over a number of angles.

        Args:
            angle_times (list): points and times for each angle to simulate.
            contrasts (list): SLDs of contrasts to simulate.

        Returns:
            numpy.ndarray: Fisher information matrix.

        """
        return self.__conditions_info(angle_times, contrasts, None,
                                      inst_or_path)

    def underlayer_info(self, angle_times, contrasts, underlayers,
                        inst_or_path='OFFSPEC'):
        """Calculates the Fisher information matrix for the lipid sample with
           `underlayers`, and `contrasts` measured over a number of angles.

        Args:
            angle_times (list): points and times for each angle to simulate.
            contrasts (list): SLDs of contrasts to simulate.
            underlayers (list): thickness and SLD of each underlayer to add.

        Returns:
            numpy.ndarray: Fisher information matrix.

        """
        return self.__conditions_info(angle_times, contrasts, underlayers,
                                      inst_or_path)

    def __conditions_info(self, angle_times, contrasts, underlayers,
                          inst_or_path='OFFSPEC'):
        """Calculates the Fisher information object for the lipid sample
           with given conditions.

        Args:
            angle_times (list): points and times for each angle to simulate.
            contrasts (list): SLDs of contrasts to simulate.
            underlayers (list): thickness and SLD of each underlayer to add.
            inst_or_path (str): instrument or path to direct beam file.

        Returns:
            Fisher: Fisher information matrix object

        """
        # Iterate over each contrast to simulate.
        qs, counts, models = [], [], []

        for contrast in contrasts:
            # Simulate data for the contrast.
            sample = self._using_conditions(contrast, underlayers)
            contrast_point = (contrast + 0.56) / (6.35 + 0.56)
            background_level = (2e-6 * contrast_point
                                + 4e-6 * (1 - contrast_point))
            model = (
                PolarisedReflectModel(
                    sample,
                    scales=1.0,
                    bkgs=background_level,
                    dq=2,
                )
                if self.is_magnetic()
                else ReflectModel(
                    sample,
                    scale=1.0,
                    bkg=background_level,
                    dq=2,
                )
            )
            data = SimulateReflectivity(model, angle_times,
                                        inst_or_path).simulate()
            qs.append(data[0])
            counts.append(data[3])
            models.append(model)

        # Exclude certain parameters if underlayers are being used.
        if underlayers is None:
            return Fisher(qs, self.params, counts, models)
        else:
            return Fisher(qs, self.underlayer_params, counts, models)

    @abstractmethod
    def _using_conditions(self):
        """Creates a structure describing the given measurement conditions."""
        pass

    def sld_profile(self,
                    save_path: str,
                    filename: str = 'sld_profile',
                    ylim: Optional[tuple] = None,
                    legend: bool = True) -> None:
        """Plots the SLD profile of the lipid sample.

        Args:
            save_path (str): path to directory to save SLD profile to.
            filename (str): file name to use when saving the SLD profile.
            ylim (tuple): limits to place on the SLD profile y-axis.
            legend (bool): whether to include a legend in the SLD profile.

        """
        fig = plt.figure()
        ax = fig.add_subplot(111)

        # Plot the SLD profile for each measured contrast.
        for structure in self.structures:
            # Set SLD plotting limits based on slab thickness and roughness
            # zmax: total slab thickness + 4 × last slab roughness, +5% margin
            # zmin: extend below zero by 4 × first slab roughness, -5%
            # This ensures full SLD profile is visible with extra margin
            # Use 500 steps for good resolution across the thickness range

            zmax = 1.05 * (structure.slabs()[:, 0].sum()
                           + 4 * int(structure.slabs()[-1, 3]))
            zmin = (0 - 4 * int(structure.slabs()[1, 3])) - 0.05 * zmax
            zsteps = 500
            ax.plot(*structure.sld_profile(np.linspace(zmin, zmax, zsteps)))

        x_label = '$\mathregular{Distance\ (\AA)}$'
        y_label = '$\mathregular{SLD\ (10^{-6} \AA^{-2})}$'
        ax.set_xlabel(x_label, fontsize=11, weight='bold')
        ax.set_ylabel(y_label, fontsize=11, weight='bold')

        # Limit the y-axis if specified.
        if ylim:
            ax.set_ylim(*ylim)

        # Add a legend if specified.
        if legend and len(self.structures) > 1:
            ax.legend(self.labels, loc='upper left')

        # Save the plot.
        save_path = os.path.join(save_path, self.name)
        save_plot(fig, save_path, filename)

    def reflectivity_profile(self,
                             save_path: str,
                             filename: str = 'reflectivity_profile') -> None:
        """Plots the reflectivity profile of the lipid sample.

        Args:
            save_path (str): path to directory to save profile to.
            filename (str): file name to use when saving the profile.

        """
        fig = plt.figure()
        ax = fig.add_subplot(111)

        # Iterate over each measured contrast.
        colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for i, objective in enumerate(self.objectives):
            # Get the measured data and calculate the model reflectivity.
            q, r, dr = objective.data.x, objective.data.y, objective.data.y_err
            r_model = objective.model(q)

            # Offset the data, for clarity.
            offset = 10 ** (-2 * i)
            r *= offset
            dr *= offset
            r_model *= offset

            # Add the offset in the label.
            label = self.labels[i]
            if offset != 1:
                label += ' $\\mathregular{(x10^{-' + str(2 * i) + '})}$'

            # Plot the measured data and the model reflectivity.
            ax.errorbar(q, r, dr,
                        marker='o', ms=3, lw=0, elinewidth=1, capsize=1.5,
                        color=colours[i], label=label)
            ax.plot(q, r_model, color=colours[i], zorder=20)

        x_label = '$\\mathregular{Q\\ (Å^{-1})}$'
        y_label = 'Reflectivity (arb.)'

        ax.set_xlabel(x_label, fontsize=11, weight='bold')
        ax.set_ylabel(y_label, fontsize=11, weight='bold')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(1e-10, 3)
        ax.set_title('Reflectivity profile')
        if len(self.structures) > 1:
            ax.legend()

        # Save the plot.
        save_path = os.path.join(save_path, self.name)
        save_plot(fig, save_path, filename)

    def nested_sampling(self,
                        contrasts: list,
                        angle_times: list,
                        save_path: str,
                        filename: str,
                        underlayers=None,
                        dynamic=False,
                        inst_or_path='OFFSPEC') -> None:
        """Runs nested sampling on simulated data of the lipid sample.

        Args:
            contrasts (list): SLDs of contrasts to simulate.
            angle_times (list): points and times for each angle to simulate.
            save_path (str): path to directory to save corner plot to.
            filename (str): file name to use when saving corner plot.
            underlayers (list): thickness and SLD of each underlayer to add.
            dynamic (bool): whether to use static or dynamic nested sampling.
            inst_or_path (str): instrument or path to direct beam file.
        """
        # Create objectives for each contrast to sample with.
        objectives = []
        for contrast in contrasts:
            # Simulate an experiment using the given contrast.
            sample = self._using_conditions(contrast, underlayers)
            contrast_point = (contrast + 0.56) / (6.35 + 0.56)
            background_level = 2e-6 * contrast_point + 4e-6 * (
                1 - contrast_point)

            model = (
                PolarisedReflectModel(
                    sample,
                    scales=1.0,
                    bkgs=background_level,
                    dq=2,
                )
                if self.is_magnetic()
                else ReflectModel(
                    sample,
                    scale=1.0,
                    bkg=background_level,
                    dq=2,
                )
            )
            data = SimulateReflectivity(model, angle_times,
                                        inst_or_path).simulate()
            
            # filter zeros as nested sampling doesn't deal with these well
            data = data[:, (data[1] != 0)]

            dataset = refnx.dataset.ReflectDataset(
                [data[0], data[1], data[2]]
            )
            objectives.append(refnx.analysis.Objective(model, dataset))

        # Combine objectives into a single global objective.
        global_objective = refnx.analysis.GlobalObjective(objectives)

        # Exclude certain parameters if underlayers are being used.
        if underlayers is None:
            global_objective.varying_parameters = lambda: self.params
        else:
            global_objective.varying_parameters = (
                lambda: self.underlayer_params
            )

        # Sample the objective using nested sampling.
        sampler = Sampler(global_objective)
        fig = sampler.sample(dynamic=dynamic)

        # Save the sampling corner plot.
        save_path = os.path.join(save_path, self.name)
        save_plot(fig, save_path, 'nested_sampling_' + filename)
