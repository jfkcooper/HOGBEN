"""Base classes for different sample types """

import os
from abc import ABC, abstractmethod
from typing import Optional

import matplotlib.pyplot as plt

import refnx.dataset
import refnx.reflect
import refnx.analysis
from refnx.reflect import ReflectModel
from refnx.reflect.structure import Slab
from refnx._lib import flatten

from hogben.simulate import SimulateReflectivity
from hogben.utils import Fisher, Sampler, save_plot

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
                    if isinstance(layer, MagneticSLD):
                        up_structure[i] = layer.spin_up
                        down_structure[i] = layer.spin_down
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
            for layer in structure:
                if isinstance(layer, MagneticSLD):
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
        return [refnx.reflect.ReflectModel(structure,
                                           scale=scale,
                                           bkg=bkg,
                                           dq=dq)
                for structure, scale, bkg, dq
                in zip(self.structures, self.scale, self.bkg, self.dq)]

    def simulate_reflectivity(self, angle_times,
                              inst_or_path='OFFSPEC') -> None:
        """
        Plot a simulated reflectivity curve given a set of `angle_times` and
        the neutron instrument.

        Args:
            angle_times (list): points and times for each angle.
            inst_or_path (str): either the name of an instrument already in ,
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
                        elinewidth=1, capsize=1.5, label='Simulated Data'
                                                         f'{label}')
            if max(q) > current_xmax:
                current_xmax = max(q)

        ax.set_xlabel('$\mathregular{Q\ (Å^{-1})}$',
                      weight='bold')
        ax.set_ylabel('Reflectivity (arb.)', weight='bold')
        ax.set_yscale('log')
        ax.set_title('Reflectivity Profile')
        ax.set_xlim(0, 1.05 * current_xmax)
        ax.legend()

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

    def angle_info(self, angle_times, contrasts):
        """Calculates the Fisher information matrix for the lipid sample
           measured over a number of angles.

        Args:
            angle_times (list): points and times for each angle to simulate.
            contrasts (list): SLDs of contrasts to simulate.

        Returns:
            numpy.ndarray: Fisher information matrix.

        """
        return self.__conditions_info(angle_times, contrasts, None)

    def contrast_info(self, angle_times, contrasts):
        """Calculates the Fisher information matrix for the lipid sample
           with contrasts measured over a number of angles.

        Args:
            angle_times (list): points and times for each angle to simulate.
            contrasts (list): SLDs of contrasts to simulate.

        Returns:
            numpy.ndarray: Fisher information matrix.

        """
        return self.__conditions_info(angle_times, contrasts, None)

    def underlayer_info(self, angle_times, contrasts, underlayers):
        """Calculates the Fisher information matrix for the lipid sample with
           `underlayers`, and `contrasts` measured over a number of angles.

        Args:
            angle_times (list): points and times for each angle to simulate.
            contrasts (list): SLDs of contrasts to simulate.
            underlayers (list): thickness and SLD of each underlayer to add.

        Returns:
            numpy.ndarray: Fisher information matrix.

        """
        return self.__conditions_info(angle_times, contrasts, underlayers)

    def __conditions_info(self, angle_times, contrasts, underlayers):
        """Calculates the Fisher information object for the lipid sample
           with given conditions.

        Args:
            angle_times (list): points and times for each angle to simulate.
            contrasts (list): SLDs of contrasts to simulate.
            underlayers (list): thickness and SLD of each underlayer to add.

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
            model = ReflectModel(sample)
            model.bkg = background_level
            model.dq = 2
            data = SimulateReflectivity(model, angle_times).simulate()
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
            ax.plot(*structure.sld_profile(self.distances))

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
                        dynamic=False) -> None:
        """Runs nested sampling on simulated data of the lipid sample.

        Args:
            contrasts (list): SLDs of contrasts to simulate.
            angle_times (list): points and times for each angle to simulate.
            save_path (str): path to directory to save corner plot to.
            filename (str): file name to use when saving corner plot.
            underlayers (list): thickness and SLD of each underlayer to add.
            dynamic (bool): whether to use static or dynamic nested sampling.

        """
        # Create objectives for each contrast to sample with.
        objectives = []
        for contrast in contrasts:
            # Simulate an experiment using the given contrast.
            sample = self._using_conditions(contrast, underlayers)
            contrast_point = (contrast + 0.56) / (6.35 + 0.56)
            background_level = 2e-6 * contrast_point + 4e-6 * (
                1 - contrast_point)

            model = ReflectModel(sample)
            model.bkg = background_level
            model.dq = 2
            data = SimulateReflectivity(model, angle_times).simulate()
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


class MagneticSLD(Slab):
    """
    A class to represent a layer with a magnetic SLD component.

    This class extends the `Slab` class from `refnx.reflect.structure` to
    include properties for magnetic Scattering Length Density (SLD) as well.

    Attributes:
        SLD_n (float): Nuclear scattering length density.
        SLD_m (float): Magnetic scattering length density.
        thickness (float): Thickness of the layer.
        roughness (float): Roughness of the layer.
        name (str): Name of the layer.
    """

    def __init__(self,
                 SLDn: float = 0,
                 SLDm: float = 0,
                 thick: float = 0,
                 rough: float = 0,
                 vfsolv: float = 0,
                 interface: refnx.reflect.interface = None,
                 name: str = 'Magnetic Layer'):
        """
        Initialize a MagneticSLD object.

        Parameters:
            SLDn (float): Nuclear scattering length density. Default is 0.
            SLDm (float): Magnetic scattering length density. Default is 0.
            thick (float): Thickness of the layer. Default is 0.
            rough (float): Roughness of the layer. Default is 0.
            vfsolv (float): Volume fraction of the solvent, between 0 and 1.
                            Default is 0.
            name (str): Name of the layer. Default is "Magnetic Layer".
            interface (`refnx.reflect.Interface`):
                The type of interfacial roughness associated with the Slab.
                If `None`, then the default interfacial roughness is an Error
                function (also known as Gaussian roughness).
        """
        self.SLDn = SLDn
        self.SLDm = SLDm
        self.thick = thick
        self.rough = rough
        self.vfsolv = vfsolv
        self.interface = interface
        self.name = name
        super().__init__(thick=self.thick, sld=self.SLDn,
                         rough=self.rough, name=self.name, vfsolv=self.vfsolv,
                         interface=self.interface)

    @property
    def spin_up(self):
        """
        Calculate the spin-up scattering length density.

        Returns:
            SLD: An SLD object representing the spin-up component with the
            appropriate thickness and roughness.
        """
        SLD_value = self.SLDn + self.SLDm
        return Slab(thick=self.thick, sld=SLD_value, rough=self.rough,
                    vfsolv=self.vfsolv, name='Spin up',
                    interface=self.interface)

    @property
    def spin_down(self):
        """
        Calculate the spin-down scattering length density.

        Returns:
            SLD: An SLD object representing the spin-down component with the
            appropriate thickness and roughness.
        """
        SLD_value = self.SLDn - self.SLDm
        return Slab(thick=self.thick, sld=SLD_value, rough=self.rough,
                    vfsolv=self.vfsolv, name='Spin down',
                    interface=self.interface)

    def __call__(self, thick=None, rough=None, vfsolv=None):
        """
        Update the thickness and roughness of the layer.

        Parameters:
            thick (float): New thickness of the layer. If None, the
            current thickness is retained.
            rough (float): New roughness of the layer. If None, the
            current roughness is retained.
            vfsolv (float): New volume fraction of the solvent, between 0 and
                            1. Default is 0.

        Returns:
            MagneticSLD: The updated MagneticSLD object.
        """
        self.thick = thick if thick is not None else self.thick
        self.rough = rough if rough is not None else self.rough
        self.vfsolv = vfsolv if vfsolv is not None else self.vfsolv
        super().__init__(thick=self.thick, sld=self.SLDn,
                         rough=self.rough, name=self.name, vfsolv=self.vfsolv,
                         interface=self.interface)
        return self
