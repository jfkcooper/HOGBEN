"""
Contains class and methods related to the Sample class.
"""

import os
from typing import Any, Optional

import matplotlib.pyplot as plt
import numpy as np

import refnx.dataset
import refnx.reflect
import refnx.analysis

from hogben.simulate import SimulateReflectivity
from hogben.utils import Fisher, Sampler, save_plot
from hogben.models.base import BaseSample
from refnx.analysis import Objective, GlobalObjective
from refnx.reflect import ReflectModel

plt.rcParams['figure.figsize'] = (9, 7)
plt.rcParams['figure.dpi'] = 600


class Sample(BaseSample):
    """Wrapper class for a standard refnx reflectometry sample.

    Attributes:
        structures (list[refnx.reflect.Structure]): List of structures in the
                                                    sample
        name (str): name of the sample.
        polarised (bool): Whether the sample is polarised or not.
        labels (list): List of the labels corresponding to each structure
        scale (list): List of the scale corresponding to each structure
        bkg (list): List of the backgrounds corresponding to each structure
        dq (list): List of the resolutions corresponding to each structure
        params (list): The varying parameters in the sample.
    """

    def __init__(self, structures, **settings):
        """
        Initializes a sample given a structure, and sets the sample name and
        parameters

        Args:
            structures: Sample structure defined in the refnx model
            **settings The settings which can be applied to the experiment:
                labels (list): The labels for each structure
                scale (list|float): The scale factor used for the structures
                bkg (list|float): The backgrounds for the structures
                dq (list|float): The resolutions for the structures
                polarised (bool): Whether the sample is polarised
        """
        super().__init__()
        if isinstance(structures, refnx.reflect.Structure):
            structures = [structures]
        self.structures = structures
        self.polarised = settings.get('polarised', self.is_magnetic())
        self.name = ', '.join(
            {structure.name for structure in self.structures})
        self._labels = None
        self.labels = settings.get('labels', self._labels)
        self.scale = settings.get('scale', 1)
        self.bkg = settings.get('bkg', 5e-6)
        self.dq = settings.get('dq', 2)

    def _validate_and_set(self, attribute, value_list):
        """
        Validates the length of the value list against the number of structures
        and sets the attribute accordingly.

        Args:
            attribute (str): The name of the attribute to set.
            value_list (float or list of floats): The value(s) to set for the
                                                  attribute.

        Raises:
            ValueError: If the length of the value list does not match the
                        number of structures.
        """
        if isinstance(value_list, list):
            if self.is_magnetic() and self.polarised:
                # Duplicate each item for magnetic samples (up+down)
                new_list = []
                for item in value_list:
                    new_list.append(item)
                    new_list.append(item)
                value_list = new_list
            if len(value_list) != len(self.structures):
                raise ValueError(
                    f'The length of `{attribute}` must be equal to the number '
                    f'of structures in the sample when using a list!'
                )
            else:
                setattr(self, f'_{attribute}', value_list)
        else:
            setattr(self, f'_{attribute}', [value_list] * len(self.structures))

    @property
    def bkg(self):
        """
        Gets the background levels for the structures.

        Returns:
            list of floats: The background levels.
        """
        return self._bkg

    @bkg.setter
    def bkg(self, value):
        """
        Sets the background levels for the structures.

        Args:
            value (float or list of floats): The background level(s) to set.
        """
        self._validate_and_set('bkg', value)

    @property
    def dq(self):
        """
        Gets the dq values for the structures.

        Returns:
            list of floats: The dq values.
        """
        return self._dq

    @dq.setter
    def dq(self, value):
        """
        Sets the dq values for the structures.

        Args:
            value (float or list of floats): The dq value(s) to set.
        """
        self._validate_and_set('dq', value)

    @property
    def scale(self):
        """
        Gets the scale factors for the structures.

        Returns:
            list of floats: The scale factors.
        """
        return self._scale

    @scale.setter
    def scale(self, value):
        """
        Sets the scale factors for the structures.

        Args:
            value (float or list of floats): The scale factor(s) to set.
        """
        self._validate_and_set('scale', value)

    @property
    def params(self) -> list:
        """List of all varying parameters of the sample"""
        return self.get_param_by_attribute('vary')

    def _vary_structure(self, bound_size: float = 0.2) -> list:
        """Varies the SLD and thickness of each layer in the sample structures.

        Args:
            structure (refnx.reflect.Structure): structure to vary.
            bound_size (float): size of bounds to place on varying parameters.

        Returns:
            list: varying parameters of sample.

        """
        for structure in self.structures:
            params = []
            # Vary the SLD and thickness of each component (layer).
            for component in structure[1:-1]:
                sld = component.sld.real
                if sld.bounds.lb == -float('inf'):
                    sld.bounds.lb = sld.value * (1 - bound_size)
                if sld.bounds.ub == float('inf'):
                    sld.bounds.ub = sld.value * (1 + bound_size)

                sld_bounds = (
                    sld.bounds.lb,
                    sld.bounds.ub,
                )
                sld.setp(vary=True, bounds=sld_bounds)
                params.append(sld)

                thick = component.thick
                thick_bounds = (
                    thick.value * (1 - bound_size),
                    thick.value * (1 + bound_size),
                )
                thick.setp(vary=True, bounds=thick_bounds)
                params.append(thick)

            return params

    @property
    def labels(self) -> list:
        """
        Returns a list of all refnx `ReflectModel` models that are
        associated with each structure of the sample.
        """
        # Use provided labels if present, else generate labels automatically
        if self._labels:
            labels = self._labels
        else:
            labels = []
            # Check if the different structures have the same solvent
            solvent_sld = self.structures[0][-1].sld.real.value
            same_solvent = all(structure[-1].sld.real.value == solvent_sld
                               for structure in self._structures)

            for index, structure in enumerate(self._structures):
                if len(self._structures) == 1:
                    label = ''
                elif same_solvent:
                    label = f'Structure {index}'
                else:
                    sld_value = structure[-1].sld.real.value
                    label = f'Solvent SLD:{"{:.3g}".format(sld_value)}'
                labels.append(label)

        # Add spin-labels for polarised and magnetic samples
        if self.is_magnetic() and self.polarised:
            # Duplicate each label per spin state
            labels = [item for item in labels for _ in range(2)]
            # Add spin-state for every structure
            labels = [
                f'Spin-up, {label}' if index % 2 == 0
                else f'Spin-down, {label}'
                for index, label in enumerate(labels)
            ]
        return labels

    @labels.setter
    def labels(self, labels: list) -> None:
        """
        Returns a list of all refnx `ReflectModel` models that are
        associated with each structure of the sample.
        """

        # Don't try to set the labels if no labels were specified
        if labels is None:
            return

        if (isinstance(labels, list)
                and all(isinstance(label, str) for label in labels)):
            if len(labels) == len(self._structures):
                self._labels = labels
            else:
                raise ValueError(
                    'The amount of labels must be equal to the number '
                    'of structures in the sample!'
                )
        else:
            raise TypeError(
                'The labels need to be given in the form of a list of'
                ' strings!'
            )


    def angle_info(self,
                   angle_times: list[tuple],
                   contrasts: Any | None = None) -> Fisher:
        """Calculates the Fisher information matrix for a sample measured
           over a number of angles.

        Args:
            angle_times (list): points and times for each angle to simulate.

        Returns:
            Fisher: Fisher information object

        """
        # Return the Fisher information matrix calculated from simulated data.
        models = self.get_models()
        qs, counts = [], []
        for model in models:
            data = SimulateReflectivity(model, angle_times).simulate()
            qs.append(data[0])
            counts.append(data[3])
        return Fisher(qs, self.params, counts, models)

    def sld_profile(self, save_path=None):
        """Plots the SLD profile of the sample.

        Args:
            save_path (str): path to directory to save SLD profile to.
        """
        # Create a figure and axes based on the 'single' parameter
        fig, ax = plt.subplots()
        for i, (z, slds) in enumerate(self._get_sld_profile()):
            # Create a new subplot for each profile if not 'single'
            label = f'{self.labels[i]}'

            # Create a new subplot for each profile if not 'single'
            ax.set_xlim(min(z), max(z))
            ax.plot(z, slds, label=label)

            ax.set_xlabel(r'$\mathregular{Distance\ (\AA)}$')
            ax.set_ylabel(r'$\mathregular{SLD\ (10^{-6} \AA^{-2})}$')
            ax.set_title('SLD Profile')
            if len(self.structures) > 1:
                ax.legend()

        # Save the plot.
        if save_path:
            save_path = os.path.join(save_path, self.name)
            save_plot(fig, save_path, 'sld_profile')

    def _get_sld_profile(self):
        """
        Obtains the SLD profile of the sample, in terms of z (depth) vs SLD

        Returns:
            numpy.ndarray: depth
            numpy.ndarray: SLD values
        """
        return [structure.sld_profile() for structure in self.structures]

    def reflectivity_profile(self,
                             save_path: str = None,
                             q_min: float = 0.005,
                             q_max: float = 0.4,
                             points: int = 500,
                             scale: float = 1,
                             bkg: float = 1e-7,
                             dq: float = 2,
                             ) -> None:
        """Plots the reflectivity profile of the sample.

        Args:
            save_path (str): path to directory to save reflectivity profile to.
            q_min (float): minimum Q value to plot.
            q_max (float): maximum Q value to plot.
            points (int): number of points to plot.
            scale (float): experimental scale factor.
            bkg (float): level of instrument background noise.
            dq (float): instrument resolution.

        """
        fig, ax = plt.subplots()
        profiles = self._get_reflectivity_profile(q_min, q_max, points, scale,
                                                  bkg, dq)
        for i, (q, r) in enumerate(profiles):

            # Plot Q versus model reflectivity.
            ax.plot(q, r, label=self.labels[i])

            ax.set_xlabel(r'$\mathregular{Q\ (Å^{-1})}$')
            ax.set_ylabel('Reflectivity (arb.)')
            ax.set_title('Reflectivity profile')
            ax.set_yscale('log')

            if len(self.structures) > 1:
                ax.legend()

        # Save the plot.
        if save_path:
            save_path = os.path.join(save_path, self.name)
            save_plot(fig, save_path, 'reflectivity_profile')

    def _get_reflectivity_profile(self,
                                  q_min: float,
                                  q_max: float,
                                  points: int,
                                  scale: float,
                                  bkg: float,
                                  dq: float) -> list:
        """
        Obtains the reflectivity profile of the sample, in terms of q
        vs r

        Returns:
            numpy.ndarray: q values at each reflectivity point
            numpy.ndarray: model reflectivity values
        """
        profiles = []
        # Geometrically-space Q points over the specified range.
        q = np.geomspace(q_min, q_max, points)

        for structure in self.structures:
            model = ReflectModel(structure, scale=scale,
                                 bkg=bkg, dq=dq)
            r = SimulateReflectivity(model).reflectivity(q)
            profiles.append((q, r))
        return profiles

    def nested_sampling(self,
                        angle_times: list,
                        save_path: str,
                        filename: str,
                        dynamic: bool = False) -> None:
        """Runs nested sampling on simulated data of the sample.

        Args:
            angle_times (list: points and times for each angle to simulate.
            save_path (str): path to directory to save corner plot to.
            filename (str): file name to use when saving corner plot.
            dynamic (bool): whether to use static or dynamic nested sampling.

        """
        # Simulate data for the sample.
        objectives = []
        for structure in self.structures:
            model = refnx.reflect.ReflectModel(structure)
            data = SimulateReflectivity(model, angle_times).simulate()
            # filter zeros as nested sampling doesn't deal with these well
            data = data[:, (data[1] != 0)]
            objective = Objective(model, data)
            objectives.append(objective)
        global_objective = GlobalObjective(objectives)
        global_objective.varying_parameters = lambda: self.params

        # Sample the objective using nested sampling.
        sampler = Sampler(global_objective)
        fig = sampler.sample(dynamic=dynamic)

        # Save the sampling corner plot.
        save_path = os.path.join(save_path, self.name)
        save_plot(fig, save_path, filename + '_nested_sampling')


def simple_sample():
    """Defines a 2-layer simple sample.

    Returns:
        samples.Sample: structure in format for design optimisation.

    """
    air = refnx.reflect.SLD(0, name='Air')
    layer1 = refnx.reflect.SLD(4, name='Layer 1')(thick=100, rough=2)
    layer2 = refnx.reflect.SLD(8, name='Layer 2')(thick=150, rough=2)
    substrate = refnx.reflect.SLD(2.047, name='Substrate')(thick=0, rough=2)

    structure = air | layer1 | layer2 | substrate
    structure.name = 'simple_sample'
    return Sample(structure)


def many_param_sample():
    """Defines a 5-layer sample with many parameters.

    Returns:
        samples.Sample: structure in format for design optimisation.

    """
    air = refnx.reflect.SLD(0, name='Air')
    layer1 = refnx.reflect.SLD(2.0, name='Layer 1')(thick=50, rough=6)
    layer2 = refnx.reflect.SLD(1.7, name='Layer 2')(thick=15, rough=2)
    layer3 = refnx.reflect.SLD(0.8, name='Layer 3')(thick=60, rough=2)
    layer4 = refnx.reflect.SLD(3.2, name='Layer 4')(thick=40, rough=2)
    layer5 = refnx.reflect.SLD(4.0, name='Layer 5')(thick=18, rough=2)
    substrate = refnx.reflect.SLD(2.047, name='Substrate')(thick=0, rough=2)

    structure = air | layer1 | layer2 | layer3 | layer4 | layer5 | substrate
    structure.name = 'many_param_sample'
    return Sample(structure)


def thin_layer_sample_1():
    """Defines a 2-layer sample with thin layers.

    Returns:
        samples.Sample: structure in format for design optimisation.

    """
    air = refnx.reflect.SLD(0, name='Air')
    layer1 = refnx.reflect.SLD(4, name='Layer 1')(thick=200, rough=2)
    layer2 = refnx.reflect.SLD(6, name='Layer 2')(thick=6, rough=2)
    substrate = refnx.reflect.SLD(2.047, name='Substrate')(thick=0, rough=2)

    structure = air | layer1 | layer2 | substrate
    structure.name = 'thin_layer_sample_1'
    return Sample(structure)


def thin_layer_sample_2():
    """Defines a 3-layer sample with thin layers.

    Returns:
        samples.Sample: structure in format for design optimisation.

    """
    air = refnx.reflect.SLD(0, name='Air')
    layer1 = refnx.reflect.SLD(4, name='Layer 1')(thick=200, rough=2)
    layer2 = refnx.reflect.SLD(5, name='Layer 2')(thick=30, rough=6)
    layer3 = refnx.reflect.SLD(6, name='Layer 3')(thick=6, rough=2)
    substrate = refnx.reflect.SLD(2.047, name='Substrate')(thick=0, rough=2)

    structure = air | layer1 | layer2 | layer3 | substrate
    structure.name = 'thin_layer_sample_2'
    return Sample(structure)


def similar_sld_sample_1():
    """Defines a 2-layer sample with layers of similar SLD.

    Returns:
        samples.Sample: structure in format for design optimisation.

    """
    air = refnx.reflect.SLD(0, name='Air')
    layer1 = refnx.reflect.SLD(0.9, name='Layer 1')(thick=80, rough=2)
    layer2 = refnx.reflect.SLD(1.0, name='Layer 2')(thick=50, rough=6)
    substrate = refnx.reflect.SLD(2.047, name='Substrate')(thick=0, rough=2)

    structure = air | layer1 | layer2 | substrate
    structure.name = 'similar_sld_sample_1'
    return Sample(structure)


def similar_sld_sample_2():
    """Defines a 3-layer sample with layers of similar SLD.

    Returns:
        samples.Sample: structure in format for design optimisation.

    """
    air = refnx.reflect.SLD(0, name='Air')
    layer1 = refnx.reflect.SLD(3.0, name='Layer 1')(thick=50, rough=2)
    layer2 = refnx.reflect.SLD(5.5, name='Layer 2')(thick=30, rough=6)
    layer3 = refnx.reflect.SLD(6.0, name='Layer 3')(thick=35, rough=2)
    substrate = refnx.reflect.SLD(2.047, name='Substrate')(thick=0, rough=2)

    structure = air | layer1 | layer2 | layer3 | substrate
    structure.name = 'similar_sld_sample_2'
    return Sample(structure)


def run_main(save_path: Optional[str] = '../results') -> None:
    """
    Runs the main function of the module, retrieves an SLD and
    reflectivity profile for each defined structure, and saves it in the
    results directory by default.

    Args:
        save_path: The directory where the SLD and reflectivity profiles
        are saved
    """
    # Plot the SLD and reflectivity profiles of all structures in this file.
    for structure in [simple_sample, many_param_sample,
                      thin_layer_sample_1, thin_layer_sample_2,
                      similar_sld_sample_1, similar_sld_sample_2]:
        sample = structure()
        sample.sld_profile(save_path)
        sample.reflectivity_profile(save_path)

        # Close the plots.
        plt.close('all')


if __name__ == '__main__':
    run_main()
