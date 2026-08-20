import os

import numpy as np

from dynesty import NestedSampler, DynamicNestedSampler
from dynesty import plotting as dyplot
from dynesty import utils as dyfunc

import refnx.reflect
import refnx.analysis

from hogben.simulate import SimulateReflectivity


class Sampler:
    """Contains code for running nested sampling on refnx models.

    Attributes:
        objective (refnx.analysis.Objective): objective to sample.
        params (list): varying model parameters.
        ndim (int): number of varying model parameters.
        sampler_static (dynesty.NestedSampler): static nested sampler.
        sampler_dynamic (dynesty.DynamicNestedSampler): dynamic nested sampler.

    """

    def __init__(self, objective):
        """
        Initialise the sample given an objective to the sample

        Args:
            objective: objective to the sample
        """
        self.objective = objective

        # Use log-likelihood and prior transform methods of refnx objective
        self.params = objective.varying_parameters()
        logl = objective.logl
        prior_transform = objective.prior_transform

        self.ndim = len(self.params)
        self.sampler_static = NestedSampler(logl, prior_transform, self.ndim)
        self.sampler_dynamic = DynamicNestedSampler(logl, prior_transform,
                                                    self.ndim)

    def sample(self, verbose=True, dynamic=False):
        """Samples an Objective/FitProblem using nested sampling.

        Args:
            verbose (bool): whether to display sampling progress.
            dynamic (bool): whether to use static or dynamic nested sampling.

        Returns:
            matplotlib.pyplot.Figure or float: corner plot.

        """
        # Run either static or dynamic nested sampling.
        if dynamic:
            # Weighting is entirely on the posterior (0 weight on evidence).
            self.sampler_dynamic.run_nested(print_progress=verbose,
                                            wt_kwargs={'pfrac': 1.0})
            results = self.sampler_dynamic.results

        else:
            self.sampler_static.run_nested(print_progress=verbose)
            results = self.sampler_static.results

        # Calculate the parameter means.
        weights = np.exp(results.logwt - results.logz[-1])
        mean, _ = dyfunc.mean_and_cov(results.samples, weights)

        # Set the parameter values to the estimated means.
        for i, param in enumerate(self.params):
            param.value = mean[i]

        # Return the corner plot
        return self.__corner(results)

    def __corner(self, results):
        """Calculates a corner plot from given nested sampling `results`.

        Args:
            results (dynesty.results.Results): full output of a sampling run.

        Returns:
            matplotlib.pyplot.Figure: nested sampling corner plot.

        """
        # Get the corner plot from dynesty package.
        fig, _ = dyplot.cornerplot(results, color='blue', quantiles=None,
                                   show_titles=True, max_n_ticks=3,
                                   truths=np.zeros(self.ndim),
                                   truth_color='black')

        # Label the axes with parameter labels.
        axes = np.reshape(np.array(fig.get_axes()), (self.ndim, self.ndim))
        for i in range(1, self.ndim):
            for j in range(self.ndim):
                if i == self.ndim - 1:
                    axes[i, j].set_xlabel(self.params[j].name)
                if j == 0:
                    axes[i, j].set_ylabel(self.params[i].name)

        axes[self.ndim - 1, self.ndim - 1].set_xlabel(self.params[-1].name)
        return fig


class Fisher():
    """Calculates the Fisher information matrix for multiple `models`
    containing parameters `xi`. The model describes the experiment,
    including the sample, and is defined using `refnx`. The
    lower and upper bounds of each parameter in the model are transformed
    into a standardized range from 0 to 1, which is used to calculate the
    Fisher information matrix. Each parameter in the Fisher information
    matrix is scaled using an importance parameter. By default,
    the importance parameter is set to 1 for all parameters, and can be set
    by changing the `importance` attribute of the parameter when setting up
    the model. For example the relative importance of the thickness in
    "layer1" can be set to 2 using `layer1.thickness.importance = 2` in
     `refnx`.

    Attributes:
        qs: The Q points for each model.
        xi: The varying model parameters.
        counts: incident neutron counts corresponding to each Q value.
        models: models to calculate gradients with.
        step: step size to take when calculating gradient.
        fisher_information: The Fisher information matrix
        min_eigenval: The minimum eigenvalue of the Fisher information matrix
    """

    def __init__(self,
                 qs: list[np.ndarray],
                 xi: list[refnx.analysis.Parameter],
                 counts: list[int],
                 models: list[refnx.reflect.ReflectModel],
                 step: float = 0.005):
        """Initialize the Fisher matrix class.

        Args:
            qs: The Q points for each model.
            xi: The varying model parameters.
            counts: incident neutron counts corresponding to each Q value.
            models: models to calculate gradients with.
            step: step size to take when calculating gradient.
        """
        self.qs = qs
        self.xi = xi
        self.counts = counts
        self.models = models
        self.step = step

    @classmethod
    def from_sample(cls,
                    sample,
                    angle_times,
                    inst_or_path='OFFSPEC'):
        """
        This class method constructs a Fisher object using a HOGBEN sample.

        Args:
        sample (BaseSample): The sample object to calculate the FI for.
        angle_times (list [tuple]): points and times for each angle to
                                    simulate.
        inst_or_path (str): The instrument or path to use for the simulation.
                            Defaults to 'OFFSPEC'.

        Returns:
        Fisher
            A new instance of the Fisher class
        """
        qs, counts = [], []
        models = sample.get_models()

        if not isinstance(angle_times[0], list):
            angle_times = [angle_times for _ in models]

        for model, angle_time in zip(models, angle_times):
            sim = SimulateReflectivity(model, angle_time,
                                       inst_or_path=inst_or_path)
            data = sim.simulate()
            qs.append(data[0])
            counts.append(data[3])
        xi = sample.get_param_by_attribute('vary')
        return cls(qs, xi, counts, models)

    @property
    def fisher_information(self) -> np.ndarray:
        """Calculate and return the Fisher information matrix.

        Returns:
            numpy.ndarray: The Fisher information matrix.
        """
        return self._calculate_fisher_information()

    @property
    def min_eigenval(self) -> float:
        """Calculate and return the minimum eigenvalue of the Fisher
        information matrix.

        Returns:
            float: The minimum eigenvalue.
        """
        return np.linalg.eigvalsh(self.fisher_information)[0]

    @property
    def n(self) -> int:
        """The total number of datapoints.

        Returns:
            int: total number of datapoints.
        """
        return sum(len(q) for q in self.qs)

    @property
    def m(self) -> int:
        """The total number of parameters.

        Returns:
            int: total number of parameters.
        """
        return len(self.xi)

    def _calculate_fisher_information(self) -> np.ndarray:
        """Calculates the Fisher information matrix using the class attributes.

        Returns:
            numpy.ndarray: The Fisher information matrix.
        """
        if self.n == 0:
            return np.zeros((self.m, self.m))
        J = self._get_gradient_matrix()
        # Calculate the reflectance for each model for the given Q values.
        r = np.concatenate([SimulateReflectivity(model).reflectivity(q)
                            for q, model in list(zip(self.qs, self.models))])
        # Calculate the Fisher information matrix using equations from
        # the paper.
        M = np.diag(np.concatenate(self.counts) / r, k=0)
        g = np.dot(np.dot(J.T, M), J)
        # Perform unit scaling if there's at least one parameter
        if len(self.xi) >= 1:
            g = self._scale_units(g)  # Scale by unit bounds
            g = self._scale_importance(g)  # Scale by importance

        return g

    def _scale_units(self, g: np.ndarray) -> np.ndarray:

        """Scale the values of the fisher information matrix for each parameter
        from interval [lb, ub] to the interval [0, 1]

        Args:
            g: The Fisher information matrix.

        Returns:
            numpy.ndarray: The scaled Fisher information matrix.
        """
        lb, ub = self._get_bounds()
        H = np.diag(1 / (ub - lb))  # Get unit scaling Jacobian.
        return np.dot(np.dot(H.T, g), H)  # Perform unit scaling.

    def _scale_importance(self, g: np.ndarray) -> np.ndarray:
        """Scale the Fisher information matrix using importance scaling.

        Args:
            g: The Fisher information matrix.

        Returns:
            numpy.ndarray: The scaled Fisher information matrix.
        """
        importance_array = [param.importance if hasattr(param, 'importance')
                            else 1 for param in self.xi]
        # Importance gets squarerooted and the matrix multiplied
        # by it twice to keep it symetric.
        importance = np.sqrt(np.diag(importance_array))
        return np.dot(np.dot(importance, g), importance)

    def _get_gradient_matrix(self) -> np.ndarray:
        """Calculate the gradient matrix.

        Returns:
            numpy.ndarray: The gradient matrix.
        """
        J = np.zeros((self.n, self.m))
        for i, parameter in enumerate(self.xi):
            old = parameter.value

            # Calculate reflectance for each model for first part of gradient.
            x1 = parameter.value = old * (1 - self.step)
            y1 = np.concatenate([SimulateReflectivity(model).reflectivity(q)
                                 for q, model in list(zip(self.qs,
                                                          self.models))]
                                )

            # Calculate reflectance for each model for second part of gradient.
            x2 = parameter.value = old * (1 + self.step)
            y2 = np.concatenate([SimulateReflectivity(model).reflectivity(q)
                                 for q, model in list(zip(self.qs,
                                                          self.models))]
                                )
            parameter.value = old  # Reset the parameter.
            J[:, i] = (y2 - y1) / (x2 - x1)  # Calculate the gradient.
        return J

    def _get_bounds(self) -> tuple[np.ndarray, np.ndarray]:
        """Get the bounds from the refnx parameters.

        Returns:
            tuple: The lower and upper bounds of the parameters.
        """
        if isinstance(self.xi[0], refnx.analysis.Parameter):
            lb = np.array([param.bounds.lb for param in self.xi])
            ub = np.array([param.bounds.ub for param in self.xi])
        else:
            raise RuntimeError('Invalid sample given')
        return lb, ub


def save_plot(fig, save_path, filename):
    """Saves a figure to a given directory.

    Args:
        fig (matplotlib.pyplot.Figure): figure to save.
        save_path (str): path to directory to save figure to.
        filename (str): name of file to save plot as.

    """
    # Create the directory if not present.
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    file_path = os.path.join(save_path, filename + '.png')
    fig.savefig(file_path, dpi=600)
