"""Functions for visualizing and optimizing angle and contrast for magnetic
experiments."""

import matplotlib.pyplot as plt
import os

import numpy as np

from scipy.optimize import differential_evolution

from hogben.models.magnetic import SampleYIG
from hogben.utils import save_plot

plt.rcParams['figure.figsize'] = (9, 7)
plt.rcParams['figure.dpi'] = 600


def _magnetism_results_visualise(save_path):
    """Visualises the choice of YIG and Pt layer thicknesses for the
       magnetic YIG sample.

    Args:
        save_path (str): path to directory to save results to.

    """
    sample = SampleYIG()
    sample.pt_mag.value = 0.01638

    # Number of points and counting times for each angle to simulate.
    angle_times = [(0.5, 100, 20), (1.0, 100, 40), (2.0, 100, 80)]

    # Range of YIG and Pt layer thicknesses to calculate over.
    yig_thick_range = np.linspace(400, 900, 75)
    pt_thick_range = np.concatenate(
        (np.linspace(21.5, 30, 50), np.linspace(30, 100, 50))
    )

    # Iterate over each YIG and Pt layer thickness being considered.
    x, y, infos = [], [], []
    n = len(pt_thick_range) * len(yig_thick_range)  # Number of calculations.
    for i, yig_thick in enumerate(yig_thick_range):
        # Display progress.
        if i % 5 == 0:
            print('>>> {0}/{1}'.format(i * len(pt_thick_range), n))

        for pt_thick in pt_thick_range:
            # Calculate the Fisher information using current thicknesses.
            g = sample.underlayer_info(angle_times, yig_thick, pt_thick)

            infos.append(g[0, 0])
            x.append(yig_thick)
            y.append(pt_thick)

    # Create plot of YIG and Pt layer thicknesses versus Fisher information.
    fig = plt.figure(figsize=[10, 8])
    ax = fig.add_subplot(111, projection='3d')

    # Create the surface plot and add colour bar.
    surface = ax.plot_trisurf(x, y, infos, cmap='plasma')
    fig.colorbar(surface, fraction=0.046, pad=0.04)

    x_label = r'$\mathregular{YIG\ Thickness\ (\AA)}$'
    y_label = r'$\mathregular{Pt\ Thickness\ (\AA)}$'
    z_label = 'Fisher Information'

    ax.set_xlabel(x_label, fontsize=11, weight='bold')
    ax.set_ylabel(y_label, fontsize=11, weight='bold')
    ax.set_zlabel(z_label, fontsize=11, weight='bold')
    ax.ticklabel_format(axis='z', style='sci', scilimits=(0, 0))

    # Save the plot.
    save_path = os.path.join(save_path, sample.name)
    save_plot(fig, save_path, 'underlayer_choice')

    # Save different views of the 3D plot.
    # Iterate from 0 to 360 degrees in increments of 10.
    save_path = os.path.join(save_path, 'underlayer_choice')
    for i in range(0, 360, 10):
        # Set the "camera" view for the 3D plot.
        ax.view_init(elev=40, azim=i)
        # Save the view.
        save_plot(fig, save_path, 'underlayer_choice_{}'.format(i))

    # Return the YIG and Pt thicknesses with greatest Fisher information.
    maximum = np.argmax(infos)
    return x[maximum], y[maximum]


def _magnetism_results_optimise(save_path):
    """Optimises the choice of YIG and Pt layer thicknesses for the
       magnetic YIG sample.

    Args:
        save_path (str): path to directory to save results to.

    """
    sample = SampleYIG()
    sample.pt_mag.value = 0.01638

    # Number of points and counting times for each angle to simulate.
    angle_times = [(0.5, 100, 20), (1.0, 100, 40), (2.0, 100, 80)]

    # Itervals of YIG and Pt layer thicknesses to optimise over.
    yig_thick_bounds = (200, 900)
    pt_thick_bounds = (21.1, 100)
    bounds = [yig_thick_bounds, pt_thick_bounds]

    # Arguments for optimisation function.
    args = [sample, angle_times]

    # Optimise the YIG and Pt layer thicknesses and save the results.
    save_path = os.path.join(save_path, sample.name)
    file_path = os.path.join(save_path, 'optimised_underlayers.txt')
    with open(file_path, 'w') as file:
        # Optimise thicknesses using differential evolution.
        res = differential_evolution(
            _optimisation_func,
            bounds,
            args=args,
            polish=False,
            tol=0.001,
            updating='deferred',
            workers=-1,
            disp=False,
        )

        yig_thick, pt_thick = res.x[0], res.x[1]

        # Write the unoptimised function value to the .txt file.
        val = _optimisation_func([None, None], sample, angle_times)
        file.write('-------------- Unoptimised  --------------\n')
        file.write('Fisher Information: {}\n\n'.format(-val))

        # Write the thicknesses and function value to the .txt file.
        file.write('-------------- Optimised --------------\n')
        file.write('YIG Thickness: {}\n'.format(round(yig_thick, 1)))
        file.write('Pt Thickness: {}\n'.format(round(pt_thick, 1)))
        file.write('Fisher Information: {}\n\n'.format(-res.fun))


def _optimisation_func(x, sample, angle_times):
    """Defines the function for optimising YIG and Pt layer thicknesses for
       the magnetic YIG sample.

    Args:
        x (list): YIG and Pt layer thicknesses.
        sample (magnetic.SampleYIG): magnetic model being optimised.
        angle_times (list): points and times for each angle to simulate.

    Returns:
        float: negative of the Fisher information in the Pt layer magnetic
               SLD using the given conditions `x`.

    """
    # Calculate the Fisher information matrix using the given thicknesses.
    g = sample.underlayer_info(angle_times, x[0], x[1])
    return -g[0, 0]


def _magnetism_results_ratios(save_path):
    """Calculates log ratio of likelihoods between two models, one with an
      induced moment in the YIG sample Pt layer and one with no moment, as
      a function of measurement time, to determine what level of statistics
      are required for a differentiable difference between the two models.
      This process is performed for twice, once with an optimised design and
      once with a sub-optimal design, to illustrate how the optimised design
      reduces the time to confidently discern the magnetic moment.

    Args:
        save_path (str): path to directory to save results to.

    """
    # Define the range of times to consider (1 to 100 hours here).
    times = np.linspace(40, 4000, 250)

    # Points and split of total counting time for each angle to simulate.
    angle_splits = [(0.5, 100, 1 / 7), (1.0, 100, 2 / 7), (2.0, 100, 4 / 7)]

    # Calculate log ratio of likelihoods with optimal and sub-optimal designs.
    _calc_log_ratios(26, times, angle_splits, save_path)
    _calc_log_ratios(80, times, angle_splits, save_path)

    # Create the plot of counting time versus log ratio of likelihoods.
    fig = plt.figure(figsize=(6, 7))
    ax = fig.add_subplot(111)

    # Iterate over the results for the two designs.
    save_path = os.path.join(save_path, 'YIG_sample')
    for pt_thick in [26, 80]:
        file_path = os.path.join(
            save_path, 'log_ratios_{}.csv'.format(pt_thick)
        )

        # Load and plot the calculated ratios.
        data = np.loadtxt(file_path, delimiter=',')
        times, ratios = data[:, 0], data[:, 1]

        ax.plot(1.5 * times, ratios, label='{}Å Pt Thickness'.format(pt_thick))

    ax.set_xlabel('Counting Time (min.)', fontsize=11, weight='bold')
    ax.set_ylabel('Log Ratio of Likelihoods', fontsize=11, weight='bold')
    ax.legend()

    # Save the plot.
    save_plot(fig, save_path, 'log_ratios')


def _calc_log_ratios(pt_thick, times, angle_splits, save_path):
    """Calculates log ratio of likelihoods between two models, one with an
      induced moment in the YIG sample Pt layer and one with no moment, as
      a function of measurement time.

    Args:
        pt_thick (float): Pt layer thickness to use.
        times (numpy.ndarray): times to calculate the ratios over.
        angle_splits (list): points and time splits for each angle to simulate.
        save_path (str): path to directory to save results to.

    """
    # Create a .txt file to save the results to.
    save_path = os.path.join(save_path, 'YIG_sample')
    file_path = os.path.join(save_path, 'log_ratios_{}.csv'.format(pt_thick))
    with open(file_path, 'w') as file:
        # Iterate over each time being considered.
        for total_time in times:
            ratios = []
            # Get the ratio for 100 simulated data sets using the time.
            for _ in range(100):
                # Define the number of points and times for each angle.
                """angle_times = [(angle, points, split * total_time)
                    for angle, points, split in angle_splits]"""

                # Simulate data for the YIG sample with a 1 uB/atom magnetic
                # moment in the Pt layer.
                sample = SampleYIG()
                sample.pt_mag.value = 0.01638

                # structure = sample.using_conditions(pt_thick=pt_thick)
                models = []
                """simulate_magnetic(structure, angle_times,
                                                     scale=1, bkg=5e-7, dq=2,
                                                     pp=True, pm=False,
                                                     mp=False, mm=True)
                """
                # Calculate the log-likelihood of a model containing the
                # Pt layer magnetic moment.
                logl_1 = _logl(models)

                # Calculate the log-likelihood of a model without the
                # Pt layer magnetic moment.
                sample.pt_mag.value = 0
                logl_2 = _logl(models)

                # Record the ratio of likelihoods.
                ratio = logl_1 - logl_2
                ratios.append(ratio)

            # Calculate and save median ratio.
            median_ratio = np.median(ratios)
            file.write('{0},{1}\n'.format(total_time, median_ratio))
            print(median_ratio)


def _logl(models):
    """Calculates the log-likelihood for a given list of `models`
       corresponding to simulated spin states.

    Args:
        models (list): models to calculate log-likelihood for.

    Returns:
        float: log-likelihood of given `models`.

    """
    # Extract the Q, R, dR and model R for the simulated spin states.
    q, r, dr, r_model = [], [], [], []
    for model in models:
        probe = model.probe.xs[model.probe.spin_state]
        q.append(probe.Q)
        r.append(probe.R)
        dr.append(probe.dR)
        r_model.append([])  # reflectivity(probe.Q, model))

    # Combine the data from each spin state.
    r = np.concatenate(r)
    dr = np.concatenate(dr)
    r_model = np.concatenate(r_model)

    # Calculate the log-likelihood over all the models.
    return -0.5 * np.sum(((r - r_model) / dr) ** 2 + np.log(2 * np.pi * dr**2))


if __name__ == '__main__':
    save_path = './results'

    _magnetism_results_visualise(save_path)
    _magnetism_results_optimise(save_path)
    _magnetism_results_ratios(save_path)
