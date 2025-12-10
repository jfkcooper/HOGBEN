"""
Functions for visualizing and optimizing angle choice for a given sample.
"""

import os
import time

import matplotlib.pyplot as plt
import numpy as np

from hogben.optimise import Optimiser
from hogben.visualise import contrast_choice_single, contrast_choice_double
from hogben.utils import save_plot
from models.bilayers import BilayerDMPC, BilayerDPPC


def _contrast_results_visualise(save_path):
    """Visualises the choice of contrasts for a bilayer sample.

    Args:
        save_path (str): path to directory to save results to.

    """

    # Choose sample here.
    bilayer = BilayerDMPC()

    # Number of points and counting times for each angle to simulate.
    angle_times = [(0.7, 100, 10), (2.3, 100, 40)]

    # Visualise single contrast choices assuming different initial
    # measurements.
    contrast_range = np.linspace(-0.56, 6.36, 500)
    contrast_choice_single(bilayer, contrast_range, [], angle_times,
                           save_path, 'initial')
    contrast_choice_single(bilayer, contrast_range, [6.36], angle_times,
                           save_path, 'D2O')
    contrast_choice_single(bilayer, contrast_range, [-0.56], angle_times,
                           save_path, 'H2O')
    contrast_choice_single(bilayer, contrast_range, [-0.56, 6.36], angle_times,
                           save_path, 'D2O_H2O')

    # Investigate contrast pair choices assuming no prior measurement.
    contrast_range = np.linspace(-0.56, 6.36, 75)
    contrast_choice_double(bilayer, contrast_range, angle_times, save_path)

    # Run nested sampling on simulated data.
    bilayer.nested_sampling([6.36, 6.36], angle_times, save_path, 'D2O_D2O')
    bilayer.nested_sampling([-0.56, 6.36], angle_times, save_path, 'H2O_D2O')


def _contrast_results_optimise(save_path):
    """Optimises the choice of contrasts for a bilayer sample.

    Args:
        save_path (str): path to directory to save results to.

    """
    # Choose sample here.
    bilayer = BilayerDMPC()

    # Time budget for experiment.
    total_time = 1000  # A large time improves DE convergence.

    # Points and proportion of total counting time for each angle.
    angle_splits = [(0.7, 100, 0.2), (2.3, 100, 0.8)]

    # Interval containing contrast SLDs to consider.
    contrast_bounds = (-0.56, 6.36)

    # Create a new .txt file for the results.
    save_path = os.path.join(save_path, bilayer.name)
    with open(os.path.join(save_path, 'optimised_contrasts.txt'), 'w') as file:
        optimiser = Optimiser(bilayer)  # Optimiser for the experiment.

        # Optimise the experiment using 1-4 contrasts.
        for i, num_contrasts in enumerate([1, 2, 3, 4]):
            # Display progress.
            print('>>> {0}/{1}'.format(i, 4))

            # Time how long the optimisation takes.
            start = time.time()
            results = optimiser.optimise_contrasts(num_contrasts, angle_splits,
                                                   total_time, contrast_bounds,
                                                   verbose=False)
            end = time.time()

            contrasts, splits, val = results
            splits = np.array(splits) * 100  # Convert to percentages.

            # Round the optimisation function value to 4 significant figures.
            val = np.format_float_positional(val, precision=4, unique=False,
                                             fractional=False, trim='k')

            # Write the conditions, objective value and computation time.
            file.write(
                '----------- {} Contrasts -----------\n'.format(num_contrasts))
            file.write('Contrasts: {}\n'.format(list(np.round(contrasts, 2))))
            file.write('Splits (%): {}\n'.format(list(np.round(splits, 1))))
            file.write('Objective value: {}\n'.format(val))
            file.write(
                'Computation time: {}\n\n'.format(
                    round(
                        end - start, 1)))


def _figure_2(save_path):
    """Creates figure 2 for the paper. The figure shows the minimum eigenvalue
       versus contrast SLD for the DMPC and DPPC/Ra LPS bilayer models, for a
       third contrast choice, assuming D2O and H2O were previously measured.

    Args:
        save_path (str): path to directory to save figure to.

    """
    # Load the two bilayer models.
    sample_1, sample_2 = BilayerDMPC(), BilayerDPPC()

    # Define the total time and how this time is split.
    total_time = 1000
    angle_splits = [(0.7, 100, 0.2), (2.3, 100, 0.8)]

    # Define how the counting time is split between the H2O and D2O contrasts.
    d2o_split_1 = 0.515
    d2o_split_2 = 0.540
    h2o_split_1 = 0.220
    h2o_split_2 = 0.245

    # Define how much the third contrast will be measured for.
    nxt_split_1 = 1 - d2o_split_1 - h2o_split_1
    nxt_split_2 = 1 - d2o_split_2 - h2o_split_2

    # Counting times for each angle when measuring D2O.
    d2o_angle_times_1 = [(angle, points, total_time * split * d2o_split_1)
                         for angle, points, split in angle_splits]

    d2o_angle_times_2 = [(angle, points, total_time * split * d2o_split_2)
                         for angle, points, split in angle_splits]

    # Counting times for each angle when measuring H2O.
    h2o_angle_times_1 = [(angle, points, total_time * split * h2o_split_1)
                         for angle, points, split in angle_splits]

    h2o_angle_times_2 = [(angle, points, total_time * split * h2o_split_2)
                         for angle, points, split in angle_splits]

    # Counting times for each angle when measuring the third contrast.
    nxt_angle_times_1 = [(angle, points, total_time * split * (nxt_split_1))
                         for angle, points, split in angle_splits]

    nxt_angle_times_2 = [(angle, points, total_time * split * (nxt_split_2))
                         for angle, points, split in angle_splits]

    # Information from the D2O and H2O contrasts for each model.
    g_init_1 = \
        (sample_1.angle_info(d2o_angle_times_1, [6.36]).fisher_information
         + sample_1.angle_info(h2o_angle_times_1, [-0.56]).fisher_information)

    g_init_2 = \
        (sample_2.angle_info(d2o_angle_times_2, [6.36]).fisher_information
         + sample_2.angle_info(h2o_angle_times_2, [-0.56]).fisher_information)
    # Calculate the minimum eigenvalue for each choice of third contrast.
    min_eigs_1, min_eigs_2 = [], []
    contrast_range = np.linspace(-0.56, 6.36, 500)
    for new_contrast in contrast_range:
        # Combine Fisher information matrices from initial and third contrasts.
        g_new_1 = sample_1.contrast_info(
            nxt_angle_times_1, [new_contrast]).fisher_information
        g_new_2 = sample_2.contrast_info(
            nxt_angle_times_2, [new_contrast]).fisher_information
        min_eigs_1.append(np.linalg.eigvalsh(g_init_1 + g_new_1)[0])
        min_eigs_2.append(np.linalg.eigvalsh(g_init_2 + g_new_2)[0])

    # Create the plot of third contrast versus minimum eigenvalue.
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()

    # Plot the data for the two models on separate y-axes.
    line_1 = ax1.plot(contrast_range, min_eigs_1, color='b')
    line_2 = ax2.plot(contrast_range, min_eigs_2, color='g')

    x_label = r'$\mathregular{Contrast\ SLD\ (10^{-6} \AA^{-2})}$'
    y_label = 'Minimum Eigenvalue'

    ax1.set_xlabel(x_label, fontsize=11, weight='bold')
    ax1.set_ylabel(y_label, fontsize=11, weight='bold', color='b')
    ax2.set_ylabel(y_label, fontsize=11, weight='bold', color='g')
    ax1.legend(line_1 + line_2, ['DMPC', 'DPPC/RaLPS'], loc=0)

    # Save the plot.
    save_plot(fig, save_path, 'figure_2')


if __name__ == '__main__':
    save_path = './results'
    _contrast_results_visualise(save_path)
    _contrast_results_optimise(save_path)

    _figure_2('../figures')
