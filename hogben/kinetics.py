"""Functions for visualizing and optimizing angle and contrast for kinetic
experiments."""

import numpy as np
import os

import matplotlib.pyplot as plt

from scipy.optimize import differential_evolution

from hogben.models.monolayers import MonolayerDPPG
from hogben.utils import save_plot

plt.rcParams['figure.figsize'] = (6, 4)
plt.rcParams['figure.dpi'] = 600


def _kinetics_results_visualise(save_path):
    """Visualises the choice of measurement angle and contrast SLD for the
       DPPG monolayer model where the lipid area per molecule is decreasing
       over time.

    Args:
        save_path (str): path to directory to save results to.

    """
    # Choose whether the monolayer tailgroups are hydrogenated or deuterated.
    monolayer = MonolayerDPPG(deuterated=False)

    # Angles, contrasts and lipid area per molecule values to consider.
    angle_range = np.linspace(0.2, 4.0, 15)
    contrast_range = np.linspace(-0.56, 6.36, 15)
    apm_range = np.linspace(54.1039, 500, 20)

    # Number of points and measurement time to use when simulating data.
    points = 100
    time = 100

    # Iterate over each contrast and angle being considered.
    x, y, infos = [], [], []
    n = len(angle_range) * len(contrast_range)  # Number of calculations.
    for i, contrast_sld in enumerate(contrast_range):
        # Display progress.
        if i % 5 == 0:
            print('>>> {0}/{1}'.format(i * len(angle_range), n))

        for angle in angle_range:
            # Record the "true" lipid APM value.
            apm = monolayer.lipid_apm.value

            # Split the time budget based on number of APM values.
            angle_times = [(angle, points, time / len(apm_range))]

            # Calculate the lipid APM Fisher information for each value.
            information = 0
            for new_apm in apm_range:
                # Use the new APM value, corresponding to a degrading sample.
                monolayer.lipid_apm.value = new_apm
                g = monolayer.contrast_info(angle_times, [contrast_sld])
                information += g[0, 0]

            monolayer.lipid_apm.value = apm  # Reset the APM parameter.
            infos.append(information)
            x.append(contrast_sld)
            y.append(angle)

    # Create the plot of angle and contrast SLD versus Fisher information.
    fig = plt.figure(figsize=[10, 8])
    ax = fig.add_subplot(111, projection='3d')

    # Create the surface plot and add colour bar.
    surface = ax.plot_trisurf(x, y, infos, cmap='plasma')
    fig.colorbar(surface, fraction=0.046, pad=0.04)

    x_label = r'$\mathregular{Contrast\ SLD\ (10^{-6} \AA^{-2})}$'
    y_label = 'Angle (°)'
    z_label = 'Fisher Information'

    ax.set_xlabel(x_label, fontsize=11, weight='bold')
    ax.set_ylabel(y_label, fontsize=11, weight='bold')
    ax.set_zlabel(z_label, fontsize=11, weight='bold')
    ax.ticklabel_format(axis='z', style='sci', scilimits=(0, 0))

    # Save the plot.
    save_path = os.path.join(save_path, monolayer.name)
    tg_type = 'deuterated' if monolayer.deuterated else 'hydrogenated'
    filename = 'angle_contrast_choice_' + tg_type
    save_plot(fig, save_path, filename)

    # Save different views of the 3D plot.
    # Iterate from 0 to 360 degrees in increments of 10.
    save_path = os.path.join(save_path, 'angle_contrast_choice', tg_type)
    for i in range(0, 360, 120):
        # Set the "camera" view for the 3D plot.
        ax.view_init(elev=40, azim=i)

        # Save the view.
        filename_i = filename + '_' + str(i)
        save_plot(fig, save_path, filename_i)

    # Return the contrast and angle with greatest Fisher information.
    maximum = np.argmax(infos)
    return x[maximum], y[maximum]


def _kinetics_results_optimise(save_path):
    """Optimises the choice of measurement angle and contrast SLD for the
       DPPG monolayer model where the lipid area per molecule is decreasing
       over time.

    Args:
        save_path (str): path to directory to save results to.

    """
    # Counting time and number of points to use when simulating.
    time = 1000  # A large time improves DE convergence.
    points = 100

    # Intervals containing angles and contrasts to consider.
    angle_bounds = (0.2, 4.0)
    contrast_bounds = (-0.56, 6.36)

    # Lipid area per molecule values to calculate with.
    apm_range = np.linspace(54.1039, 500, 20)

    # Bounds for optimisation function variables.
    bounds = [angle_bounds, contrast_bounds]
    # Arguments for optimisation function.
    args = [apm_range, points, time]

    # Optimise both the hydrogenated and deuterated tailgroup models.
    file_path = os.path.join(save_path, 'DPPG_monolayer', 'optimised.txt')
    with open(file_path, 'w') as file:
        for i, label in enumerate(['Hydrogenated', 'Deuterated']):
            monolayer = MonolayerDPPG(deuterated=bool(i))

            # Optimise angle and contrast.
            res = differential_evolution(_optimisation_func, bounds,
                                         args=[monolayer] + args, polish=False,
                                         tol=0.001, updating='deferred',
                                         workers=-1, disp=False)

            angle, contrast = res.x[0], res.x[1]

            # Write the angle, contrast and function value to the .txt file.
            file.write('-------------- ' + label + ' --------------\n')
            file.write('Angle: {}\n'.format(round(angle, 2)))
            file.write('Contrast: {}\n'.format(round(contrast, 2)))
            file.write('Objective value: {}\n\n'.format(res.fun))


def _optimisation_func(x, sample, apm_range, points, time):
    """Defines the function for optimising the DPPG monolayer measurement
       angle and contrast.

    Args:
        x (list): angle and contrast to calculate the function with.
        sample (monolayers.MonolayerDPPG): monolayer model being optimised.
        apm_range (numpy.ndarray): lipid area per molecule values.
        points (int): number of data points to use for each angle.
        time (float): total time budget for experiment.

    Returns:
        float: negative of the Fisher information in the lipid area per
               molecule parameter using the conditions `x`.

    """
    # Define the points and counting times for each angle to simulate.
    angle, contrast_sld = x[0], x[1]
    angle_times = [(angle, points, time / len(apm_range))]

    # Calculate the information in the experiment using the given conditions.
    information = 0
    apm = sample.lipid_apm.value
    for new_apm in apm_range:
        sample.lipid_apm.value = new_apm
        information += sample.contrast_info(angle_times, [contrast_sld])[0, 0]

    # Return negative of the Fisher information as the algorithm is minimising.
    sample.lipid_apm.value = apm
    return -information


if __name__ == '__main__':
    save_path = './results'

    _kinetics_results_visualise(save_path)
    _kinetics_results_optimise(save_path)
