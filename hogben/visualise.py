"""Methods to visualise the optimised experimental conditions using a saved
plot"""

import os
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.animation import FuncAnimation, PillowWriter
from itertools import combinations

from hogben.utils import save_plot, Fisher
from hogben.models.base import (
    BaseLipid,
    BaseSample,
    VariableAngle,
    VariableContrast,
    VariableUnderlayer,
)

plt.rcParams['figure.figsize'] = (9, 7)
plt.rcParams['figure.dpi'] = 600


def angle_choice(
    sample: BaseSample,
    initial_angle_times: list,
    angle_range: np.ndarray,
    points_new: int,
    time_new: type,
    save_path: str,
    filename: str,
    contrasts: Optional[list] = None,
):
    """Plots the minimum eigenvalue of the Fisher information matrix
       as a function of measurement angle.

    Args:
        sample (base.BaseSample): sample to investigate the angle choice for.
        initial_angle_times (list): points and times for each measured angle.
        angle_range (numpy.ndarray): range of angles to investigate.
        points_new (int): number of points to use for the new data.
        time_new (type): counting time to use for the new data.
        save_path (str): path to directory to save plot to.
        filename (str): filename to use when saving the plot.
        contrasts (list): contrasts to simulate, if applicable.

    Returns:
        float: angle with the largest minimum eigenvalue.

    """
    # Check that the angle can be varied for the sample.
    assert isinstance(sample, VariableAngle)

    # Set contrasts to empty list if not provided
    contrasts = [] if contrasts is None else contrasts

    # Iterate over each angle to consider.
    min_eigs = []
    for i, angle_new in enumerate(angle_range):
        # Display progress.
        if i % 100 == 0:
            print('>>> {0}/{1}'.format(i, len(angle_range)))

        # Get the information for the new angle.
        new_angle_times = \
            (initial_angle_times + [(angle_new, points_new, time_new)])
        # Calculate the total Fisher information at the new angle together
        # with the initial angle
        fisher_new = sample.angle_info(new_angle_times, contrasts)
        min_eigs.append(fisher_new.min_eigenval)

    # Plot measurement angle versus minimum eigenvalue.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(angle_range, min_eigs)
    ax.set_xlabel('Angle (°)', fontsize=11, weight='bold')
    ax.set_ylabel('Minimum Eigenvalue', fontsize=11, weight='bold')

    # Save the plot.
    save_path = os.path.join(save_path, sample.name)
    save_plot(fig, save_path, 'angle_choice_' + filename)

    # Return the angle with largest minimum eigenvalue.
    return angle_range[np.argmax(min_eigs)]


def scan_parameters(sample: BaseSample,
                    params: list,
                    angle_times: list) -> None:
    """
    Scans parameter values over FI and generate graphs.

    Args:
        sample (object): The sample data.
        params (list): List of parameters to scan.
        angle_times (list): List of angle times.

    Returns:
        None
    """
    param_values, eigenval_list = [], []
    for param in params:
        lb, ub = param.bounds.lb, param.bounds.ub
        old_value = param.value
        param_range = np.linspace(lb, ub, 100)
        eigenvals = []
        for value in param_range:
            param.value = value
            eigenvals.append(
                Fisher.from_sample(sample, angle_times).min_eigenval)
        param_values.append(param_range)
        eigenval_list.append(eigenvals)
        param.value = old_value
    ymax = np.max(eigenval_list)

    # Plot a scan for each parameter
    for i, param in enumerate(params):
        fig, ax = plt.subplots()

        ax.set_ylabel('Minimum eigenvalue')
        ax.set_xlabel('Parameter value')

        ax.set_ylim(0, 1.05 * ymax)
        ax.set_xlim(min(param_values[i]), max(param_values[i]))

        ax.set_title(param.name)
        ax.plot(param_values[i], eigenval_list[i], label=param.name)

def angle_choice_with_time(
    sample: BaseSample,
    initial_angle: float,
    angle_range: type,
    time_range: type,
    points: int,
    new_time: float,
    save_path: str,
    contrasts: Optional[list] = None,
):
    """Investigates how the second choice of angle for a `sample` changes
       as the counting time of the initial angle is increased.

    Args:
        sample (base.BaseSample): sample to investigate the angle choice for.
        initial_angle (float): angle initially measured.
        angle_range (type): range of angles to investigate.
        time_range (type): range of times to investigate.
        points (int): number of points to simulate for each angle.
        new_time (float): counting time for the second angle.
        save_path (str): path to directory to save plot to.
        contrasts (list): contrasts to simulate, if applicable.

    """
    # Check that the angle can be varied for the sample.
    assert isinstance(sample, VariableAngle)

    # Set contrasts to empty list if not provided
    contrasts = [] if contrasts is None else contrasts

    # Create plot of angle versus minimum eigenvalue.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Angle (°)', fontsize=11, weight='bold')
    ax.set_ylabel('Minimum Eigenvalue', fontsize=11, weight='bold')
    ax.set_xlim(angle_range[0], angle_range[-1])

    # Create the line that will have data added to it.
    (line,) = ax.plot([], [], lw=3)

    def init():
        """Initialiser function for the line"""
        line.set_data([], [])
        return line,

    # Annimation function for the line.
    def animate(i):
        """Describes how the function should be animated"""
        # Display progress.
        if i % 5 == 0:
            print('>>> {0}/{1}'.format(i, len(time_range)))

        # Get the information from the initial angle at the current time.
        angle_times_init = [(initial_angle, points, time_range[i])]

        # Iterate over each angle to consider.
        min_eigs = []
        for new_angle in angle_range:
            # Combine the information from the first and second angles.
            angle_times_new = \
                angle_times_init + [(new_angle, points, new_time)]
            # Calculate the total Fisher information using the first and
            # second angle together
            fisher_new = sample.angle_info(angle_times_new, contrasts)
            min_eigs.append(fisher_new.min_eigenval)

        # Update the data of the line.
        ax.set_ylim(min(min_eigs), max(min_eigs))
        line.set_data(angle_range, min_eigs)

        return line,

    # Define the animation.
    anim_length = 4000  # Animation length in milliseconds.
    frames = len(time_range)  # Number of frames for animation.
    anim = FuncAnimation(fig, animate, init_func=init, blit=True,
                         frames=frames, interval=anim_length // frames)
    plt.close()

    # Save the animation as a .gif file.
    writergif = PillowWriter()
    save_path = os.path.join(save_path, sample.name,
                             'angle_choice_with_time.gif')

    anim.save(save_path, writer=writergif)
    return anim


def contrast_choice_single(sample: BaseLipid,
                           contrast_range: np.ndarray,
                           initial_contrasts: list,
                           angle_times: list,
                           save_path: str,
                           filename: str
                           ) -> float:
    """Plots the minimum eigenvalue of the Fisher information matrix
       as a function of a single contrast SLD.

    Args:
        sample (base.BaseLipid): sample to investigate the contrast choice for.
        contrast_range (numpy.ndarray): range of contrast SLDs to investigate.
        initial_contrasts (list): SLDs of contrasts already measured.
        angle_times (list): points and counting time of each angle to simulate.
        save_path (str): path to directory to save plot to.
        filename (str): filename to use when saving the plot.

    Returns:
        float: contrast with the largest minimum eigenvalue.

    """
    # Check that the contrast can be varied for the sample.
    assert isinstance(sample, VariableContrast)

    # Iterate over each contrast to consider.
    min_eigs = []
    for i, new_contrast in enumerate(contrast_range):
        # Display progress.
        if i % 100 == 0:
            print('>>> {0}/{1}'.format(i, len(contrast_range)))

        # Get the information from the new contrast, and calculate the total
        # Fisher information for the initial and new contrast together.
        new_contrast = initial_contrasts + [new_contrast]
        fisher_new = sample.contrast_info(angle_times, new_contrast)
        min_eigs.append(fisher_new.min_eigenval)

    # Plot contrast SLD versus minimum eigenvalue.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(contrast_range, min_eigs)
    x_label = '$\mathregular{Contrast\ SLD\ (10^{-6} \AA^{-2})}$'
    ax.set_xlabel(x_label, fontsize=11, weight='bold')
    ax.set_ylabel('Minimum Eigenvalue', fontsize=11, weight='bold')

    # Save the plot.
    save_path = os.path.join(save_path, sample.name)
    save_plot(fig, save_path, 'contrast_choice_single_' + filename)

    # Return the contrast SLD with largest minimum eigenvalue.
    return contrast_range[np.argmax(min_eigs)]


def contrast_choice_double(sample, contrast_range, angle_times, save_path):
    """Plots the minimum eigenvalue of the Fisher information matrix
       as a function of two contrast SLDs, assuming no prior measurement.

    Args:
        sample (base.BaseLipid): sample to investigate contrast choice for.
        contrast_range (numpy.ndarray): range of contrast SLDs to investigate.
        angle_times (list): points and counting time of each angle to simulate.
        save_path (str): path to directory to save plot to.

    Returns:
        tuple: contrast SLD pair with the largest minimum eigenvalue.

    """
    # Check that the contrast can be varied for the sample.
    assert isinstance(sample, VariableContrast)

    # Get all possible unique pairs of contrast SLDs.
    contrasts = np.asarray(list(combinations(contrast_range, 2)))

    # Iterate over each contrast pair to consider.
    min_eigs = []
    for i, contrast_pair in enumerate(contrasts):
        # Display progress.
        if i % 500 == 0:
            print('>>> {0}/{1}'.format(i, len(contrasts)))

        # Calculate the minimum eigenvalue of the Fisher information matrix.
        fisher = sample.contrast_info(angle_times, contrast_pair)
        min_eigs.append(fisher.min_eigenval)

    # Duplicate the data so that the plot is not half-empty (or half-full?).
    # This is valid since the choice of contrast is commutative.
    # I.e., the order in which contrasts are measured is not important.
    x = np.concatenate([contrasts[:, 0], contrasts[:, 1]])
    y = np.concatenate([contrasts[:, 1], contrasts[:, 0]])
    min_eigs.extend(min_eigs)

    # Create a 3D plot of contrast pair versus minimum eigenvalue.
    fig = plt.figure(figsize=[12, 9])
    ax = fig.add_subplot(111, projection='3d')

    # Create the surface plot and add colour bar.
    surface = ax.plot_trisurf(x, y, min_eigs, cmap='plasma')
    fig.colorbar(surface, fraction=0.046, pad=0.04)

    x_label = '$\mathregular{Contrast \ 1 \ SLD \ (10^{-6} \AA^{-2})}$'
    y_label = '$\mathregular{Contrast \ 2 \ SLD \ (10^{-6} \AA^{-2})}$'
    z_label = 'Minimum Eigenvalue'

    ax.set_xlabel(x_label, fontsize=11, weight='bold')
    ax.set_ylabel(y_label, fontsize=11, weight='bold')
    ax.set_zlabel(z_label, fontsize=11, weight='bold')
    ax.ticklabel_format(axis='z', style='sci', scilimits=(0, 0))

    # Save the plot.
    save_path = os.path.join(save_path, sample.name)
    save_plot(fig, save_path, 'contrast_choice_double')

    # Save different views of the 3D plot.
    # Iterate from 0 to 360 degrees in increments of 10.
    save_path = os.path.join(save_path, 'contrast_choice_double')
    for i in range(0, 360, 60):
        # Set the "camera" view for the 3D plot.
        ax.view_init(elev=40, azim=i)

        # Save the view.
        filename_i = 'contrast_choice_double_{}'.format(i)
        save_plot(fig, save_path, filename_i)

    # Return the contrast SLD pair with largest minimum eigenvalue.
    maximum = np.argmax(min_eigs)
    return x[maximum], y[maximum]


def underlayer_choice(
    sample: BaseLipid,
    thickness_range: np.ndarray,
    sld_range: np.ndarray,
    contrasts: list,
    angle_times: list,
    save_path: str,
    filename: str = '',
):
    """Plots the minimum eigenvalue of the Fisher information matrix
       as a function of a sample's underlayer thickness and SLD.

    Args:
        sample (base.BaseLipid): sample to investigate underlayer choice for.
        thickness_range (numpy.ndarray): range of underlayer thicknesses.
        sld_range (numpy.ndarray): range of underlayer SLDs.
        contrasts (list): SLDs of contrasts to simulate.
        angle_times (list): points and counting time of each angle to simulate.
        save_path (str): path to directory to save plot to.
        filename (str): file name to use when saving the plot.

    Returns:
        tuple: underlayer thickness and SLD with largest minimum eigenvalue.

    """
    # Check that the underlayer can be varied for the sample.
    assert isinstance(sample, VariableUnderlayer)

    # Iterate over each underlayer thickness to investigate.
    x, y, min_eigs = [], [], []
    n = len(thickness_range) * len(sld_range)  # Number of calculations.
    for i, thick in enumerate(thickness_range):
        # Display progress.
        if i % 5 == 0:
            print('>>> {0}/{1}'.format(i * len(sld_range), n))

        # Iterate over each underlayer SLD to investigate.
        for sld in sld_range:
            # Calculate the minimum eigenvalue.
            fisher = \
                sample.underlayer_info(angle_times, contrasts, [(thick, sld)])
            min_eigs.append(fisher.min_eigenval)
            x.append(thick)
            y.append(sld)

    # Create 3D plot of underlayer thickness and SLD versus minimum eigenvalue.
    fig = plt.figure(figsize=[10, 8])
    ax = fig.add_subplot(111, projection='3d')

    # Create the surface plot and add colour bar.
    surface = ax.plot_trisurf(x, y, min_eigs, cmap='plasma')
    fig.colorbar(surface, fraction=0.046, pad=0.04)

    x_label = '$\mathregular{Underlayer\ Thickness\ (\AA)}$'
    y_label = '$\mathregular{Underlayer\ SLD\ (10^{-6} \AA^{-2})}$'
    z_label = 'Minimum Eigenvalue'

    ax.set_xlabel(x_label, fontsize=11, weight='bold')
    ax.set_ylabel(y_label, fontsize=11, weight='bold')
    ax.set_zlabel(z_label, fontsize=11, weight='bold')
    ax.ticklabel_format(axis='z', style='sci', scilimits=(0, 0))

    # Save the plot.
    save_path = os.path.join(save_path, sample.name)
    filename = 'underlayer_choice' + ('_' + filename if filename != '' else '')
    save_plot(fig, save_path, filename)

    # Save different views of the 3D plot.
    # Iterate from 0 to 360 degrees in increments of 10.
    save_path = os.path.join(save_path, 'underlayer_choice', filename)
    for i in range(0, 360, 10):
        # Set the "camera" view for the 3D plot.
        ax.view_init(elev=40, azim=i)

        # Save the view.
        filename_i = filename + '_' + str(i)
        save_plot(fig, save_path, filename_i)

    # Return the underlayer thickness and SLD with largest minimum eigenvalue.
    maximum = np.argmax(min_eigs)
    return x[maximum], y[maximum]
