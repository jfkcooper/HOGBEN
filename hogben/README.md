## Directories
* [data](/hogben/data) - Directbeam files and experimentally-measured data for a selection of samples of varying complexity.
* [models](/hogben/models) - Model definitions and fits for the aforementioned samples.
* [results](/hogben/results) - Results for the aforementioned samples.

## Code
* [angles.py](/hogben/angles.py) - Optimises and visualises the choice of measurement angle(s) for a collection of samples of varying complexity.
* [contrasts.py](/hogben/contrasts.py) - Optimises and visualises the choice of contrast for the [DMPC](/hogben/models/bilayers) and [DPPC/RaLPS](/hogben/models/bilayers) bilayer models.
* [kinetics.py](/hogben/kinetics.py) - Optimises and visualises the choice of measurement angle and contrast for the [DPPG](/hogben/models/monolayers) monolayer model degrading over time.
* [magnetism.py](/hogben/magnetism.py) - Optimises and visualises the sample design of the magnetic [YIG](/hogben/models/magnetic) sample.
* [optimise.py](/hogben/optimise.py) - Contains code for optimising the choice of measurement angle(s), counting time(s), contrast(s) and underlayer(s).
* [simulate.py](/hogben/simulate.py) - Contains code for simulating experiments using a [directbeam](/hogben/data/directbeams) file of incident neutron flux as a function of wavelength.
* [underlayers.py](/hogben/underlayers.py) - Optimises and visualises the choice of underlayer thickness(es) and SLD(s) for the [DMPC](/hogben/models/bilayers) and [DPPC/RaLPS](/hogben/models/bilayers) bilayer models.
* [utils.py](/hogben/utils.py) - Contains miscellaneous code for calculating the Fisher information, nested sampling, and saving plots.
* [visualise.py](/hogben/visualise.py) - Contains code for visualising the choice of measurement angle(s), counting time(s), contrast(s) and underlayer(s).
