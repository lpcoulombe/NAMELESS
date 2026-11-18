# NAMELESS
Data reduction pipeline for JWST time series observations. 

This framework is currently only adapted for NIRISS SOSS observations, with plans to expand it to other instrument modes in the future.

If you use this package in your work, please cite Coulombe et al. 2023 (ui.adsabs.harvard.edu/abs/2023Natur.620..292C/abstract) and Coulombe et al. 2025a (https://scixplorer.org/abs/2025NatAs...9..512C/abstract).

Installation instructions:
  1. Clone this directory and create your personal branch
  2. Create a conda environment in which you will run your NAMELESS reductions. The pipeline has been tested on Python 3.13.
  3. pip install setuptools.
  4. In the NAMELESS/ directory, run ./install.sh in your terminal to install the package.

How to use it:
  1. In the NAMELESS/Observations/ directory, create a folder that will contain the _uncal.fits files for your observations.
  2. In the NAMELESS/ folder, copy the example_config.yml file and create one for the specific reduction you want to do.
  3. Once you're done playing with the .yml file, save it and run in the terminal: python run_nameless.py yourplanet_config.yml
  4. All settings in the .yml file are explained, and feel free to vary them to see how they impact the reduction!

The end result of a full reduction is a .fits file containing the 2D stellar spectra and other relevant outputs.

If you encounter any issues using the code, please raise them via GitHub. 

Similarly, if you implement new functionalities while working with the pipeline. Please feel free to create pull requests so that I can implement them in the main branch.

Happy "reductioning"! 
- LP
