				#####################
				#		    		#
				# BlackHawk scripts #
				#		    		#
				#####################



The scripts in this folder are designed to plot data computed by BlackHawk.

author: Jérémy Auffinger, jeremy.auffinger.blackhawk@protonmail.com
last modified: 15 September 2023

#########################################################################################

This folder contains 2 files:

- plot_inst.py to plot BlackHawk_inst data

- plot_tot.py to plot Black_hawk_tot data

#########################################################################################

To use the .py scripts, you need an installed version of Python (version 3.4.1 used for
the scripts), see here:

	https://www.python.org/

Follow the instructions of the Python webpage for the installation. Then, open the .py
scripts and launch the different blocks to plot data.

#########################################################################################

The script plot_inst.py contains 8 blocks:

1) Necessary importations: this block contains the commands importing the plotting and
math environments for Python. It also defines some geometric options for the plotting
window, which you can change at will. Launch it before all the rest.
2) Folder definition: this block defines the paths toward the computed data and the
localization of your favorite figure folder.
3) Recovering data: this block reads the BlackHawk data into arrays.
4) Plotting BH spectrum: this block plots the BHs mass funciton, that is to say the
comoving density as a function of mass. The mass bins are represented by extended
horizontal bars (this may fail for a Dirac distribution).
5) Plotting options (primary data): this block defines which primary particle spectra
will be plotted.
6) Plotting primary spectra: this block plots the desired primary particles spectra, that
is to say their instantaneous emission as a function of energy.
7) Plotting options (secondary data): this block defines which secondary particle spectra
will be plotted.
8) Plotting secondary spectra: this block plots the desired secondary particles spectra,
that is to say their instantaneous emission as a function of energy. The variable "epoch"
chooses between the BBN and present epoch (they have different secondary particles).

#########################################################################################

The script plot_tot.py contains 10 blocks:

1) Necessary importations: this block contains the commands importing the plotting and
math environments for Python. It also defines some geometric options for the plotting
window, which you can change at will. Launch it before all the rest.
2) Folder definition: this block defines the paths toward the computed data and the
localization of your favorite figure folder.
3) Recovering spectrum data: this block reads the BHs mass function and stores it into an
array.
4) Plotting BH spectrum: this block plots the BHs mass funciton, that is to say the
comoving density as a function of mass. The mass bins are represented by extended
horizontal bars (this may fail for a Dirac distribution).
5) Plotting options (primary data): this block defines which primary particle spectra
will be plotted. Please make sure that the corresponding "primary_spectra" have been
written by BlackHawk.
6) Plotting primary spectra (fixed time): this block plots the desired primary particles
spectra at a fixed time, that is to say their instantaneous emission as a function of
energy.
7) Plotting primary spectra (fixed energy): this block plots the desired primary particles
spectra at a fixed energy, that is to say their time-dependant emission.
8) Plotting options (secondary data): this block defines which secondary particle spectra
will be plotted.
9) Plotting secondary spectra (fixed time): this block plots the desired secondary
particles spectra at a fixed time, that is to say their instantaneous emission as a
function of energy.
10) Plotting secondary spectra (fixed energy): this block plots the desired secondary
particles spectra at a fixed energy, that is to say their time-dependent emission.

#########################################################################################

If you have any issue using these scripts please feel free to contact the author.























