				#####################
				#		            #
				# BlackHawk scripts #
				#		            #
				#####################



The script in this folder is designed to compute the stacked redshifted spectrum of
Hawking radiated particles at the time of evaporation of the considered BH.

authors : Jérémy Auffinger, jeremy.auffinger.blackhawk@protonmail.com
last modified : 07 April 2024

#########################################################################################

This folder contains 2 files:

- stack.c is the C file used to compute the stacked spectrum

- Makefile is a compilation file

#########################################################################################

To use the C script, simply execute the command:

	make

This will create a executable stack.x. Then, to use the executable, execute:

	./stack.x PATH nb_init_ener domination CMB today

where:

- PATH is the path to your BlackHawk time-dependent results, e.g. ../BlackHawk/results/test/
if you ran BlackHawk with an destination_folder = test

- nb_init_ener is the number of energies in the spectrum considered (primary or secondary)

- domination determines the cosmology: 0 for full radiation domination, 1 for full matter domination, 2 for standard alternate domination of matter and radiation with exchange at M-D equality

- CMB determines whether the BH emission before CMB is taken into account

- today determines whether the spectrum is computed at full BH evaporation or at present time

#########################################################################################

This script also includes several other (built_in) parameters:

- nb_fin_ener = 500 is the number of energies in the outcome of this script. Try to keep the
number of final energies lower than the number of initial energies for good interpolation.

- Emin and Emax are the final energies used in the computation. Emax should always be around
the Planck energy as the final stages of evaporation produce very energetic particles. Emin
should be somewhat lower than the minimum energy in the initial spectra as this energy will
be diluted by the redshift.

#########################################################################################

Please note that by default this script computes the stacked graviton spectrum.
This can be changed by replacing the char array "name_spectrum" at line 112, by any other
spectrum you would like to stack.

#########################################################################################

The outcome of this script comes as a file "results_*.txt" in which the first column if the
energy E and the second column is the numder density dn/dE of particles. The result is
given in cm^(-3).GeV^(-1). The * is either RD, MD or alternate if the domination parameter is set to
0 (radiation domination), 1 (matter domination) or 2 (standard cosmology).
