				#####################
				#		            #
				# BlackHawk scripts #
				#		            #
				#####################



The scripts in this folder are designed to compute the greybody factors and their
asymptotical fits.

author : Alexandre Arbey, alexandre.arbey@ens-lyon.fr & Jérémy Auffinger, jeremy.auffinger.blackhawk@protonmail.com
last modified : 07 April 2024

#########################################################################################

This folder contains 4 types of files:

- spin_*.m are Mathematica packages used to compute the greybody factors for particles of
spin *

- exploitation.nb is a Mathematica script used to obtain the fitting parameters for Kerr BHs

- formatting.c is a script formatting the tables

- Makefile is a compilation file

The 3 subfolders /charged, /higher and /LQG contain equivalent scripts for the associated
modified metrics. The formatting and exploitation scripts are fused together and are coded
in Python.

#########################################################################################

To use the spin_*.m packages, you will need Mathematica (version 11 or later). To obtain
this, ask your informatics department for a University license, or use whatever way you
want. Once you have installed Mathematica, you can directly launch the .m packages by
typing:

	math -script ./spin_*.m

We recommand that you add a "> output" at the end of this command due to the extensive
debug output of these packages.
These scripts will generate test_*_fM.txt and test_*_gM.txt files. In these files, each
line corresponds to a value of a* and each column to a value of x = 2*E*M.

#########################################################################################

The parameters of these scripts are the following:

- nbmodes is the number of values of l that will be used

- nbener is the number os values of x that will be tabulated. As M = 0.5, it corresponds
to the number of energies E.

- nba, nbe, nbn are the number of secondary parameter values that will be tabulated

Do not modify the constitution of the tables as well as their boundaries unless
you are sure to handle the subtelties of the numerical computation at stake. Delicate
equilibrium in the special cases was attained with these parameters. We recommand the
user to be carefull about the parameters "prox" and "far" that determine at what point
(close to the horizon) we start the integration, and at what point (at infinity) we end it.
It may be wise to choose at least far = 300/x to obtain goo numerical convergence.

#########################################################################################

Carefull : some manual checking may be useful in order to see if the computed greybody
factors are continuous and do not present spurrious pikes/values (numerical errors).
Once you have obtained the test_*_fM.txt and test_*_gM.txt tables, you can use the script
exploitation.nb in order to numerically compute the parameters of the low- and high-
energy asymptotical fits. Launch this script as a normal Mathematica notebook (information
may be found on the Mathematica web page).
You should make sure that you launch it in the same directory as the greybody factors
tables.
A first part is devoted to the computation of the high-energy fits. You may change the
parameters nbener and nba accordingly with the .m packages.
The files read by the ReadList command should be the test_*_fM.txt and test_gM_.txt
computed before, one after another.
Do not change the details of the fit parameters finding unless you are sure of what you
do. If you want to check that the fit is good, uncomment the Plot lines and plot the
computed values vs. the fitted ones.
A second part is devoted to the computation of the low-energy fits. Do not change the
details of the fit parameters unless you are sure of what you do. If you want to check
that the fit is good, uncomment the Plot lines and plot the computed values vs. the
fitted ones.
The name of the output file in which OpenWrite writes should be of the form fits*_fM.txt
or fits*_gM.txt.

#########################################################################################

Once you have computed both the values of the Kerr greybody factors and their asymptotical
fits, you may use the formatting.c script to turn them into a BlackHawk readable format.
First compile the formatting.c script by typing:

	make

This will generate an executable file formatting.x. Please make sure that the nb_a and
nb_x parameters are set correspondingly with those of the already computed tables. If
you have changed the format of the x table or the a* table, make sure you modify this
script as well. we recommand that you launch this executable in the same folder than
the previous one, or to link or copy the Mathematica tables to it. Launch by typing:

	./formating.x

This will create 8 new files spin_*_fM.txt, spin_*_gM.txt, spin_*_fits_fM.txt and
spin_*_fits_gM.txt that BlackHawk can read.

#########################################################################################

If you have any issue using these scripts please feel free to contact the author.








