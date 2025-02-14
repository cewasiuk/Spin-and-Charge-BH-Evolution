				#####################
				#		            #
				# BlackHawk scripts #
				#		            #
				#####################



The scripts in this folder are designed to compute the f(M) and g(M) tables, for different
metrics.

Authors : Alexandre Arbey, alexandre.arbey@ens-lyon.fr & Jérémy Auffinger, jeremy.auffinger.blackhawk@protonmail.com
last modified : 07 April 2024

#########################################################################################

This folder contains 2 files:

- fM.c is the script computing the tables

- Makefile is a compilation file

#########################################################################################

To use the fM.c script, first compile it by typing:

	make

This will create an executable file fM.x. This executable needs the tables computed by
the "gamma_factors" scripts, so you have to unzip them from the archive tables.zip.
To launch this script, simply type:

	./fM.x type

where "type" switches between the different metrics. Choose 0 for the Kerr metric, this
will create two tables f(M,a*) and g(M,a*). Choose 1 for the polymerized metric, this
will create one table f(M,epsilon).

#########################################################################################

In these output files, the first column is the black holes mass in GeV. The following
columns are the values of f(M) (and g(M)) in GeV units for each a*/epsilon tabulated (see
first line).

#########################################################################################

The parameters of this script are the following:

- nb_particles = 15 the number of initial (SM) particles
- grav decides whether the graviton is taken into account (1) or not (0)
- ADD_DM decides whether dark matter is taken into account (1) or not (0)
- spin_DM is the DM spin (between 0, 1, 2, 0.5 and 1.5)
- m_DM is the DM mass in GeV
- dof_DM is the number of DM degrees of freedom
- nb_a is the number of tabulated a* (Kerr)
- nb_epsilon is the number of tabulated epsilon (polymerized)
- a0 is the minimal area in LQG (polymerized BHs)
- nb_x is the number of tabulated x = 2*E*M
- nb_particle_spins is the number of spins available for each metric (4 for the polymerized
BHs, 5 for the Kerr BHs)
- nb_masses is the number of masses for which the values will be computed
- Mmax is the maximum BH mass in GeV

Please make sur that the parameters of fM.c correspond to those of the "gamma_factors"
tables in order to avoid segmentation faults.

#########################################################################################

An additional feature is available: one can add degrees of freedom to the SM (+graviton).
These correspond to the lines marked by a "// ADD" comment in the different routines.

#########################################################################################

Be careful that the spin 3/2 tables have only been computed for the Schwarzschild metric
with a* = 0, not for the Kerr metric for a* != 0 nor for the polymerized metric.

#########################################################################################

We would like to stress on the fact that for a sufficiently small BH mass M such that the
horizon area gets close to the minimum area in LQG, the form of the metrics and the
paradigm of Hawking radiation may not hold.

#########################################################################################

If you have any issue using these scripts please feel free to contact the authors.








