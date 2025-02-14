				#####################
				#		            #
				# BlackHawk scripts #
				#		            #
				#####################



The scripts in this folder are designed to compute the HERWIG hadronization tables
for final particles corresponding to BBN epoch (lifetime > 10^(-8) s).

authors : Alexandre Arbey, alexandre.arbey@ens-lyon.fr & Jérémy Auffinger, jeremy.auffinger.blackhawk@protonmail.com
last modified : 07 April 2024

#########################################################################################

This folder contains 4 types of files and folders:

- /* are folders containing scripts for each initial particle

- BH_TABLE_GENERATOR.cc is a file used to perform a Rivet analysis

- formatting.c is a C file used to format the computed tables

- Makefile is a compilation file

In each one of the subfolders, there are two files:

- main.cpp is a C++ file used to launch HERWIG

- LEP.in is an input file for HERWIG

#########################################################################################

To use the main.cpp scripts, you will need HERWIG (version 7 or later), that you can find
here:

	https://herwig.hepforge.org/downloads.html

Follow the instructions of the HERWIG webpage for the installation. You also need to link
the Rivet analysis to HERWIG, and for this step we refer the user to the HERWIG tutorials.
If you want to use the bootstrap installation method, you can simply follow the following
steps:

mkdir herwig
cd herwig
wget https://herwig.hepforge.org/downloads/herwig-bootstrap
chmod +x herwig-bootstrap
./herwig-bootstrap .
cd ..

To activate the HERWIG environment, type:

	source herwig/bin/activate

Build and link the Rivet shared library with:

rivet-buildplugin RivetBH_TABLE_GENERATOR.so BH_TABLE_GENERATOR.cc
export RIVET_ANALYSIS_PATH=$PWD

Compile all the main.cpp scripts in the subfolders with:
g++ main.cpp -o main.x

This will create executable files main.x.

To launch the scripts, type:

	./main.x

We recommand to add a "> output" to this command because of the extensive HERWIG output.
Several scripts can be launched in parallel.

#########################################################################################

These scripts will create numerous folders in the particle subfolders (one per final
energy). In each of these sub_folder, a version of LEP.in will be copied with its beam
energy fixed, and used by HERWIG to compute the collisions. this will create a LEP.yoda
file containing all the collision informations.

#########################################################################################

To be usable by BlackHawk, these LEP.yoda files need to be formatted in another way: that 
is the purpose of formatting.c. To use this script, simply return to the root directory
and compile the formatting.c script by typing:

	make

This will create an executable formatting.x. Launch it by typing:

	./formatting.x

This will read the yoda files and create the new *.txt tables.

#########################################################################################

The parameters of these scripts are the following:

- nb_initial_energies / nb_init_en : number of initial energies for the initial particles

- nb_final_energies / nb_fin_en : number of final energies for the final particles

- nb_initial_particles / nb_init_part : number of initial particle types

	14 : bottom, charm, down, electron, gluons, higgs, muon, photon, strange, tau, top,
	     up, W+-, Z0

- nb_final_particles / nb_fin_part : number of final particle types

	11 : photon, electron, muon, nu_e, nu_mu, nu_tau, pi+-, K+-, K0 long, proton,
	     neutron

- Emin_init / Emax_init / Emin / Emax : initial energy boundaries (in GeV)
	
	Emin_init should not be less than 25 GeV as it is a limit for HERWIG.
	Emax_init could be up to some 10^2 TeV.
	For more informations see the HERWIG docs.

- Emin_fin / Emax_fin : final energy boundaries (in GeV)

Please make sur that the parameters of the main.cpp and formating.c correspond in order to
avoid segmentation faults.

##########################################################################################

If you have any issue using these scripts please feel free to contact the authors.
