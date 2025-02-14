				#####################
				#		    		#
				# BlackHawk scripts #
				#		    		#
				#####################

This script launches BlackHawk instances as many times as specified, with parameters fixed
by the user following a decision tree.
author : Jérémy Auffinger, j.auffinger@ipnl.in2p3.fr
last modified : 4 January 2021

#########################################################################################

This folder contains 3 files:

- a file BH_launcher.c that contains all the routines for the program BH_launcher

- a file Makefile to compute BH_launcher.c into BH_launcher.x

- a template file parameters.txt that contains blank values for all the BlackHawk parameters

#########################################################################################

Before compiling BH_launcher, please make sure that the path to BlackHawk is correct.
Before using it, please make sure that the relevant BlackHawk programs are already compiled.
If you have moved BH_laucnher out of its default position, then you may have to modify the
paths inside BH_launcher routines. To compile BH_launcher.c into BH_launcher.x, just type

make

into the command from this directory. To use BH_launcher.x, just type

./BH_launcher.x

without any argument: the script is interactive. It will ask the user about several parameters
and fill the blank values of the BlackHawk parameters into new parameter files which names are
user-defined. Then, BH_launcher launches the BlackHawk_inst and/or BlackHawk_tot programs with
the parameters defined by the user as many times as specified.

#########################################################################################

The output of BH_launcher consists in four types of files:

- parameters files are generated following the user's instructions in the present folder

- result files from BlackHawk are stored as usual in the BlackHawk/results/ folder

- nohup files showing the advancement of each BlackHawk run are generated in the main BlackHawk/
folder

- a file containing the number and name of each run is generated in the present folder

#########################################################################################

We advise that you compile BlackHawk with the parameter

#define HARDTABLES

in BlackHawk/src/include.h to save execution time. We point out that this script was designed
for BlackHawk v2.0 and subsequent versions. If you have any problem with that script, please
contact me.