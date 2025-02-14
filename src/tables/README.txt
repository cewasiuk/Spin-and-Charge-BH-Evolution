				#####################
				#		    		#
				# BlackHawk scripts #
				#		    		#
				#####################



This script is designed to convert the hadronization tables in C hardcoded tables in
order to make BlackHawk execution faster.

author : Alexandre Arbey, alexandre.arbey@ens-lyon.fr and Jérémy Auffinger, jeremy.auffinger.blackhawk@protonmail.com
last modified : 15 September 2023

#########################################################################################

This folder contains 2 files:

- a source file convert_tables.c used to convert the tables.

- a parameters file convert.txt containing the hadronization informations.

#########################################################################################

In order to convert the tables:
1) copy convert_tables.c in the BlackHawk main directory and go there
2) open src/include.h and comment #define HARDTABLES
3) go back to the main directory and type:

	make convert_tables

4) execute it:

	./convert_tables.x

5) remove convert_tables.c and convert_tables.x from the main directory

6) uncomment #define HARDTABLES in src/include.h
7) Type:

	make distclean

8) Type:

	make

#########################################################################################

This script will generate a file hadronization_tables.h in src/tables/hadronization_tables.
This file contains a hardcoded version of the tables, which will be compiled with the other
source files of BlackHawk. The routine which converts the tables is contained in src/secondary.c
and has the prototype:
	
	void convert_hadronization_tables(double ****tables,double *initial_energies,double
		*final_energies,struct param *parameters);
