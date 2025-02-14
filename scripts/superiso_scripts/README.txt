				#####################
				#		            #
				# BlackHawk scripts #
				#		            #
				#####################



The codes in this folder are designed to reinterpret the dark matter indirect detection
data implemented in SuperIso Relic in the context of primordial black holes.

authors : Alexandre Arbey, alexandre.arbey@ens-lyon.fr
last modified : 08 May 2019

#########################################################################################

To use blackholes.c, you need to download SuperIso Relic (version 4.1 or later), that you
can find here:

	http://superiso.in2p3.fr/relic/

Uncompress the superiso_relic package, and copy blackholes.c into it.

After having run the configure scripts, edit the Makefile of the main directory and set:

	RELIC := 3

Compile SuperIso Relic with shared library and make blackholes.

The usage is:

	./blackholes.x secondary_spectrum mass_BH

where secondary_spectrum is the name of the secondary spectrum file in BlackHawk format, 
and mass_BH the mass of the black hole in grams.

The output of the code will say if the point is valid or excluded by FERMI-LAT or AMS-02
results.

##########################################################################################

If you have any issue using these scripts please feel free to contact the author.























