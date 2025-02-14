// This is the source file where the HAZMA hadronization
// table is read.
// Last modification: 15 September 2023
// Authors: Jérémy Auffinger jeremy.auffinger.blackhawk@protonmail.com & Alexandre Arbey alexandre.arbey@ens-lyon.fr

#include "include.h"

#ifdef HARDTABLES
#include "./tables/hadronization_tables/hadronization_tables_hazma.h"

void read_hadronization_hazma(double ****tables,double *initial_energies,double *final_energies,struct param *parameters){
	// This function reads the hadronization tables contained in the folder './hazma_tables'

	for(int j = 0;j<parameters->nb_init_en;j++) initial_energies[j]=initial_energies_hazma[j];
	for(int k = 0;k<parameters->nb_fin_en;k++) final_energies[k]=final_energies_hazma[k];
	for(int i = 0;i<parameters->nb_fin_part;i++) for(int j = 0;j<parameters->nb_init_en;j++) for(int k = 0;k<parameters->nb_fin_en;k++) 
	{
		tables[i][j][k][0]=tables_photon_hazma[i][j][k];
		tables[i][j][k][1]=tables_electron_hazma[i][j][k];
	}

	return;
}

#endif

