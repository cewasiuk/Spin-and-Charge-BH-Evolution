// This is the source file where the PYTHIA hadronization
// table is read.
// Last modification: 15 September 2023
// Authors: Jérémy Auffinger jeremy.auffinger.blackhawk@protonmail.com & Alexandre Arbey alexandre.arbey@ens-lyon.fr

#include "include.h"

#ifdef HARDTABLES
#include "./tables/hadronization_tables/hadronization_tables_pythia.h"

void read_hadronization_pythia(double ****tables,double *initial_energies,double *final_energies,struct param *parameters){
	// This function reads the hadronization tables contained in the folder './pythia_tables'

	for(int j = 0;j<parameters->nb_init_en;j++) initial_energies[j]=initial_energies_pythia[j];
	for(int k = 0;k<parameters->nb_fin_en;k++) final_energies[k]=final_energies_pythia[k];
	for(int i = 0;i<parameters->nb_fin_part;i++) for(int j = 0;j<parameters->nb_init_en;j++) for(int k = 0;k<parameters->nb_fin_en;k++) 
	{
		tables[i][j][k][0]=tables_photon_pythia[i][j][k];
		tables[i][j][k][1]=tables_gluon_pythia[i][j][k];
		tables[i][j][k][2]=tables_higgs_pythia[i][j][k];
		tables[i][j][k][3]=tables_W_pythia[i][j][k];
		tables[i][j][k][4]=tables_Z_pythia[i][j][k];
		tables[i][j][k][5]=tables_e_pythia[i][j][k];
		tables[i][j][k][6]=tables_mu_pythia[i][j][k];
		tables[i][j][k][7]=tables_tau_pythia[i][j][k];
		tables[i][j][k][8]=tables_u_pythia[i][j][k];
		tables[i][j][k][9]=tables_d_pythia[i][j][k];
		tables[i][j][k][10]=tables_c_pythia[i][j][k];
		tables[i][j][k][11]=tables_s_pythia[i][j][k];
		tables[i][j][k][12]=tables_t_pythia[i][j][k];
		tables[i][j][k][13]=tables_b_pythia[i][j][k];
	}

	return;
}

#endif

