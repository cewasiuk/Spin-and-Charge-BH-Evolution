// This is the source file where the PYTHIA new hadronization
// table is read.
// Last modification: 07 April 2024
// Authors: Alexandre Arbey alexandre.arbey@ens-lyon.fr & Jérémy Auffinger jeremy.auffinger.blackhawk@protonmail.com

#include "include.h"

#ifdef HARDTABLES
#include "./tables/hadronization_tables/hadronization_tables_hdmspectra.h"

void read_hadronization_hdmspectra(double ****tables,double *initial_energies,double *final_energies,struct param *parameters){
	// This function reads the hadronization tables contained in the folder './hdmspectra_tables'

	for(int j = 0;j<parameters->nb_init_en;j++) initial_energies[j]=initial_energies_hdmspectra[j];
	for(int k = 0;k<parameters->nb_fin_en;k++) final_energies[k]=final_energies_hdmspectra[k];
	for(int i = 0;i<parameters->nb_fin_part;i++) for(int j = 0;j<parameters->nb_init_en;j++) for(int k = 0;k<parameters->nb_fin_en;k++) 
	{
		tables[i][j][k][0]=tables_photon_hdmspectra[i][j][k];
		tables[i][j][k][1]=tables_gluon_hdmspectra[i][j][k];
		tables[i][j][k][2]=tables_higgs_hdmspectra[i][j][k];
		tables[i][j][k][3]=tables_W_hdmspectra[i][j][k];
		tables[i][j][k][4]=tables_Z_hdmspectra[i][j][k];
		tables[i][j][k][5]=tables_e_hdmspectra[i][j][k];
		tables[i][j][k][6]=tables_mu_hdmspectra[i][j][k];
		tables[i][j][k][7]=tables_tau_hdmspectra[i][j][k];
		tables[i][j][k][8]=tables_nue_hdmspectra[i][j][k];
		tables[i][j][k][9]=tables_numu_hdmspectra[i][j][k];
		tables[i][j][k][10]=tables_nutau_hdmspectra[i][j][k];
		tables[i][j][k][11]=tables_u_hdmspectra[i][j][k];
		tables[i][j][k][12]=tables_d_hdmspectra[i][j][k];
		tables[i][j][k][13]=tables_c_hdmspectra[i][j][k];
		tables[i][j][k][14]=tables_s_hdmspectra[i][j][k];
		tables[i][j][k][15]=tables_t_hdmspectra[i][j][k];
		tables[i][j][k][16]=tables_b_hdmspectra[i][j][k];
	}

	return;
}

#endif

