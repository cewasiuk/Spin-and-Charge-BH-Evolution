// This is the source file where the methods computing
// the secondary spectra are implemented.
// Last modification: 07 April 2024
// Authors: Jérémy Auffinger jeremy.auffinger.blackhawk@protonmail.com & Alexandre Arbey alexandre.arbey@ens-lyon.fr

#include "include.h"

void read_hadronization_tables(double ****tables,double *initial_energies,double *final_energies,struct param *parameters){
	// This function reads the hadronization tables contained in the folder '/pythia_tables'
	// (if hadronization_choice = 0), '/herwig_tables' (if hadronization_choice = 1),
	// '/pythia_tables_new' (if hadronization_choice = 2), '/hazma_tables' (if
	// hadronization_choice = 3) or '/hdmspectra_tables' (if
	// hadronization_choice = 4). If HARDTABLES is defined, it directly translates 
	// the header files hadronization_tables_*.h into C tables. It fills the tabulated 
	// initial_energies[], final_energies[], and branching ratios tables[][][][].

#ifdef HARDTABLES

	if(parameters->hadronization_choice==0)
	{
		read_hadronization_pythia(tables,initial_energies,final_energies,parameters);
	}
	else if(parameters->hadronization_choice==1)
	{
		read_hadronization_herwig(tables,initial_energies,final_energies,parameters);
	}
	else if(parameters->hadronization_choice==2)
	{
		read_hadronization_pythianew(tables,initial_energies,final_energies,parameters);
	}
	else if(parameters->hadronization_choice == 3){
		read_hadronization_hazma(tables,initial_energies,final_energies,parameters);
	}
	else if(parameters->hadronization_choice ==4){
		read_hadronization_hdmspectra(tables,initial_energies,final_energies,parameters);
	}
	else
	{
		printf("\n\t [read_hadronization_tables] : ERROR WRONG HADRONIZATION CHOICE !\n");
		fflush(stdout);
		exit(0);
	}

#else
	
	char **table_names = (char **)malloc(parameters->nb_fin_part*sizeof(char *)); // contains the names of the table files
	for(int i = 0;i<parameters->nb_fin_part;i++){
		table_names[i] = (char *)malloc(64*sizeof(char));
	}
	switch(parameters->hadronization_choice){
		case 0:{ // PYTHIA tables
			sprintf(table_names[0],"%s","./src/tables/hadronization_tables/pythia_tables/photon.txt"); // branching ratios to produce photons
			sprintf(table_names[1],"%s","./src/tables/hadronization_tables/pythia_tables/electron.txt"); // branching ratios to produce electrons
			sprintf(table_names[2],"%s","./src/tables/hadronization_tables/pythia_tables/muon.txt"); // branching ratios to produce muons
			sprintf(table_names[3],"%s","./src/tables/hadronization_tables/pythia_tables/nue.txt"); // branching ratios to produce electron neutrinos
			sprintf(table_names[4],"%s","./src/tables/hadronization_tables/pythia_tables/numu.txt"); // branching ratios to produce muon neutrinos
			sprintf(table_names[5],"%s","./src/tables/hadronization_tables/pythia_tables/nutau.txt"); // branching ratios to produce tau neutrinos
			sprintf(table_names[6],"%s","./src/tables/hadronization_tables/pythia_tables/pipm.txt"); // branching ratios to produce pi+- pions
			sprintf(table_names[7],"%s","./src/tables/hadronization_tables/pythia_tables/K0L.txt"); // branching ratios to produce K0L kaons
			sprintf(table_names[8],"%s","./src/tables/hadronization_tables/pythia_tables/Kpm.txt"); // branching ratios to produce K+- kaons
			sprintf(table_names[9],"%s","./src/tables/hadronization_tables/pythia_tables/proton.txt"); // branching ratios to produce protons
			sprintf(table_names[10],"%s","./src/tables/hadronization_tables/pythia_tables/neutron.txt"); // branching ratios to produce neutrons
			break;
		}
		case 1:{ // HERWIG tables
			sprintf(table_names[0],"%s","./src/tables/hadronization_tables/herwig_tables/photon.txt"); // branching ratios to produce photons
			sprintf(table_names[1],"%s","./src/tables/hadronization_tables/herwig_tables/electron.txt"); // branching ratios to produce electrons
			sprintf(table_names[2],"%s","./src/tables/hadronization_tables/herwig_tables/muon.txt"); // branching ratios to produce muons
			sprintf(table_names[3],"%s","./src/tables/hadronization_tables/herwig_tables/nue.txt"); // branching ratios to produce electron neutrinos
			sprintf(table_names[4],"%s","./src/tables/hadronization_tables/herwig_tables/numu.txt"); // branching ratios to produce muon neutrinos
			sprintf(table_names[5],"%s","./src/tables/hadronization_tables/herwig_tables/nutau.txt"); // branching ratios to produce tau neutrinos
			sprintf(table_names[6],"%s","./src/tables/hadronization_tables/herwig_tables/pipm.txt"); // branching ratios to produce pi+- pions
			sprintf(table_names[7],"%s","./src/tables/hadronization_tables/herwig_tables/K0L.txt"); // branching ratios to produce K0L kaons
			sprintf(table_names[8],"%s","./src/tables/hadronization_tables/herwig_tables/Kpm.txt"); // branching ratios to produce K+- kaons
			sprintf(table_names[9],"%s","./src/tables/hadronization_tables/herwig_tables/proton.txt"); // branching ratios to produce protons
			sprintf(table_names[10],"%s","./src/tables/hadronization_tables/herwig_tables/neutron.txt"); // branching ratios to produce neutrons
			break;
		}
		case 2:{ // PYTHIA (today) tables
			sprintf(table_names[0],"%s","./src/tables/hadronization_tables/pythia_tables_new/photon.txt"); // branching ratios to produce photons
			sprintf(table_names[1],"%s","./src/tables/hadronization_tables/pythia_tables_new/electron.txt"); // branching ratios to produce electrons
			sprintf(table_names[2],"%s","./src/tables/hadronization_tables/pythia_tables_new/nue.txt"); // branching ratios to produce electron neutrinos
			sprintf(table_names[3],"%s","./src/tables/hadronization_tables/pythia_tables_new/numu.txt"); // branching ratios to produce muon neutrinos
			sprintf(table_names[4],"%s","./src/tables/hadronization_tables/pythia_tables_new/nutau.txt"); // branching ratios to produce tau neutrinos
			sprintf(table_names[5],"%s","./src/tables/hadronization_tables/pythia_tables_new/proton.txt"); // branching ratios to produce protons
			break;
		}
		case 3:{ // Hazma tables
			sprintf(table_names[0],"%s","./src/tables/hadronization_tables/hazma_tables/photon.txt"); // branching ratios to produce photons
			sprintf(table_names[1],"%s","./src/tables/hadronization_tables/hazma_tables/electron.txt"); // branching ratios to produce electrons
			break;
		}
		case 4:{ // HDMSpectra tables
			sprintf(table_names[0],"%s","./src/tables/hadronization_tables/hdmspectra_tables/photon.txt"); // branching ratios to produce photons
			sprintf(table_names[1],"%s","./src/tables/hadronization_tables/hdmspectra_tables/electron.txt"); // branching ratios to produce electrons
			sprintf(table_names[2],"%s","./src/tables/hadronization_tables/hdmspectra_tables/nue.txt"); // branching ratios to produce electron neutrinos
			sprintf(table_names[3],"%s","./src/tables/hadronization_tables/hdmspectra_tables/numu.txt"); // branching ratios to produce muon neutrinos
			sprintf(table_names[4],"%s","./src/tables/hadronization_tables/hdmspectra_tables/nutau.txt"); // branching ratios to produce tau neutrinos
			sprintf(table_names[5],"%s","./src/tables/hadronization_tables/hdmspectra_tables/proton.txt"); // branching ratios to produce protons
			break;
		}
		default:{
			printf("\n\t [read_hadronization_tables] : ERROR WRONG HADRONIZATION CHOICE !\n");
			fflush(stdout);
			exit(0);
			break;
		}
	}
	FILE *table_file;
	if(parameters->full_output){
		printf("\n\n\t reading:\t");
		fflush(stdout);
	}
	for(int i = 0;i<parameters->nb_fin_part;i++){
		table_file = fopen(table_names[i],"r+");
		if(!table_file){
			printf("\n\t [read_hadronization_tables] : ERROR COULD NOT READ TABLE '%s' !\n",table_names[i]);
			fflush(stdout);
			exit(0);
		}
		rewind(table_file);
		char dummy[64];
		for(int j = 0;j<parameters->nb_init_part+2;j++){ // reading the first line of no interest
			fscanf(table_file,"%s",dummy);
		}
		for(int j = 0;j<parameters->nb_init_en;j++){
			for(int k = 0;k<parameters->nb_fin_en;k++){
				fscanf(table_file,"%lf",&(initial_energies[j]));
				fscanf(table_file,"%lf",&(final_energies[k]));
				for(int l = 0;l<parameters->nb_init_part;l++){
					fscanf(table_file,"%lf",&(tables[i][j][k][l]));
				}
			}
			for(int k = 0;k<parameters->nb_fin_en;k++){
				for(int l = 0;l<parameters->nb_init_part;l++){
					if(k != parameters->nb_fin_en-1){
						if(parameters->hadronization_choice != 3){
							tables[i][j][k][l] = fabs(tables[i][j][k][l])/(final_energies[k+1] - final_energies[k]); // differential branching ratios dN/dE'
						}
						else{
							tables[i][j][k][l] = fabs(tables[i][j][k][l]);
						}
					}
					else{
						tables[i][j][k][l] = fabs(tables[i][j][k][l])/(final_energies[k] - final_energies[k-1]); // differential branching ratios dN/dE'
					}
				}
			}
		}
		if(parameters->full_output){ // prints the state of reading
			if(i == 0){
				printf("%s\n",table_names[i]);
			}
			else{
				printf("\t\t\t%s\n",table_names[i]);
			}
			fflush(stdout);
		}
		fclose(table_file);
	}
	if(parameters->full_output){
		printf("\n");
		fflush(stdout);
	}
	free2D_char(table_names,parameters->nb_fin_part);

#endif

	return;
}

void convert_hadronization_tables(double ****tables,double *initial_energies,double *final_energies,struct param *parameters){
	// This function converts the hadronization tables (.txt files) into C arrays in src/tables/hadro_*.h.
	
	FILE *table_file;
	
	char genname[64],filename[500];
	
	if(parameters->hadronization_choice == 0){
		sprintf(genname,"pythia");
	}
	else if(parameters->hadronization_choice == 1){
		sprintf(genname,"herwig");
	}
	else if(parameters->hadronization_choice == 2){
		sprintf(genname,"pythianew");
	}
	else if(parameters->hadronization_choice == 3){
		sprintf(genname,"hazma");
	}
	else if(parameters->hadronization_choice == 4){
		sprintf(genname,"hdmspectra");
	}
	else
	{
		printf("\n\t [convert_hadronization_tables] : ERROR WRONG HADRONIZATION CHOICE !\n");
		fflush(stdout);
		exit(0);
	}
	sprintf(filename,"./src/tables/hadronization_tables/hadronization_tables_%s.h",genname);
	table_file=fopen(filename,"w+");
	if(!table_file){
		printf("\n\t [convert_hadronization_tables] : ERROR COULD NOT OPEN FILE '%s' !\n",filename);
		fflush(stdout);
		exit(0);
	}
	rewind(table_file);
	
	int i,j,k;
	
	fprintf(table_file,"const double initial_energies_%s[%d]={",genname,parameters->nb_init_en);
	for(j = 0;j<parameters->nb_init_en;j++) fprintf(table_file,"%.3e,",initial_energies[j]);
	fprintf(table_file,"};\n\n");
	
	fprintf(table_file,"const double final_energies_%s[%d]={",genname,parameters->nb_fin_en);
	for(k = 0;k<parameters->nb_fin_en;k++) fprintf(table_file,"%.3e,",final_energies[k]);
	fprintf(table_file,"};\n\n");
		
	fprintf(table_file,"const double tables_photon_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
	for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
	if(tables[i][j][k][0]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][0]);
	fprintf(table_file,"};\n\n");
	
	switch(parameters->hadronization_choice){
		case 0: case 1: case 2:{
			fprintf(table_file,"const double tables_gluon_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][1]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][1]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_higgs_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][2]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][2]);
			fprintf(table_file,"};\n\n");
			
			fprintf(table_file,"const double tables_W_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][3]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][3]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_Z_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][4]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][4]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_e_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][5]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][5]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_mu_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][6]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][6]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_tau_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][7]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][7]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_u_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][8]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][8]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_d_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][9]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][9]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_c_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][10]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][10]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_s_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][11]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][11]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_t_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][12]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][12]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_b_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][13]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][13]);
			fprintf(table_file,"};\n\n");	
			break;
		}
		case 3:{
			fprintf(table_file,"const double tables_electron_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][1]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][13]);
			fprintf(table_file,"};\n\n");
			break;
		}
		case 4:{
			fprintf(table_file,"const double tables_gluon_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][1]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][1]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_higgs_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][2]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][2]);
			fprintf(table_file,"};\n\n");
			
			fprintf(table_file,"const double tables_W_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][3]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][3]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_Z_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][4]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][4]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_e_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][5]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][5]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_mu_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][6]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][6]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_tau_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][7]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][7]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_nue_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][8]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][8]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_numu_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][9]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][9]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_nutau_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][10]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][10]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_u_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][11]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][11]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_d_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][12]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][12]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_c_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][13]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][13]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_s_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][14]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][14]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_t_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][15]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][15]);
			fprintf(table_file,"};\n\n");

			fprintf(table_file,"const double tables_b_%s[%d][%d][%d]={",genname,parameters->nb_fin_part,parameters->nb_init_en,parameters->nb_fin_en);
			for(i = 0;i<parameters->nb_fin_part;i++) for(j = 0;j<parameters->nb_init_en;j++) for(k = 0;k<parameters->nb_fin_en;k++)
			if(tables[i][j][k][16]==0.) fprintf(table_file,"0,"); else fprintf(table_file,"%.3e,",tables[i][j][k][16]);
			fprintf(table_file,"};\n\n");
			break;
		}
	}
	
	fclose(table_file);

	return;
}

void read_writing_instructions(int *write_primary,int *write_secondary,struct param *parameters){
	// This function reads the writing instructions in the files "write_primary.txt"
	// and "write_secondary_*.txt" and fills the writing arrays write_*[] with 0's
	// and 1's.
	
	char file_name[500];
	sprintf(file_name,"./src/tables/write_primary.txt");
	FILE *write = fopen(file_name,"r");
	if(!write){
		printf("\n\t [read_writing_instructions] : ERROR COULD NOT OPEN FILE '%s' !\n",file_name);
		fflush(stdout);
		exit(0);
	}
	rewind(write);
	char dummy[64];
	for(int i = 0;i<parameters->particle_number;i++){
		fscanf(write,"%s",dummy);
		fscanf(write,"%s",dummy);
		fscanf(write,"%i",&(write_primary[i]));
	}
	if(parameters->grav){
		fscanf(write,"%s",dummy);
		fscanf(write,"%s",dummy);
		fscanf(write,"%i",&(write_primary[15]));
		if(parameters->add_DM){
			fscanf(write,"%s",dummy);
			fscanf(write,"%s",dummy);
			fscanf(write,"%i",&(write_primary[16]));
		}
	}
	else if(parameters->add_DM){
		fscanf(write,"%s",dummy);
		fscanf(write,"%s",dummy);
		fscanf(write,"%i",dummy);
		fscanf(write,"%s",dummy);
		fscanf(write,"%s",dummy);
		fscanf(write,"%i",&(write_primary[15]));
	}
	fclose(write);
	
	if(parameters->primary_only == 0){
		if(parameters->hadronization_choice == 0 || parameters->hadronization_choice == 1){
			sprintf(file_name,"./src/tables/write_secondary_BBN.txt");
		}
		else if(parameters->hadronization_choice == 2 || parameters->hadronization_choice == 4){
			sprintf(file_name,"./src/tables/write_secondary_today.txt");
		}
		else if(parameters->hadronization_choice == 3){
			sprintf(file_name,"./src/tables/write_secondary_hazma.txt");
		}
		else if(parameters->hadronization_choice == 4){
			sprintf(file_name,"./src/tables/write_secondary_hdmspectra.txt");
		}
		else{
			printf("\n\t [read_writing_instructions] : ERROR WRONG HADRONIZATION CHOICE !\n");
			fflush(stdout);
			exit(0);
		}
		write = fopen(file_name,"r");
		if(!write){
			printf("\n\t [read_writing_instructions] : ERROR COULD NOT OPEN FILE '%s' !\n",file_name);
			fflush(stdout);
			exit(0);
		}
		rewind(write);
		for(int i = 0;i<parameters->nb_fin_part;i++){
			fscanf(write,"%s",dummy);
			fscanf(write,"%s",dummy);
			fscanf(write,"%i",&(write_secondary[i]));
		}
		fclose(write);
	}
	
	return;
}

void write_lines(char **file_names_primary,char **file_names_secondary,double **partial_primary_spectra,double **partial_integrated_hadronized_spectra,int *write_primary,int *write_secondary,double time,struct param *parameters){
	// This function writes the primary and secondary spectra at a given time in
	// files '*_primary_spectrum.txt' and '*_secondary_spectrum.txt'.
	
	FILE *file;
	for(int i = 0;i<parameters->particle_number+parameters->grav+parameters->add_DM;i++){
		if(write_primary[i]){
			file = fopen(file_names_primary[i],"a");
			if(!file){
				printf("\n\t [write_lines] : ERROR COULD NOT OPEN FILE '%s' !\n",file_names_primary[i]);
				fflush(stdout);
				exit(0);
			}
			fprintf(file,"%15.5e",time/time_conversion);
			for(int j = 0;j<parameters->E_number;j++){
				fprintf(file,"%15.5e",partial_primary_spectra[i][j]*rate_conversion); // conversion from GeV to CGS units, plus the 1.e-100 correction
			}
			fprintf(file,"\n");
			fclose(file);
		}
	}	
	if(parameters->primary_only == 0){
		for(int i = 0;i<parameters->nb_fin_part;i++){
			if(write_secondary[i]){
				file = fopen(file_names_secondary[i],"a");
				if(!file){
					printf("\n\t [write_lines] : ERROR COULD NOT OPEN FILE '%s' !\n",file_names_secondary[i]);
					fflush(stdout);
					exit(0);
				}
				fprintf(file,"%15.5e",time/time_conversion);
				for(int j = 0;j<parameters->nb_fin_en;j++){
					fprintf(file,"%15.5e",partial_integrated_hadronized_spectra[i][j]*rate_conversion); // conversion from GeV to CGS units, plus the 1.e-100 correction
				}
				fprintf(file,"\n");
				fclose(file);
			}
		}
	}
	
	return;
}

void total_spectra(double ***partial_hadronized_spectra,double **partial_primary_spectra,double **partial_integrated_hadronized_spectra,
	double ****tables,double *initial_energies,double *final_energies,double **spec_table,double *times,double ***life_masses,double **BH_masses,
	double ***life_params,double **BH_params,double *energies,double *masses_primary,double *spins,double *dof,double *masses_secondary,
	double ***gammas,double *gamma_param,double *gamma_x,double ***fits,struct param *parameters){
	// This function computes the primary and secondary spectra of all particles at all times,
	// by calling the "instantaneous" functions. It fills the partial arrays partial_primary_spectra[][],
	// partial_hadronized_spectra[][][] and partial_integrated_hadronized_spectra[][].
	
	char **file_names_primary = (char **)malloc((parameters->particle_number+parameters->grav+parameters->add_DM)*sizeof(char *)); // contains the names of the output files
	for(int i = 0;i<parameters->particle_number+parameters->grav+parameters->add_DM;i++){
		file_names_primary[i] = (char *)malloc(64*sizeof(char));
	}
	sprintf(file_names_primary[0],"./results/%s/photon_primary_spectrum.txt",parameters->destination_folder);
	sprintf(file_names_primary[1],"./results/%s/gluon_primary_spectrum.txt",parameters->destination_folder);
	sprintf(file_names_primary[2],"./results/%s/higgs_primary_spectrum.txt",parameters->destination_folder);
	sprintf(file_names_primary[3],"./results/%s/wpm_primary_spectrum.txt",parameters->destination_folder);
	sprintf(file_names_primary[4],"./results/%s/z0_primary_spectrum.txt",parameters->destination_folder);
	sprintf(file_names_primary[5],"./results/%s/neutrinos_primary_spectrum.txt",parameters->destination_folder);
	sprintf(file_names_primary[6],"./results/%s/electron_primary_spectrum.txt",parameters->destination_folder);
	sprintf(file_names_primary[7],"./results/%s/muon_primary_spectrum.txt",parameters->destination_folder);
	sprintf(file_names_primary[8],"./results/%s/tau_primary_spectrum.txt",parameters->destination_folder);
	sprintf(file_names_primary[9],"./results/%s/up_primary_spectrum.txt",parameters->destination_folder);
	sprintf(file_names_primary[10],"./results/%s/down_primary_spectrum.txt",parameters->destination_folder);
	sprintf(file_names_primary[11],"./results/%s/charm_primary_spectrum.txt",parameters->destination_folder);
	sprintf(file_names_primary[12],"./results/%s/strange_primary_spectrum.txt",parameters->destination_folder);
	sprintf(file_names_primary[13],"./results/%s/top_primary_spectrum.txt",parameters->destination_folder);
	sprintf(file_names_primary[14],"./results/%s/bottom_primary_spectrum.txt",parameters->destination_folder);
	if(parameters->grav == 1){
		sprintf(file_names_primary[15],"./results/%s/graviton_primary_spectrum.txt",parameters->destination_folder);
		if(parameters->add_DM){
			if(parameters->hadronization_choice != 3){
				sprintf(file_names_primary[16],"./results/%s/DM_primary_spectrum.txt",parameters->destination_folder);
			}
			else{
				sprintf(file_names_primary[16],"./results/%s/pion_primary_spectrum.txt",parameters->destination_folder);				
			}
		}
	}
	else{
		if(parameters->add_DM){
			sprintf(file_names_primary[15],"./results/%s/DM_primary_spectrum.txt",parameters->destination_folder);
		}
	}
	FILE *file;
	for(int i = 0;i<parameters->particle_number+parameters->grav+parameters->add_DM;i++){
		file = fopen(file_names_primary[i],"w+");
		if(!file){
			printf("\n\t [total_spectra] : ERROR COULD NOT OPEN FILE '%s' !\n",file_names_primary[i]);
			fflush(stdout);
			exit(0);
		};
		rewind(file);
		fprintf(file,"Hawking primary spectrum as a function of time.\n");
		fprintf(file,"%15s","time/energy");
		for(int l = 0;l<parameters->E_number;l++){
			fprintf(file,"%15.5e",energies[l]);
		}
		fprintf(file,"\n");
		fclose(file);
	}
	
	char **file_names_secondary; // contains the output file names
	if(parameters->primary_only == 0){
		file_names_secondary = (char **)malloc(parameters->nb_fin_part*sizeof(char *));
		for(int i = 0;i<parameters->nb_fin_part;i++){
			file_names_secondary[i] = (char *)malloc(64*sizeof(char));
		}
		switch(parameters->hadronization_choice){
			case 0: case 1:{ // PYTHIA or HERWIG tables
				sprintf(file_names_secondary[0],"./results/%s/photon_secondary_spectrum.txt",parameters->destination_folder);
				sprintf(file_names_secondary[1],"./results/%s/electron_secondary_spectrum.txt",parameters->destination_folder);
				sprintf(file_names_secondary[2],"./results/%s/muon_secondary_spectrum.txt",parameters->destination_folder);
				sprintf(file_names_secondary[3],"./results/%s/neutrino_e_secondary_spectrum.txt",parameters->destination_folder);
				sprintf(file_names_secondary[4],"./results/%s/neutrino_mu_secondary_spectrum.txt",parameters->destination_folder);
				sprintf(file_names_secondary[5],"./results/%s/neutrino_tau_secondary_spectrum.txt",parameters->destination_folder);
				sprintf(file_names_secondary[6],"./results/%s/pipm_secondary_spectrum.txt",parameters->destination_folder);
				sprintf(file_names_secondary[7],"./results/%s/K0L_secondary_spectrum.txt",parameters->destination_folder);
				sprintf(file_names_secondary[8],"./results/%s/Kpm_secondary_spectrum.txt",parameters->destination_folder);
				sprintf(file_names_secondary[9],"./results/%s/proton_secondary_spectrum.txt",parameters->destination_folder);
				sprintf(file_names_secondary[10],"./results/%s/neutron_secondary_spectrum.txt",parameters->destination_folder);
				break;
			}
			case 2: case 4:{ // PYTHIA (today) tables or HDMSpectra tables
				sprintf(file_names_secondary[0],"./results/%s/photon_secondary_spectrum.txt",parameters->destination_folder);
				sprintf(file_names_secondary[1],"./results/%s/electron_secondary_spectrum.txt",parameters->destination_folder);
				sprintf(file_names_secondary[2],"./results/%s/nu_e_secondary_spectrum.txt",parameters->destination_folder);
				sprintf(file_names_secondary[3],"./results/%s/nu_mu_secondary_spectrum.txt",parameters->destination_folder);
				sprintf(file_names_secondary[4],"./results/%s/nu_tau_secondary_spectrum.txt",parameters->destination_folder);
				sprintf(file_names_secondary[5],"./results/%s/proton_secondary_spectrum.txt",parameters->destination_folder);
				break;
			}
			case 3:{ // Hazma tables
				sprintf(file_names_secondary[0],"./results/%s/photon_secondary_spectrum.txt",parameters->destination_folder);
				sprintf(file_names_secondary[1],"./results/%s/electron_secondary_spectrum.txt",parameters->destination_folder);
				break;
			}
			default:{
				printf("\n\t [total_spectra] : ERROR WRONG HADRONIZATION CHOICE !\n");
				fflush(stdout);
				exit(0);
				break;
			}
		}
		for(int i = 0;i<parameters->nb_fin_part;i++){
			file = fopen(file_names_secondary[i],"w+");
			if(!file){
				printf("\n\t [total_spectra] : ERROR COULD NOT OPEN FILE '%s' !\n",file_names_secondary[i]);
				fflush(stdout);
				exit(0);
			}
			rewind(file);
			fprintf(file,"Hawking secondary spectrum as a function of time.\n");
			fprintf(file,"%15s","time/energy");
			for(int j = 0;j<parameters->nb_fin_en;j++){
				fprintf(file,"%15.5e",final_energies[j]);
			}
			fprintf(file,"\n");
			fclose(file);
		}
	}
	
	int full_output = parameters->full_output; // We don't print the "full output" of the next functions but save the parameter for later
	parameters->full_output = 0;
	
	int *write_primary = (int *)malloc((parameters->particle_number+parameters->grav+parameters->add_DM)*sizeof(int));
	int *write_secondary = (int *)malloc(parameters->nb_fin_part*sizeof(int));
	
	read_writing_instructions(write_primary,write_secondary,parameters);
	
	int *compute_primary = (int *)malloc((parameters->particle_number+parameters->grav+parameters->add_DM)*sizeof(int)); // this table decides whether the corresponding spectra are computed or not
	compute_primary[0] = 1; // photons
	compute_primary[1] = 1; // gluons
	compute_primary[2] = 1; // higgs boson
	compute_primary[3] = 1; // W+- boson
	compute_primary[4] = 1; // Z0 boson
	compute_primary[5] = 1; // neutrinos
	compute_primary[6] = 1; // electron
	compute_primary[7] = 1; // muon
	compute_primary[8] = 1; // tau
	compute_primary[9] = 1; // up quark
	compute_primary[10] = 1; // down quark
	compute_primary[11] = 1; // charm quark
	compute_primary[12] = 1; // strange quark
	compute_primary[13] = 1; // top quark
	compute_primary[14] = 1; // bottom quark
	if(parameters->grav){
		compute_primary[15] = 1; // graviton
		if(parameters->add_DM){
			compute_primary[16] = 1;
		}
	}
	else if(parameters->add_DM){
		compute_primary[15] = 1;
	}
	int *compute_secondary = (int *)malloc(parameters->nb_fin_part*sizeof(int)); // decides whether we compute the contribution from this particle or not
	switch(parameters->hadronization_choice){
		case 0: case 1:{ // BBN epoch particles
			compute_secondary[0] = 1; // photon
			compute_secondary[1] = 1; // electron
			compute_secondary[2] = 1; // muon
			compute_secondary[3] = 1; // nu_e
			compute_secondary[4] = 1; // nu_mu
			compute_secondary[5] = 1; // nu_tau
			compute_secondary[6] = 1; // pi+-
			compute_secondary[7] = 1; // K0 long
			compute_secondary[8] = 1; // K+-
			compute_secondary[9] = 1; // proton
			compute_secondary[10] = 1; // neutron	
			break;
		}
		case 2: case 4:{ // present epoch particles
			compute_secondary[0] = 1; // photon
			compute_secondary[1] = 1; // electron
			compute_secondary[2] = 1; // nu_e
			compute_secondary[3] = 1; // nu_mu
			compute_secondary[4] = 1; // nu_tau
			compute_secondary[5] = 1; // proton
			break;
		}
		case 3:{ // Hazma particles
			compute_secondary[0] = 1; // photon
			compute_secondary[1] = 1; // electron
			break;
		}
		default:{
			printf("\n\t [total_spectra] : ERROR WRONG HADRONIZATION CHOICE !\n");
			fflush(stdout);
			exit(0);
			break;
		}
	}
	
	int **counters_a = (int **)malloc(parameters->BH_number*sizeof(int *)); // contains the running position in the tabulated spins
	int ***counters_x = (int ***)malloc(parameters->BH_number*sizeof(int **)); // contains the running position in the tabulated energies
	for(int i = 0;i<parameters->BH_number;i++){
		counters_a[i] = (int *)malloc(parameters->param_number*sizeof(int));
		counters_x[i] = (int **)malloc(parameters->param_number*sizeof(int *));
		for(int j = 0;j<parameters->param_number;j++){
			counters_x[i][j] = (int *)malloc(parameters->E_number*sizeof(int));
		}
	}
	
	if(full_output){
		printf("\n\n\t computing :\t");
		fflush(stdout);
	}
	for(int i = 0;i<parameters->nb_fin_times;i++){ // for each time we do the instantaneous computation
		
		for(int j = 0;j<parameters->BH_number;j++){
			for(int k = 0;k<parameters->param_number;k++){
				BH_masses[j][k] = life_masses[j][k][i];
				BH_params[j][k] = life_params[j][k][i];
			}
		}
		
		instantaneous_primary_spectrum(partial_primary_spectra,BH_masses,BH_params,spec_table,energies,gammas,gamma_param,gamma_x,fits,dof,spins,masses_primary,counters_a,counters_x,compute_primary,parameters);
		
		if(parameters->primary_only == 0){
			hadronize_instantaneous(partial_hadronized_spectra,tables,initial_energies,final_energies,partial_primary_spectra,energies,masses_secondary,compute_secondary,parameters);
			
			integrate_initial_energies_instantaneous(partial_hadronized_spectra,partial_integrated_hadronized_spectra,energies,final_energies,compute_secondary,parameters);
			
			if(compute_secondary[0]){
				add_photons_instantaneous(partial_primary_spectra,partial_integrated_hadronized_spectra,energies,final_energies,parameters);
			}
			if(((parameters->hadronization_choice == 0 || parameters->hadronization_choice == 1) && compute_secondary[3]+compute_secondary[4]+compute_secondary[5]) 
			|| ((parameters->hadronization_choice == 2 || parameters->hadronization_choice == 4) && compute_secondary[2]+compute_secondary[3]+compute_secondary[4])){
				add_neutrinos_instantaneous(partial_primary_spectra,partial_integrated_hadronized_spectra,energies,final_energies,parameters);
			}
			if(compute_secondary[1]){
				add_electrons_instantaneous(partial_primary_spectra,partial_integrated_hadronized_spectra,energies,final_energies,parameters);
			}
		}
		
		write_lines(file_names_primary,file_names_secondary,partial_primary_spectra,partial_integrated_hadronized_spectra,write_primary,write_secondary,times[i],parameters);
		
		if(full_output == 1 && (i+1)%(parameters->nb_fin_times/6) == 0){
			printf("%i/%i ",i+1,parameters->nb_fin_times);
			fflush(stdout);
		}
		else if(full_output == 1 && i == parameters->nb_fin_times-1){
			printf("%i/%i",i+1,parameters->nb_fin_times);
			fflush(stdout);
		}
	}
	if(full_output){
		printf("\n\n");
		fflush(stdout);
	}
	
	parameters->full_output = full_output;
	free2D_int(counters_a,parameters->BH_number);
	free3D_int(counters_x,parameters->BH_number,parameters->param_number);
	free1D_int(compute_primary);
	free1D_int(compute_secondary);
	free1D_int(write_primary);
	free1D_int(write_secondary);
	free2D_char(file_names_primary,parameters->particle_number+parameters->grav+parameters->add_DM);
	if(parameters->primary_only == 0){
		free2D_char(file_names_secondary,parameters->nb_fin_part);
	}
	return;
}

double contribution_instantaneous(int j,int counter,int k,double **instantaneous_primary_spectra,double ****tables,double *initial_energies,
	double *final_energies,int particle_type,struct param *parameters){
	// This function computes the contribution to the secondary spectra of particle particle_type
	// at initial energy j, closest tabulated initial energy counter and final energy k by
	// summing over the channels of production (i) the initial spectra times the branching ratios.
	
	double temp = 0.;
	switch(parameters->hadronization_choice){
		case 0: case 1: case 2:{ // PYTHIA and HERWIG tables
			for(int i = 0;i<parameters->nb_init_part;i++){ // photons count as 0 anyway; neutrinos are skipped
				temp += instantaneous_primary_spectra[i+(i>4)][j]*tables[particle_type][counter][k][i];
			}
			break;
		}
		case 3:{ // Hazma tables, we must detail the contributions
			temp += instantaneous_primary_spectra[15+parameters->grav][j]*tables[particle_type][counter][k][0]; // contribution from the pi0
			temp += 2.*instantaneous_primary_spectra[15+parameters->grav][j]*tables[particle_type][counter][k][1]; // contribution from the pi+- (we computed the pion spectrum for a single dof)
			temp += instantaneous_primary_spectra[7][j]*tables[particle_type][counter][k][2]; // contribution from the muons
			break;
		}
		case 4:{ // HDMSpectra tables
			for(int i = 0;i<14;i++){ // neutrinos are detailed below
				temp += instantaneous_primary_spectra[i+(i>4)][j]*tables[particle_type][counter][k][i+3*(i>7)]; // neutrinos treated separately
			}
			temp += instantaneous_primary_spectra[5][j]*tables[particle_type][counter][k][8]; // nu_e contribution
			temp += instantaneous_primary_spectra[5][j]*tables[particle_type][counter][k][9]; // nu_mu contribution
			temp += instantaneous_primary_spectra[5][j]*tables[particle_type][counter][k][10]; // nu_tau contribution
			break;
		}
		default:{
			printf("\n\t [contribution_instantaneous] : ERROR WRONG HADRONIZATION CHOICE !\n");
			fflush(stdout);
			exit(0);
		}
	}
	return temp;
}

void hadronize_instantaneous(double ***instantaneous_hadronized_spectra,double ****tables,double *initial_energies,double *final_energies,
	double **instantaneous_primary_spectra,double *energies,double *masses_secondary,int *compute,struct param *parameters){
	// This function computes the Hawking instantaneous secondary spectra by calling the function
	// contribution_instantaneous and fills the array instantaneous_hadronized_spectra[][][].
	// The initial energies are still separated.
	
	int lim_init_en; // lowest energy with trustworthy hadronization results
	int lim_fin_en; // start of the results at this energy
	switch(parameters->hadronization_choice){
		case 0: case 2:{ // PYTHIA tables, PYTHIA new tables
			lim_init_en = 0;
			lim_fin_en = 304;
			break;
		}
		case 1:{ // HERWIG tables
			lim_init_en = 0;
			lim_fin_en = 64;
			break;
		}
		case 3: case 4:{
			// extrapolation is not available with the Hazma and HDMSpectra tables
			lim_init_en = 0;
			lim_fin_en = 0;
			break;
		}
		default:{
			printf("\n\t [hadronize_instantaneous] : ERROR WRONG HADRONIZATION CHOICE !\n");
			fflush(stdout);
			exit(0);
			break;
		}
	}
	int counter;
	int counter2;
	double contrib1;
	double contrib2;
	if(parameters->full_output){
		printf("\n\n\t hadronize:\t");
		fflush(stdout);
	}
	for(int i = 0;i<parameters->nb_fin_part;i++){
		if(compute[i]){
			counter = 0;
			for(int j = 0;j<parameters->E_number;j++){
				if(parameters->hadronization_choice == 3 || parameters->hadronization_choice == 4 || (energies[j] >= initial_energies[lim_init_en] && energies[j] <= initial_energies[parameters->nb_init_en-1])){ // we use the computed tables
					if(fabs(energies[j] - initial_energies[counter]) < machine_precision){ // if the initial energy is in the tables
						for(int k = 0;k<parameters->nb_fin_en;k++){
							instantaneous_hadronized_spectra[i][j][k] = contribution_instantaneous(j,counter,k,instantaneous_primary_spectra,tables,initial_energies,final_energies,i,parameters);
						}
						counter++;
					}
					else{ // else we have to interpolate the value
						if(parameters->interpolation_method > 1 || parameters->interpolation_method < 0)
						{
							printf("\n\t [hadronize_instantaneous] : ERROR WRONG INTERPOLATION METHOD !\n");
							fflush(stdout);
							exit(0);
						}
						
						if(energies[j] < initial_energies[counter]){

							for(int k = 0;k<parameters->nb_fin_en;k++){
								contrib1 = contribution_instantaneous(j,counter-1,k,instantaneous_primary_spectra,tables,initial_energies,final_energies,i,parameters);
								contrib2 = contribution_instantaneous(j,counter,k,instantaneous_primary_spectra,tables,initial_energies,final_energies,i,parameters);
								if(parameters->interpolation_method == 0){ // linear interpolation
									instantaneous_hadronized_spectra[i][j][k] = contrib1 + (contrib2 - contrib1)/(initial_energies[counter] - initial_energies[counter-1])*(energies[j] - initial_energies[counter-1]);
								}
								else if(parameters->interpolation_method == 1){
									if(contrib1 <= 0. || contrib2 <= 0.){
										instantaneous_hadronized_spectra[i][j][k] = 0.;
									}
									else{
										instantaneous_hadronized_spectra[i][j][k] = pow(10.,log10(contrib1) + (log10(contrib2) - log10(contrib1))/(log10(initial_energies[counter]) - log10(initial_energies[counter-1]))*(log10(energies[j]) - log10(initial_energies[counter-1])));
									}
								}
							}
						}
						else{
							counter++;
							j--;
						}
					}
				}
				else if(parameters->hadronization_choice != 3 && parameters->hadronization_choice != 4 && energies[j] < initial_energies[lim_init_en]){ // we extrapolate the result of the computed tables by shifting it at lower energies
					counter2 = lim_fin_en - 1;//parameters->nb_fin_en-1;
					while(final_energies[counter2] > energies[j] && counter2 > 0){ // finds the shift
						counter2--;
					}
					for(int k = 0;k<counter2+1;k++){ // we copy the results of the table at a more energetic line
						if(final_energies[k] > masses_secondary[i] && final_energies[k] < energies[j]){ // if the energy is sufficient
							instantaneous_hadronized_spectra[i][j][k] = contribution_instantaneous(j,lim_init_en,lim_fin_en-counter2+k,instantaneous_primary_spectra,tables,initial_energies,final_energies,i,parameters);
						}
					}
				}
				else if(parameters->hadronization_choice == 3){
					printf("\n [hadronize_instantaneous] : ERROR EXTRAPOLATION IS NOT AVAILABLE WITH THE HAZMA TABLES !\n");
					fflush(stdout);
					exit(0);
				}
				else if(parameters->hadronization_choice == 4){
					printf("\n [hadronize_instantaneous] : ERROR EXTRAPOLATION IS NOT AVAILABLE WITH THE HDMSpectra TABLES !\n");
					fflush(stdout);
					exit(0);
				}
			}
		}
		if(parameters->full_output){ // prints the current state of hadronization
			printf("%i/%i\t",i+1,parameters->nb_fin_part);
			fflush(stdout);
		}
	}
	if(parameters->full_output){
		printf("\n\n");
		fflush(stdout);
	}
	return;
}

void integrate_initial_energies_instantaneous(double ***instantaneous_hadronized_spectra,double **instantaneous_integrated_hadronized_spectra,
	double *energies,double *final_energies,int *compute,struct param *parameters){
	// This function integrates the initial energy dependency to give the secondary Hawking
	// spectra, contained in the array instantaneous_integrated_hadronized_spectra[][].
	
	if(parameters->full_output){
		printf("\n\n\t integrate:\t");
		fflush(stdout);
	}
	for(int i = 0;i<parameters->nb_fin_part;i++){
		if(compute[i]){
/*#if defined(_OPENMP)
#pragma omp parallel for
#endif*/
			for(int k = 0;k<parameters->nb_fin_en;k++){
				instantaneous_integrated_hadronized_spectra[i][k] = 0.;
				for(int l = 0;l<parameters->E_number-1;l++){
					instantaneous_integrated_hadronized_spectra[i][k] += trapeze(energies[l],energies[l+1],instantaneous_hadronized_spectra[i][l][k],instantaneous_hadronized_spectra[i][l+1][k]);
				}
			}
			if(parameters->full_output){ // prints the current state of integration
				printf("%i/%i\t",i+1,parameters->nb_fin_part);
				fflush(stdout);
			}
		}
	}
	if(parameters->full_output){
		printf("\n\n");
		fflush(stdout);
	}
	return;
}

void add_photons_instantaneous(double **instantaneous_primary_spectra,double **instantaneous_integrated_hadronized_spectra,
	double *energies,double *final_energies,struct param *parameters){
	// This function adds the contribution from primary photons to the secondary spectrum.
	
	if(parameters->full_output){
		printf("\t adding primary photons...");
		fflush(stdout);
	}
	int counter = 0; // position in initial energies
	for(int k = 0;k<parameters->nb_fin_en;k++){
		while(counter < parameters->E_number && energies[counter] < final_energies[k]){
			counter++;
		}
		if(counter>0 && counter < parameters->E_number){
			if(parameters->interpolation_method == 0){ // linear interpolation
				instantaneous_integrated_hadronized_spectra[0][k] += fabs(instantaneous_primary_spectra[0][counter-1] + (instantaneous_primary_spectra[0][counter] - instantaneous_primary_spectra[0][counter-1])/(energies[counter] - energies[counter-1])*(final_energies[k] - energies[counter-1]));
			}
			else if(parameters->interpolation_method == 1){ // logarithmic interpolation
				if(instantaneous_primary_spectra[0][counter-1] <= 0. || instantaneous_primary_spectra[0][counter] <= 0.){
					instantaneous_integrated_hadronized_spectra[0][k] += 0.;
				}
				else{
					instantaneous_integrated_hadronized_spectra[0][k] += fabs(pow(10.,log10(instantaneous_primary_spectra[0][counter-1]) + (log10(instantaneous_primary_spectra[0][counter]) - log10(instantaneous_primary_spectra[0][counter-1]))/(log10(energies[counter]) - log10(energies[counter-1]))*(log10(final_energies[k]) - log10(energies[counter-1]))));
				}
			}
			else{
				printf("\n\t [add_photons_instantaneous] : ERROR WRONG INTERPOLATION METHOD !\n");
				fflush(stdout);
				exit(0);
			}
		}
	}
	if(parameters->full_output){
		printf("\t DONE");
		fflush(stdout);
	}
	return;
}

void add_neutrinos_instantaneous(double **instantaneous_primary_spectra,double **instantaneous_integrated_hadronized_spectra,
	double *energies,double *final_energies,struct param *parameters){
	// This function adds the contribution from primary neutrinos to the secondary spectra.
	
	if(parameters->full_output){
		printf("\n\n\t adding primary neutrinos...");
		fflush(stdout);
	}
	int counter = 0; // position in initial energies
	double temp;
	for(int k = 0;k<parameters->nb_fin_en;k++){
		while(counter < parameters->E_number && energies[counter] < final_energies[k]){
			counter++;
		}
		if(counter>0 && counter < parameters->E_number){
			if(parameters->interpolation_method == 0){ // linear interpolation
				temp = fabs(instantaneous_primary_spectra[5][counter-1] + (instantaneous_primary_spectra[5][counter] - instantaneous_primary_spectra[5][counter-1])/(energies[counter] - energies[counter-1])*(final_energies[k] - energies[counter-1]));
			}
			else if(parameters->interpolation_method == 1){ // logarithmic interpolation
				if(instantaneous_primary_spectra[5][counter-1] == 0. || instantaneous_primary_spectra[5][counter] == 0.){
					temp = 0.;
				}
				else{
					temp = fabs(pow(10.,log10(instantaneous_primary_spectra[5][counter-1]) + (log10(instantaneous_primary_spectra[5][counter]) - log10(instantaneous_primary_spectra[5][counter-1]))/(log10(energies[counter]) - log10(energies[counter-1]))*(log10(final_energies[k]) - log10(energies[counter-1]))));
				}
			}
			else{
				printf("\n\t [add_neutrinos_instantaneous] : ERROR WRONG INTERPOLATION METHOD !\n");
				fflush(stdout);
				exit(0);
			}
			switch(parameters->hadronization_choice){
				case 0: case 1:{ // PYTHIA or HERWIG tables
					instantaneous_integrated_hadronized_spectra[3][k] += temp/3.;
					instantaneous_integrated_hadronized_spectra[4][k] += temp/3.;
					instantaneous_integrated_hadronized_spectra[5][k] += temp/3.;
					break;
				}
				case 2: case 4:{ // PYTHIA today or HDMSpectra tables
					instantaneous_integrated_hadronized_spectra[2][k] += temp/3.;
					instantaneous_integrated_hadronized_spectra[3][k] += temp/3.;
					instantaneous_integrated_hadronized_spectra[4][k] += temp/3.;
					break;
				}
				default:{
					printf("\n\t [add_neutrinos_instantaneous] : ERROR WRONG HADRONIZATION CHOICE !\n");
					fflush(stdout);
					exit(0);
					break;
				}
			}
		}
	}
	if(parameters->full_output){
		printf("\t DONE");
		fflush(stdout);
	}
	return;
}

void add_electrons_instantaneous(double **instantaneous_primary_spectra,double **instantaneous_integrated_hadronized_spectra,
	double *energies,double *final_energies,struct param *parameters){
	// This function adds the contribution from primary electrons to the secondary spectra.
	
	if(parameters->full_output){
		printf("\n\n\t adding primary electrons...");
		fflush(stdout);
	}
	int counter = 0; // position in initial energies
	double temp;
	for(int k = 0;k<parameters->nb_fin_en;k++){
		while(counter < parameters->E_number && energies[counter] < final_energies[k]){
			counter++;
		}
		if(counter>0 && counter < parameters->E_number){
			if(parameters->interpolation_method == 0){ // linear interpolation
				temp = fabs(instantaneous_primary_spectra[6][counter-1] + (instantaneous_primary_spectra[6][counter] - instantaneous_primary_spectra[6][counter-1])/(energies[counter] - energies[counter-1])*(final_energies[k] - energies[counter-1]));
			}
			else if(parameters->interpolation_method == 1){ // logarithmic interpolation
				if(instantaneous_primary_spectra[6][counter-1] == 0. || instantaneous_primary_spectra[6][counter] == 0.){
					temp = 0.;
				}
				else{
					temp = fabs(pow(10.,log10(instantaneous_primary_spectra[6][counter-1]) + (log10(instantaneous_primary_spectra[6][counter]) - log10(instantaneous_primary_spectra[6][counter-1]))/(log10(energies[counter]) - log10(energies[counter-1]))*(log10(final_energies[k]) - log10(energies[counter-1]))));
				}
			}
			else{
				printf("\n\t [add_electrons_instantaneous] : ERROR WRONG INTERPOLATION METHOD !\n");
				fflush(stdout);
				exit(0);
			}
			if(final_energies[k] >= m_electron){ // just in case we have an interpolation problem
				instantaneous_integrated_hadronized_spectra[1][k] += temp;
			}
		}
	}
	if(parameters->full_output){
		printf("\t DONE\n\n");
		fflush(stdout);
	}
	return;
}

void add_FSR_instantaneous(double **instantaneous_primary_spectra,double **instantaneous_integrated_hadronized_spectra,
	double *energies,double *final_energies,double *masses_primary,struct param *parameters){
	// This function computes the final state radiation of electrons, muons and charged pions
	// if the Hazma hadronization at low energy is chosen. See https://arxiv.org/pdf/1907.11846.pdf
	// for the formulas.
	
	double Q;
	double x;
	double mu_epm;
	double mu_mupm;
	double mu_pipm;
	for(int i = 0;i<parameters->nb_fin_en;i++){
		for(int j = 0;j<parameters->E_number-1;j++){
			Q = 2.*energies[j];
			x = 2.*final_energies[i]/Q;
			mu_epm = masses_primary[6]/Q;
			mu_mupm = masses_primary[7]/Q;
			mu_pipm = m_pipm/Q;
			if(energies[j] > final_energies[i] && log((1.-x)/pow(mu_epm,2.))-1. > 0.){ // contribution from e+-
				instantaneous_integrated_hadronized_spectra[0][i] += alpha_em/(pi*Q)*(1.+pow(1.-x,2.))/x*(log((1.-x)/pow(mu_epm,2.))-1.)*(instantaneous_primary_spectra[6][j])*(energies[j+1]-energies[j]);
			}
			if(energies[j] > final_energies[i] && log((1.-x)/pow(mu_mupm,2.))-1. > 0.){ // contribution from mu+-
				instantaneous_integrated_hadronized_spectra[0][i] += alpha_em/(pi*Q)*(1.+pow(1.-x,2.))/x*(log((1.-x)/pow(mu_mupm,2.))-1.)*(instantaneous_primary_spectra[7][j])*(energies[j+1]-energies[j]);
			}
			if(energies[j] > m_pipm && energies[j] > final_energies[i] && (1.-x)*log((1.-x)/pow(mu_mupm,2.))-1. > 0.){ // contribution from pi+-, we put back manually the condition E > m_pipm as the spectrum of pions is computed for m_pi0 (factor 2 for anti-particles already taken into account by the formula, but we computed the pion spectrum for a single dof)
				instantaneous_integrated_hadronized_spectra[0][i] += alpha_em/(pi*Q)*2.*(1.-x)/x*(log((1.-x)/pow(mu_pipm,2.))-1.)*(instantaneous_primary_spectra[15+parameters->grav][j]*2.)*(energies[j+1]-energies[j]);
			}
		}
	}
	return;
}

void write_instantaneous_secondary_spectra(double **instantaneous_integrated_hadronized_spectra,double *hadronized_energies,struct param *parameters){
	// This function writes the Hawking instantaneous secondary spectra of a distribution of BH,
	// contained in the array instantaneous_integrated_hadronized_spectra[][], into the file 
	// 'instantaneous_secondary_spectra.txt'.
	
	char file_name[500];
	sprintf(file_name,"./results/%s/instantaneous_secondary_spectra.txt",parameters->destination_folder);
	FILE *file = fopen(file_name,"w+");
	if(!file){
		printf("\n\t [write_instantaneous_secondary_spectra] : ERROR COULD NOT OPEN FILE '%s' !\n",file_name);
		fflush(stdout);
		exit(0);
	}
	rewind(file);
	fprintf(file,"Hawking secondary spectra for each particle type.\n");
	switch(parameters->hadronization_choice){
		case 0: case 1:{ // PYTHIA or HERWIG hadronization tables
			fprintf(file,"%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n","energy/particle","photon","electron","muon","nu_e","nu_mu","nu_tau","pipm","K0L","Kpm","proton","neutron");
			break;
		}
		case 2: case 4:{ // PYTHIA today or HDMSpectra hadronization tables
			fprintf(file,"%15s%15s%15s%15s%15s%15s%15s\n","energy/particle","photon","electron","nu_e","nu_mu","nu_tau","proton");
			break;
		}
		case 3:{ // Hazma hadronization tables
			fprintf(file,"%15s%15s%15s\n","energy/particle","photon","electron");
			break;
		}
		default:{
			printf("\n\t [write_instantaneous_secondary_spectra] : ERROR WRONG HADRONIZATION CHOICE !\n");
			fflush(stdout);
			exit(0);
			break;
		}
	}
	for(int i = 0;i<parameters->nb_fin_en;i++){
		fprintf(file,"%15.5e",hadronized_energies[i]);
		for(int j = 0;j<parameters->nb_fin_part;j++){
			fprintf(file,"%15.5e",instantaneous_integrated_hadronized_spectra[j][i]*rate_conversion); // conversion from GeV to CGS units
		}
		fprintf(file,"\n");
	}
	fclose(file);
	return;
}

