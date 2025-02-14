// This is the source file where the methods computing
// the primary spectra are implemented.
// Last modification: 15 September 2023
// Authors: Jérémy Auffinger jeremy.auffinger.blackhawk@protonmail.com & Alexandre Arbey alexandre.arbey@ens-lyon.fr

#include "include.h"

void read_gamma_tables(double ***gammas,double *gamma_param,double *gamma_x,struct param *parameters){
	// This function reads the greybody factors tables stored in the folder '/gamma_tables'.
	// It fills the tabulated BH spins, energies and factors in the arrays
	// gamma_param[], gamma_x[] and gammas[][][] respectively.
	
	char **gamma_names = (char **)malloc(5*sizeof(char *)); // contains the names of the greybody tables files
	for(int i = 0;i<parameters->nb_gamma_spins;i++){
		gamma_names[i] = (char *)malloc(64*sizeof(char));
	}
	switch(parameters->metric){
		case 0:{ // Kerr BHs
			sprintf(gamma_names[0],"./src/tables/gamma_tables/spin_0.txt");
			sprintf(gamma_names[1],"./src/tables/gamma_tables/spin_1.txt");
			sprintf(gamma_names[2],"./src/tables/gamma_tables/spin_2.txt");
			sprintf(gamma_names[3],"./src/tables/gamma_tables/spin_0.5.txt");
			sprintf(gamma_names[4],"./src/tables/gamma_tables/spin_1.5.txt");
			break;
		}
		case 1:{ // polymerized BHs
			if(parameters->epsilon_LQG > 0.79){
				sprintf(gamma_names[0],"./src/tables/gamma_tables/LQG/spin_0_LQG_highe.txt");
				sprintf(gamma_names[1],"./src/tables/gamma_tables/LQG/spin_1_LQG_highe.txt");
				sprintf(gamma_names[2],"./src/tables/gamma_tables/LQG/spin_2_LQG_highe.txt");
				sprintf(gamma_names[3],"./src/tables/gamma_tables/LQG/spin_0.5_LQG_highe.txt");
			}
			else{
				if(parameters->a0_LQG == 0.){
					sprintf(gamma_names[0],"./src/tables/gamma_tables/LQG/spin_0.txt");
					sprintf(gamma_names[1],"./src/tables/gamma_tables/LQG/spin_1.txt");
					sprintf(gamma_names[2],"./src/tables/gamma_tables/LQG/spin_2.txt");
					sprintf(gamma_names[3],"./src/tables/gamma_tables/LQG/spin_0.5.txt");
				}
				else{
					sprintf(gamma_names[0],"./src/tables/gamma_tables/LQG/spin_0_a0.txt");
					sprintf(gamma_names[1],"./src/tables/gamma_tables/LQG/spin_1_a0.txt");
					sprintf(gamma_names[2],"./src/tables/gamma_tables/LQG/spin_2_a0.txt");
					sprintf(gamma_names[3],"./src/tables/gamma_tables/LQG/spin_0.5_a0.txt");
				}
			}
			break;
		}
		case 2:{ // charged BHs
			sprintf(gamma_names[0],"./src/tables/gamma_tables/charged/spin_0.txt");
			sprintf(gamma_names[1],"./src/tables/gamma_tables/charged/spin_1.txt");
			sprintf(gamma_names[2],"./src/tables/gamma_tables/charged/spin_2.txt");
			sprintf(gamma_names[3],"./src/tables/gamma_tables/charged/spin_0.5.txt");
			break;
		}
		case 3:{ // higher-dimensional BHs
			sprintf(gamma_names[0],"./src/tables/gamma_tables/higher/spin_0.txt");
			sprintf(gamma_names[1],"./src/tables/gamma_tables/higher/spin_1.txt");
			sprintf(gamma_names[2],"./src/tables/gamma_tables/higher/spin_2.txt");
			sprintf(gamma_names[3],"./src/tables/gamma_tables/higher/spin_0.5.txt");
			break;
		}
		default:{
			printf("\n\t [read_gamma_tables] : ERROR WRONG BH METRIC !\n");
			fflush(stdout);
			exit(0);
			break;
		}
	}
	
	char dummy[32];
	FILE *table;
	if(parameters->full_output){
		printf("\n\n\t reading:\t");
		fflush(stdout);
	}
	for(int i = 0;i<parameters->nb_gamma_spins;i++){
		table = fopen(gamma_names[i],"r");
		if(!table){
			printf("\n\t [read_gamma_tables] : ERROR COULD NOT FIND TABLE '%s' !\n",gamma_names[i]);
			fflush(stdout);
			exit(0);
		}
		rewind(table);
		fscanf(table,"%s",dummy);
		for(int j = 0;j<parameters->nb_gamma_x;j++){
			fscanf(table,"%lf",&(gamma_x[j]));
		}
		for(int j = 0;j<parameters->nb_gamma_param;j++){
			fscanf(table,"%lf",&(gamma_param[j]));
			for(int k = 0;k<parameters->nb_gamma_x;k++){
				fscanf(table,"%lf",&(gammas[i][j][k]));
			}
		}
		fclose(table);
		if(parameters->full_output){ // prints the state of reading
			if(i == 0){
				printf("%s\n",gamma_names[i]);
			}
			else{
				printf("\t\t\t%s\n",gamma_names[i]);
			}
			fflush(stdout);
		}
	}
	if(parameters->full_output){
		printf("\n");
		fflush(stdout);
	}
	free2D_char(gamma_names,parameters->nb_gamma_spins);
	return;
}

void read_asymp_fits(double ***fits,struct param *parameters){
	// This function reads the tables fitting asymptotic Kerr emissivities
	// in the folder "/gamma_tables". It fills the array fits[][][].
	
	char **table_names;
	table_names = (char **)malloc(parameters->nb_gamma_spins*sizeof(char *)); // contains the names of the fits tables files
	for(int i = 0;i<parameters->nb_gamma_spins;i++){
		table_names[i] = (char *)malloc(64*sizeof(char));
	}
	switch(parameters->metric){
		case 0:{ // Kerr BHs
			sprintf(table_names[0],"./src/tables/gamma_tables/spin_0_fits.txt");
			sprintf(table_names[1],"./src/tables/gamma_tables/spin_1_fits.txt");
			sprintf(table_names[2],"./src/tables/gamma_tables/spin_2_fits.txt");
			sprintf(table_names[3],"./src/tables/gamma_tables/spin_0.5_fits.txt");
			sprintf(table_names[4],"./src/tables/gamma_tables/spin_1.5_fits.txt");
			break;
		}
		case 1:{ // polymerized BHs
			if(parameters->epsilon_LQG > 0.79){
				sprintf(table_names[0],"./src/tables/gamma_tables/LQG/spin_0_fits_highe.txt");
				sprintf(table_names[1],"./src/tables/gamma_tables/LQG/spin_1_fits_highe.txt");
				sprintf(table_names[2],"./src/tables/gamma_tables/LQG/spin_2_fits_highe.txt");
				sprintf(table_names[3],"./src/tables/gamma_tables/LQG/spin_0.5_fits_highe.txt");
			}
			else{
				if(parameters->a0_LQG == 0.){
					sprintf(table_names[0],"./src/tables/gamma_tables/LQG/spin_0_fits.txt");
					sprintf(table_names[1],"./src/tables/gamma_tables/LQG/spin_1_fits.txt");
					sprintf(table_names[2],"./src/tables/gamma_tables/LQG/spin_2_fits.txt");
					sprintf(table_names[3],"./src/tables/gamma_tables/LQG/spin_0.5_fits.txt");
				}
				else{
					sprintf(table_names[0],"./src/tables/gamma_tables/LQG/spin_0_fits_a0.txt");
					sprintf(table_names[1],"./src/tables/gamma_tables/LQG/spin_1_fits_a0.txt");
					sprintf(table_names[2],"./src/tables/gamma_tables/LQG/spin_2_fits_a0.txt");
					sprintf(table_names[3],"./src/tables/gamma_tables/LQG/spin_0.5_fits_a0.txt");
				}
			}
			break;
		}
		case 2:{ // charged BHs
			sprintf(table_names[0],"./src/tables/gamma_tables/charged/spin_0_fits.txt");
			sprintf(table_names[1],"./src/tables/gamma_tables/charged/spin_1_fits.txt");
			sprintf(table_names[2],"./src/tables/gamma_tables/charged/spin_2_fits.txt");
			sprintf(table_names[3],"./src/tables/gamma_tables/charged/spin_0.5_fits.txt");
			break;
		}
		case 3:{ // higher dimensional BHs
			sprintf(table_names[0],"./src/tables/gamma_tables/higher/spin_0_fits.txt");
			sprintf(table_names[1],"./src/tables/gamma_tables/higher/spin_1_fits.txt");
			sprintf(table_names[2],"./src/tables/gamma_tables/higher/spin_2_fits.txt");
			sprintf(table_names[3],"./src/tables/gamma_tables/higher/spin_0.5_fits.txt");
			break;
		}
		default:{
			printf("\n\t [read_asymp_fits] : ERROR WRONG BH METRIC !\n");
			fflush(stdout);
			exit(0);
			break;
		}
	}
	
	FILE *file;
	char dummy[32];
	if(parameters->full_output){
		printf("\n\n\t reading:\t");
		fflush(stdout);
	}
	for(int i = 0;i<parameters->nb_gamma_spins;i++){
		file = fopen(table_names[i],"r");
		if(file == 0){
			printf("\n\t [read_asymp_fits] : ERROR COULD NOT OPEN FILE '%s' !\n",table_names[i]);
			fflush(stdout);
			exit(0);
		}
		rewind(file);
		for(int j = 0;j<1+parameters->nb_gamma_fits;j++){
			fscanf(file,"%s",dummy);
		}
		for(int j = 0;j<parameters->nb_gamma_param;j++){
			fscanf(file,"%s",dummy);
			for(int k = 0;k<parameters->nb_gamma_fits;k++){
				fscanf(file,"%lf",&(fits[i][j][k]));
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
		fclose(file);
	}
	if(parameters->full_output){
		printf("\n");
		fflush(stdout);
	}
	free2D_char(table_names,parameters->nb_gamma_spins);
	return;
}

double dNdtdE(double E,double M,double param,int particle_index,double ***gammas,double *gamma_param,double *gamma_x,double ***fits,double *dof,double *spins,double *masses_primary,int counter_param,int counter_x,struct param *parameters){
	// This function computes the emission rate d²N/dtdE for the elementary particle of type
	// particle_index. It uses the greybody factors tables or the asymptotical fits.
	
	if(M == 0. || (parameters->BH_remnant && M == parameters->M_remnant)){ // the BH has evaporated or is a stable remnant
		return 0.;
	}
	double T;
	double x = 2.*E*M*G;
	switch(parameters->metric){
		case 0:{ // Kerr BHs
			T = temp_Kerr(M,param);
			break;
		}
		case 1:{ // polymerized BHs
			T = temp_LQG(M,param,parameters->a0_LQG);
			break;
		}
		case 2:{ // charged BHs
			T = temp_charged(M,param);
			break;
		}
		case 3:{ // higher BHs
			T = temp_higher(M,param,parameters->M_star);
			break;
		}
		default:{
			printf("\n\t [dNdtdE] : ERROR WRONG BH METRIC !\n");
			fflush(stdout);
			exit(0);
			break;
		}
	}
	int type; // contains the spin type of the particle
	if(spins[particle_index] == 0.){
		type = 0;
	}
	if(spins[particle_index] == 1.){
		type = 1;
	}
	if(spins[particle_index] == 2.){
		type = 2;
	}
	if(spins[particle_index] == 0.5){
		type = 3;
	}
	if(spins[particle_index] == 1.5){
		type = 4;
	}
	if(E < masses_primary[particle_index]){
		return 0.;
	}
	else if(x > 0.01 && x < 5.){ // tabulated value
		switch(parameters->interpolation_method){
			case 0:{ // linear interpolation
				return dof[particle_index]*(gammas[type][counter_param-1][counter_x-1] + (gammas[type][counter_param][counter_x-1] - gammas[type][counter_param-1][counter_x-1])/(gamma_param[counter_param] - gamma_param[counter_param-1])*(param - gamma_param[counter_param-1]) + (gammas[type][counter_param-1][counter_x] - gammas[type][counter_param-1][counter_x-1])/(gamma_x[counter_x] - gamma_x[counter_x-1])*(2.*E*M*G - gamma_x[counter_x-1]))/(2.*pi);
				break;
			}
			case 1:{ // logarithmic interpolation
				if(gammas[type][counter_param-1][counter_x-1] == 0. || gammas[type][counter_param][counter_x-1] == 0. || gammas[type][counter_param-1][counter_x] == 0.){
					return 0.;
				}
				else{
					return dof[particle_index]*(pow(10.,log10(gammas[type][counter_param-1][counter_x-1]) + (log10(gammas[type][counter_param-1][counter_x]) - log10(gammas[type][counter_param-1][counter_x-1]))/(log10(gamma_x[counter_x]) - log10(gamma_x[counter_x-1]))*(log10(2.*E*M*G) - log10(gamma_x[counter_x-1]))) + (gammas[type][counter_param][counter_x-1] - gammas[type][counter_param-1][counter_x-1])/(gamma_param[counter_param] - gamma_param[counter_param-1])*(param - gamma_param[counter_param-1]))/(2.*pi);
				}
				break;
			}
			default:{
				printf("\n\t [dNdtdE] : ERROR WRONG INTERPOLATION METHOD !\n");
				fflush(stdout);
				exit(0);
				break;
			}
		}
	}
	else if(x <= 0.01){ // low-energy asymptotic value
		double a1;
		double b1;
		a1 = fits[type][counter_param-1][0] + (fits[type][counter_param][0] - fits[type][counter_param-1][0])/(gamma_param[counter_param] - gamma_param[counter_param-1])*(param - gamma_param[counter_param-1]);
		b1 = fits[type][counter_param-1][1] + (fits[type][counter_param][1] - fits[type][counter_param-1][1])/(gamma_param[counter_param] - gamma_param[counter_param-1])*(param - gamma_param[counter_param-1]);
		if(parameters->metric == 0){
			return dof[particle_index]*pow(10.,a1*log10(x) + b1)/(2.*pi);
		}
		else{
			if(type <= 2){ // bosons
				return dof[particle_index]*pow(10.,a1*log10(x) + b1)/exp_adapt(E/T)/(pi);
			}
			else{ // fermions
				return dof[particle_index]*pow(10.,a1*log10(x) + b1)/(exp(E/T) + 1.)/(pi);
			}
		}
	}
	else{ // high-energy asymptotic value x >= 5.
		if(parameters->metric == 0){
			double a1;
			double b1;
			double c1;
			double d1;
			double e1;
			a1 = fits[type][counter_param-1][2] + (fits[type][counter_param][2] - fits[type][counter_param-1][2])/(gamma_param[counter_param] - gamma_param[counter_param-1])*(param - gamma_param[counter_param-1]);
			b1 = fits[type][counter_param-1][3] + (fits[type][counter_param][3] - fits[type][counter_param-1][3])/(gamma_param[counter_param] - gamma_param[counter_param-1])*(param - gamma_param[counter_param-1]);
			c1 = fits[type][counter_param-1][4] + (fits[type][counter_param][4] - fits[type][counter_param-1][4])/(gamma_param[counter_param] - gamma_param[counter_param-1])*(param - gamma_param[counter_param-1]);
			d1 = fits[type][counter_param-1][5] + (fits[type][counter_param][5] - fits[type][counter_param-1][5])/(gamma_param[counter_param] - gamma_param[counter_param-1])*(param - gamma_param[counter_param-1]);
			e1 = fits[type][counter_param-1][6] + (fits[type][counter_param][6] - fits[type][counter_param-1][6])/(gamma_param[counter_param] - gamma_param[counter_param-1])*(param - gamma_param[counter_param-1]);
			return dof[particle_index]*pow(10.,a1*(2.*M*E*G) + b1 + c1*cos(e1*2.*E*M*G) + d1*sin(e1*2.*E*M*G))/(2.*pi);
		}
		else{
			double c1;
			c1 = fits[type][counter_param-1][2] + (fits[type][counter_param][2] - fits[type][counter_param-1][2])/(gamma_param[counter_param] - gamma_param[counter_param-1])*(param - gamma_param[counter_param-1]);
			if(type <= 2){ // bosons
				return dof[particle_index]*(27./4.*pow(x,2.)/exp_adapt(E/T)*c1)/(2.*pi);
			}
			else{ // fermions
				return dof[particle_index]*(27./4.*pow(x,2.)/(exp(E/T) + 1.)*c1)/(2.*pi);
			}
		}
	}
}

void instantaneous_primary_spectrum(double **instantaneous_primary_spectra,double **BH_masses,double **BH_params,double **spec_table,double *energies,
	double ***gammas,double *gamma_param,double *gamma_x,double ***fits,double *dof,double *spins,double *masses_primary,int **counters_param, int ***counters_x,
	int *compute,struct param *parameters){
	// This function computes the instantaneous primary spectra for all particles thanks to the function
	// dNdtdE and the BH spectrum at some timestep.
	
	for(int i = 0;i<parameters->BH_number;i++){
		for(int j = 0;j<parameters->param_number;j++){
			counters_param[i][j] = 0;
			while(counters_param[i][j] < parameters->nb_gamma_param - 1 && BH_params[i][j] >= gamma_param[counters_param[i][j]]){ // finding the nearest parameter
				counters_param[i][j]++;
			}
		}
		for(int j = 0;j<parameters->param_number;j++){
			counters_x[i][j][0] = 0;
			for(int k = 0;k<parameters->E_number;k++){
				if(k>0){
					counters_x[i][j][k] = counters_x[i][j][k-1];
				}
				while(counters_x[i][j][k] < parameters->nb_gamma_x - 1 && energies[k]*2.*BH_masses[i][j]*G >= gamma_x[counters_x[i][j][k]]){ // finding the nearest omega*r_H
					counters_x[i][j][k]++;
				}
			}
		}
	}	
	if(parameters->full_output){
		printf("\n\n\t compute:\t");
		fflush(stdout);
	}
	for(int i = 0;i<parameters->particle_number+parameters->grav+parameters->add_DM;i++){
/*#if defined(_OPENMP)
#pragma omp parallel for
#endif*/
		for(int j = 0;j<parameters->E_number;j++){
			instantaneous_primary_spectra[i][j] = 0.;
			if(compute[i]){
				for(int k = 0;k<parameters->BH_number;k++){
					for(int l = 0;l<parameters->param_number;l++){
						instantaneous_primary_spectra[i][j] += spec_table[k][l]*dNdtdE(energies[j],BH_masses[k][l],BH_params[k][l],i,gammas,gamma_param,gamma_x,fits,dof,spins,masses_primary,counters_param[k][l],counters_x[k][l][j],parameters);
					}
				}
			}
		}
		if(parameters->full_output){ // prints the current state of computation
			printf("%i/%i ",i+1,parameters->particle_number+parameters->grav+parameters->add_DM);
			fflush(stdout);
		}
	}
	if(parameters->full_output){
		printf("\n\n");
		fflush(stdout);
	}
	return;
}

void write_instantaneous_primary_spectra(double **instantaneous_primary_spectra,double *energies,struct param *parameters){
	// This function writes the instantaneous Hawking primary spectrum of a distribution of BH,
	// contained in the array instantaneous_primary_spectra[][], into the file 'instantaneous_primary_spectra.txt'.
	
	char file_name[500];
	sprintf(file_name,"./results/%s/instantaneous_primary_spectra.txt",parameters->destination_folder);
	FILE *file = fopen(file_name,"w+");
	if(!file){
		printf("\n\t [write_instantaneous_primary_spectra] : ERROR COULD NOT OPEN FILE '%s' !\n",file_name);
		fflush(stdout);
		exit(0);
	}
	rewind(file);
	fprintf(file,"Hawking primary spectra for each particle types.\n");
	if(!parameters->grav){
		if(parameters->add_DM){
			if(parameters->hadronization_choice != 3){
				fprintf(file,"%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n","energy/particle","photon","gluons","higgs","W+-","Z0","neutrinos","electron","muon","tau","up","down","charm","strange","top","bottom","DM");
			}
			else{
				fprintf(file,"%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n","energy/particle","photon","gluons","higgs","W+-","Z0","neutrinos","electron","muon","tau","up","down","charm","strange","top","bottom","pions");
			}
		}
		else{
			fprintf(file,"%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n","energy/particle","photon","gluons","higgs","W+-","Z0","neutrinos","electron","muon","tau","up","down","charm","strange","top","bottom");
		}
	}
	else{
		if(parameters->add_DM){
			fprintf(file,"%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n","energy/particle","photon","gluons","higgs","W+-","Z0","neutrinos","electron","muon","tau","up","down","charm","strange","top","bottom","graviton","DM");
		}
		else{
			fprintf(file,"%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n","energy/particle","photon","gluons","higgs","W+-","Z0","neutrinos","electron","muon","tau","up","down","charm","strange","top","bottom","graviton");
		}
	}
	for(int j = 0;j<parameters->E_number;j++){
		fprintf(file,"%15.5e",energies[j]);
		for(int i = 0;i<parameters->particle_number+parameters->grav+parameters->add_DM;i++){
			fprintf(file,"%15.5e",instantaneous_primary_spectra[i][j]*rate_conversion); // conversion from GeV to CGS units
		}
		fprintf(file,"\n");
	}
	fclose(file);
	return;
}

