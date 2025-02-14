// This is the source file where the methods computing
// the initial distribution of black holes are implemented.
// Last modification: 15 September 2023
// Authors: Jérémy Auffinger jeremy.auffinger.blackhawk@protonmail.com & Alexandre Arbey alexandre.arbey@ens-lyon.fr

#include "include.h"

void read_users_table(double *init_masses,double *init_params,double **spec_table,struct param *parameters){
	// This function reads the user's BH distribution in the file parameters->table.
	
	char table_name[100];
	sprintf(table_name,"./src/tables/users_spectra/%s",parameters->table);
	FILE *table = fopen(table_name,"r+");
	if(!table){
		printf("\n\t [read_user_table] : ERROR COULD NOT OPEN FILE '%s' !",parameters->table);
		fflush(stdout);
		exit(0);
	};
	rewind(table);
	char dummy[64];
	fscanf(table,"%s",dummy);
	for(int i = 0;i<parameters->param_number;i++){
		fscanf(table,"%lf",&(init_params[i]));
	}
	for(int i = 0;i<parameters->BH_number;i++){
		fscanf(table,"%lf",&(init_masses[i]));
		init_masses[i] = init_masses[i]*mass_conversion; // conversion from CGS to GeV
		for(int j = 0;j<parameters->param_number;j++){
			fscanf(table,"%lf",&(spec_table[i][j]));
			spec_table[i][j] = spec_table[i][j]/dens_conversion; // conversion from CGS to GeV
		}
	}
	fclose(table);
	return;
}

double nu(double M){
	// This function computes the intermediate quantity nu(M) (Tashiro & Sugiyama 2008)
	
	double spec_ind = 1.3; // spectral exponent characterizing the power spectrum of primordial fluctuations
	double gamma_ind = 6.2202728740498765; // Gamma((spec_ind-1)/2), Riemann function
	double k0 = (0.002/(3.085678e+22*leng_conversion*100.)); // typical length scale of primordial fluctuations
	double Rc = 24e-10; // amplitude of primordial fluctuations
	double zeta = 0.7; // (or 1.2) theoretical threshold for direct collapse
	double H0 = (67800./(3.085678e+22*time_conversion)); // current Hubble parameter (PDG 2017)
	double OmegaM = 0.308; // matter component in the Universe (PDG 2017)
	double zeq = 3200; // redshift of the radiation/matter equality
	double geq = 3.36; // number of entropic ultra-relativistic degrees of freedom at radiation/matter equality
	double g = 106.75; // initial number of entropic ultra-relativistic degrees of freedom (here at the end of inflation)
	double X = 1./(2.*G)*sqrt(pow(H0,2.)*OmegaM/(1.+zeq)*pow(geq/g,1./3.));
	return sqrt(2.*pow(pow(k0,2.)*M/X,(spec_ind-1.)/2.)/(Rc*gamma_ind))*zeta;
}

double M_dist(double M,struct param *parameters){
	// This function computes the  marginal distribution -dn/dM(M) of BH
	// obtained through peak-theory (case spectrum_choice = 4, Tashiro &
	// Sugiyama 2008), lognormal distribution for the mass density
	// (case spectrum_choice = 1), or for the number density (case
	// spectrum_choice = 11), power-law distribution (case spectrum_choice = 2), critical
	// collapse theory (case spectrum_choice = 3) (cases 1-3 are described
	// in Carr, Raidal, Tenkanen, Vaskonen & Veermäe 2017) or uniform
	// distribution (case spectrum_choice = 5). The case
	// spectrum_choice = 0 describes a Dirac distribution.
	
	switch(parameters->spectrum_choice){
		case 0:{ // in this case we use a Dirac distribution
			return 1./dens_conversion;
			break;
		}
		case 1:{ // in this case we use a log-normal distribution for the mass density
			return parameters->amplitude/(sqrt(2.*pi)*parameters->stand_dev*pow(M,2.))*exp(-pow(log(M/parameters->crit_mass),2.)/(2.*pow(parameters->stand_dev,2.)));
			break;
		}
		case 11:{ // in this case we use a log-normal distribution for the number density
			return parameters->amplitude/(sqrt(2.*pi)*parameters->stand_dev*pow(M,1.))*exp(-pow(log(M/parameters->crit_mass),2.)/(2.*pow(parameters->stand_dev,2.)));
			break;
		}
		case 2:{ // in this case we use a power-law distribution
			double gamma = -2.*parameters->eqstate/(1.+parameters->eqstate);
			return parameters->amplitude*pow(M,gamma-2.);
			break;
		}
		case 3:{ // in this case we use a critical collapse distribution
			return parameters->amplitude*pow(M,1.85)*exp(-pow((M/parameters->crit_mass),2.85));
			break;
		}
		case 4:{ // in this case we use the peak-theory distribution
			double spec_ind = 1.3; // spectral exponent characterizing the power spectrum of primordial fluctuations
			double H0 = (67800./(3.085678e+22*time_conversion)); // current Hubble parameter (PDG 2017)
			double OmegaM = 0.308; // matter component in the Universe (PDG 2017)
			double zeq = 3200; // redshift of the radiation/matter equality
			double geq = 3.36; // number of entropic ultra-relativistic degrees of freedom at radiation/matter equality
			double g = 106.75; // initial number of entropic ultra-relativistic degrees of freedom (here at the end of inflation)
			double X = 1./(2.*G)*sqrt(pow(H0,2.)*OmegaM/(1.+zeq)*pow(geq/g,1./3.));
			return 1./(4.*pow(pi,2.)*M)*pow(X*(spec_ind-1.)/(6.*M),3./2.)*(spec_ind-1.)/2.*pow(nu(M),4)*exp(-pow(nu(M),2.)/2.);
			break;
		}
		case 5:{ // in this case we use a uniform distribution
			return parameters->amplitude/parameters->BH_number/(pow(10.,log10(M)+0.5*(log10(parameters->Mmax) - log10(parameters->Mmin))/(parameters->BH_number-1)) - pow(10.,log10(M)-0.5*(log10(parameters->Mmax) - log10(parameters->Mmin))/(parameters->BH_number-1)));
		}
		default:{
			printf("\n\t [M_dist] : ERROR WRONG MASS SPECTRUM CHOICE !\n");
			fflush(stdout);
			exit(0);
			break;
		}
	}
}

double param_dist(double param,struct param *parameters){
	// This function computes the marginal parameter distribution of black holes.
	// spectrum_choice_param = 0 denotes a Dirac distribution.
	
	switch(parameters->spectrum_choice_param){
		case 0:{ // in this case we use a Dirac distribution
			return 1.;
			break;
		}
		case 1:{ // in this case we use a uniform distribution
			if(parameters->metric == 0){
				return 1./parameters->param_number/(parameters->amax - parameters->amin)*(parameters->param_number - 1);
			}
			else if(parameters->metric == 2){
				return 1./parameters->param_number/(parameters->Qmax - parameters->Qmin)*(parameters->param_number - 1);
			}
			break;
		}
		case 2:{ // in this case we use a gaussian distribution
			return 1./sqrt(2.*pi*parameters->stand_dev_param)*exp(-pow(param - parameters->mean_param,2.)/(2.*parameters->stand_dev_param));
			break;
		}
		default:{
			printf("\n\t [param_dist] : ERROR WRONG PARAMETER SPECTRUM CHOICE !\n");
			fflush(stdout);
			exit(0);
			break;
		}
	}
}

void spectrum(double *init_masses,double *init_params,double **spec_table,struct param *parameters){
	// This function fills the init_masses[] with a logarithmic distribution of masses,
	// fills the init_params[] with a linear distribution of spins
	// and computes the number codensity dn(M,a) of PBH in intervals [M,M+dM], [x,x+dx]
	// (where x = a^* or Q^*) in the array spec_table[][]. It adds a correction of 1.e+100 due to the
	// very small numbers at stake (unit conversion).
	
	if(parameters->BH_number != 1){
		for(int i = 0;i<parameters->BH_number;i++){
			init_masses[i] = pow(10.,log10(parameters->Mmin) + i*(log10(parameters->Mmax) - log10(parameters->Mmin))/(parameters->BH_number-1));
		}
	}
	else{
		init_masses[0] = parameters->Mmin;
	}
	if(parameters->param_number != 1){
		for(int i = 0;i<parameters->param_number;i++){
			if(parameters->metric == 0){
				init_params[i] = parameters->amin + (parameters->amax - parameters->amin)/(parameters->param_number - 1)*i;
			}
			else if(parameters->metric == 2){
				init_params[i] = parameters->Qmin + (parameters->Qmin - parameters->Qmax)/(parameters->param_number - 1)*i;
			}
		}
	}
	else{
		switch(parameters->metric){
			case 0:{ // Kerr BHs
				init_params[0] = parameters->amin;
				break;
			}
			case 1:{ // polymerized BHs
				init_params[0] = parameters->epsilon_LQG;
				break;
			}
			case 2:{ // charged BHs
				init_params[0] = parameters->Qmin;
				break;
			}
			case 3:{ // higher dimensional BHs
				init_params[0] = parameters->n;
				break;
			}
			default:{
				printf("\n\t [spectrum] : ERROR WRONG BH METRIC !\n");
				fflush(stdout);
				exit(0);
				break;
			}
		}
	}
	if(parameters->BH_number != 1 && parameters->param_number != 1){
		for(int i = 0;i<parameters->BH_number;i++){
			for(int j = 0;j<parameters->param_number;j++){
				spec_table[i][j] = 1.e+100*M_dist(init_masses[i],parameters)*(pow(10.,log10(init_masses[i])+0.5*(log10(parameters->Mmax) - log10(parameters->Mmin))/(parameters->BH_number-1)) - pow(10.,log10(init_masses[i])-0.5*(log10(parameters->Mmax) - log10(parameters->Mmin))/(parameters->BH_number-1))); // WARNING the density of black holes is here rescaled from 1e-100 to 1 !
				if(parameters->metric == 0){
					spec_table[i][j] = spec_table[i][j]*param_dist(init_params[j],parameters)*(parameters->amax - parameters->amin)/(parameters->param_number - 1);
				}
				else if(parameters->metric == 2){
					spec_table[i][j] = spec_table[i][j]*param_dist(init_params[j],parameters)*(parameters->Qmax - parameters->Qmin)/(parameters->param_number - 1);
				}
			}
		}
	}
	else if(parameters->BH_number != 1){
		for(int i = 0;i<parameters->BH_number;i++){
			spec_table[i][0] = 1.e+100*M_dist(init_masses[i],parameters)*(pow(10.,log10(init_masses[i])+0.5*(log10(parameters->Mmax) - log10(parameters->Mmin))/(parameters->BH_number-1)) - pow(10.,log10(init_masses[i])-0.5*(log10(parameters->Mmax) - log10(parameters->Mmin))/(parameters->BH_number-1)));
			spec_table[i][0] = spec_table[i][0]*param_dist(init_params[0],parameters);
		}
	}
	else if(parameters->param_number != 1){
		for(int j = 0;j<parameters->param_number;j++){
			spec_table[0][j] = M_dist(init_masses[0],parameters);
			if(parameters->metric == 0){
				spec_table[0][j] = spec_table[0][j]*param_dist(init_params[j],parameters)*(parameters->amax - parameters->amin)/(parameters->param_number - 1);
			}
			else if(parameters->metric == 2){
				spec_table[0][j] = spec_table[0][j]*param_dist(init_params[j],parameters)*(parameters->Qmax - parameters->Qmin)/(parameters->param_number - 1);
			}
		}
	}
	else{
		spec_table[0][0] = M_dist(init_masses[0],parameters);
		spec_table[0][0] = spec_table[0][0]*param_dist(init_params[0],parameters);
	}
	return;
}

void write_spectrum(double *init_masses,double *init_params,double **spec_table,struct param *parameters){
	// This function writes the BH initial masses, parameters and densities
	// contained in the arrays spec_table[], init_masses[] and init_params[] in
	// the file 'BH_spectrum.txt'.
	
	char file_name[500];
	sprintf(file_name,"./results/%s/BH_spectrum.txt",parameters->destination_folder);
	FILE *file = fopen(file_name,"w+");
	if(!file){
		printf("\n\t [write_spectrum] : ERROR COULD NOT OPEN FILE '%s' !\n",file_name);
		fflush(stdout);
		exit(0);
	}
	rewind(file);
	fprintf(file,"Initial BH comoving number density as a function of their mass and parameter.\n\n");
	switch(parameters->metric){
		case 0:{ // Kerr BHs
			fprintf(file,"%15s","mass/spin");
			break;
		}
		case 1:{ // polymerized BHs
			fprintf(file,"%15s","mass/epsilon");
			break;
		}
		case 2:{ // charged BHs
			fprintf(file,"%15s","mass/Q");
			break;
		}
		case 3:{ // higher dimensional BHs
			fprintf(file,"%15s","mass/n");
			break;
		}
		default:{
			printf("\n\t [write_spectrum] : ERROR WRONG BH METRIC !\n");
			fflush(stdout);
			exit(0);
			break;
		}
	}
	for(int i = 0;i<parameters->param_number;i++){
		fprintf(file,"%15.5e",init_params[i]);
	}
	fprintf(file,"\n");
	for(int i = 0;i<parameters->BH_number;i++){
		fprintf(file,"%15.5e",init_masses[i]/mass_conversion); // conversion from GeV to CGS units
		for(int j = 0;j<parameters->param_number;j++){
			fprintf(file,"%15.5e",spec_table[i][j]*dens_conversion); // conversion from GeV to CGS units
		}
		fprintf(file,"\n");
	}
	fclose(file);
	return;
}
