// Program that computes the Page f(M,a*) and g(M,a*) factors of a Kerr Black Hole,
// or the f(M,epsilon) factor of a polymerized BH
// Author : Jérémy Auffinger, j.auffinger@ipnl.in2p3.fr
// Last modification : 21 April 2021

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

int nb_particles = 15; // nb of SM particles
int grav = 0; // additional graviton
int ADD_DM = 1; // additional DM
double spin_DM = 0.0;


double m_DM = 0.0003;
double dof_DM = 100000.0;


double pi = 3.14159265389795;
double Mp = 1.221e+19; // Planck mass in GeV
int nb_a = 50; // nb of spin values in the Kerr greybody factor tables
int nb_epsilon = 11; // nb of epsilon values in the polymerized greybody factor tables
int nb_x = 200; // nb of energies in the greybody factor tables
int nb_particle_spins = 5; // nb of particle spins (5 for Kerr BHs, 4 for other metrics for which the spin 3/2 has not been evaluated)
int nb_masses = 1000; // nb of masses to tabulate
double mass_conversion = 5.60958884e+23;
double time_conversion = 1.519267407e+24;
double Mmax = 2.e+56; // maximum mass to tabulate, in GeV (minimum is the Planck mass)
double a0 = 0.; // the polymerized minimum area

#define G (1./pow(Mp,2.)) // Newton constant in GeV

double expbis(double x){
	if(x > 1.e-5){
		return exp(x) - 1.;
	}
	else{
		return x + pow(x,2.)/2. + pow(x,3.)/6. + pow(x,4.)/24. + pow(x,5.)/120.;
	};
};

double rplus_Kerr(double M,double a){
	// This function computes the external Kerr radius
	// of a BH of mass M and spin parameter a.
	
	return M*(1. + sqrt(1. - pow(a,2.)));
};

double temp_Kerr(double M,double a){
	// This function computes the Hawking temperature of a Kerr
	// PBH of mass M and spin parameter a.
	
	return 1./(2.*pi*G)*(rplus_Kerr(M,a) - M)/(pow(rplus_Kerr(M,a),2.) + pow(a*M,2.));
};

double P_LQG(double epsilon){
	// This function computes the polymerization fonction P.
	
	return (sqrt(1. + pow(epsilon,2.)) - 1.)/(sqrt(1. + pow(epsilon,2.)) + 1.);
};

double m_LQG(double M,double epsilon){
	// This function computes the modified mass m in polymeriezd metrics.
	
	return M/pow(1. + P_LQG(epsilon),2.);
};

double temp_LQG(double M,double epsilon){
	// This function computes the Hawking temperature of a polymerized BH
	// with parameters epsilon and a0
	
	return 4.*pow(m_LQG(M,epsilon),3.)*(1-pow(P_LQG(epsilon),2.))/(32.*pi*pow(m_LQG(M,epsilon),4.) + 2.*pi*pow(a0,2.))/G;
};

void read_gamma_tables_Kerr(double ***gammas,double *gamma_a,double *gamma_x){
	// This function reads the Kerr greybody factors tables.
	
	char **gamma_names = (char **)malloc(nb_particle_spins*sizeof(char *));
	for(int i = 0;i<nb_particle_spins;i++){
		gamma_names[i] = (char *)malloc(64*sizeof(char));
	};
	sprintf(gamma_names[0],"./spin_0_fM.txt");
	sprintf(gamma_names[1],"./spin_1_fM.txt");
	sprintf(gamma_names[2],"./spin_2_fM.txt");
	sprintf(gamma_names[3],"./spin_0.5_fM.txt");
	sprintf(gamma_names[4],"./spin_1.5_fM.txt");
	
	char dummy[32];
	FILE *table;
	for(int i = 0;i<nb_particle_spins;i++){
		table = fopen(gamma_names[i],"r");
		rewind(table);
		if(!table){
			printf("\n[read_simple_tables] : ERROR COULD NOT FIND TABLE %s !\n",gamma_names[i]);
			fflush(stdout);
			return;
		};
		fscanf(table,"%s",dummy);
		for(int j = 0;j<nb_x;j++){
			fscanf(table,"%lf",&(gamma_x[j]));
		};
		for(int j = 0;j<nb_a;j++){// WARNING the spin 3/2 case is a dummy table, it was only computed for the case a* = 0
			fscanf(table,"%lf",&(gamma_a[j]));
			for(int k = 0;k<nb_x;k++){
				fscanf(table,"%lf",&(gammas[i][j][k]));
			};
		};
		fclose(table);
	};
	for(int i = 0;i<nb_particle_spins;i++){
		free(gamma_names[i]);
	};
	return;
};

void read_gamma_tables_LQG(double ***gammas,double *gamma_epsilon,double *gamma_x){
	// This function reads the polymerized greybody factors tables.
	
	char **gamma_names = (char **)malloc(nb_particle_spins*sizeof(char *));
	for(int i = 0;i<nb_particle_spins;i++){
		gamma_names[i] = (char *)malloc(64*sizeof(char));
	};
	if(a0 == 0.){
		sprintf(gamma_names[0],"./spin_0_LQG.txt");
		sprintf(gamma_names[1],"./spin_1_LQG.txt");
		sprintf(gamma_names[2],"./spin_2_LQG.txt");
		sprintf(gamma_names[3],"./spin_0.5_LQG.txt");
	}
	else{
		sprintf(gamma_names[0],"./spin_0_LQG_a0.txt");
		sprintf(gamma_names[1],"./spin_1_LQG_a0.txt");
		sprintf(gamma_names[2],"./spin_2_LQG_a0.txt");
		sprintf(gamma_names[3],"./spin_0.5_LQG_a0.txt");
	};
	
	char dummy[32];
	FILE *table;
	for(int i = 0;i<nb_particle_spins;i++){
		table = fopen(gamma_names[i],"r");
		rewind(table);
		if(!table){
			printf("\n[read_simple_tables] : ERROR COULD NOT FIND TABLE %s !\n",gamma_names[i]);
			fflush(stdout);
			return;
		};
		fscanf(table,"%s",dummy);
		for(int j = 0;j<nb_x;j++){
			fscanf(table,"%lf",&(gamma_x[j]));
		};
		for(int j = 0;j<nb_epsilon;j++){
			fscanf(table,"%lf",&(gamma_epsilon[j]));
			for(int k = 0;k<nb_x;k++){
				fscanf(table,"%lf",&(gammas[i][j][k]));
			};
		};
		fclose(table);
	};
	for(int i = 0;i<nb_particle_spins;i++){
		free(gamma_names[i]);
	};
	return;
};

void read_fM_fits_Kerr(double ***fits){
	// This function reads the tables fitting asymptotic Kerr emissivities.
	
	char **table_names;
	table_names = (char **)malloc(nb_particle_spins*sizeof(char *));
	for(int i = 0;i<nb_particle_spins;i++){
		table_names[i] = (char *)malloc(64*sizeof(char));
	};
	sprintf(table_names[0],"./spin_0_fits_fM.txt");
	sprintf(table_names[1],"./spin_1_fits_fM.txt");
	sprintf(table_names[2],"./spin_2_fits_fM.txt");
	sprintf(table_names[3],"./spin_0.5_fits_fM.txt");
	sprintf(table_names[4],"./spin_1.5_fits_fM.txt"); // WARNING this is a dummy table: the spin 3/2 high and low energy fits have only been evaluated for a* = 0
	
	FILE *file;
	char dummy[32];
	for(int i = 0;i<nb_particle_spins;i++){
		file = fopen(table_names[i],"r");
		if(file == 0){
			printf("[read_asymp_fits] : ERROR COULD NOT OPEN FILE %s !\n",table_names[i]);
			fflush(stdout);
			return;
		};
		rewind(file);
		for(int j = 0;j<8;j++){
			fscanf(file,"%s",dummy);
		};
		for(int j = 0;j<nb_a;j++){
			fscanf(file,"%s",dummy);
			for(int k = 0;k<7;k++){
				fscanf(file,"%lf",&(fits[i][j][k]));
			};
		};
		fclose(file);
	};
	for(int i = 0;i<nb_particle_spins;i++){
		free(table_names[i]);
	};
	return;
};

void read_fits_LQG(double ***fits){
	// This function reads the tables fitting asymptotic polymerized emissivities.
	
	char **table_names;
	table_names = (char **)malloc(nb_particle_spins*sizeof(char *));
	for(int i = 0;i<nb_particle_spins;i++){
		table_names[i] = (char *)malloc(64*sizeof(char));
	};
	if(a0 == 0.){
		sprintf(table_names[0],"./spin_0_fits_LQG.txt");
		sprintf(table_names[1],"./spin_1_fits_LQG.txt");
		sprintf(table_names[2],"./spin_2_fits_LQG.txt");
		sprintf(table_names[3],"./spin_0.5_fits_LQG.txt");
	}
	else{
		sprintf(table_names[0],"./spin_0_fits_LQG_a0.txt");
		sprintf(table_names[1],"./spin_1_fits_LQG_a0.txt");
		sprintf(table_names[2],"./spin_2_fits_LQG_a0.txt");
		sprintf(table_names[3],"./spin_0.5_fits_LQG_a0.txt");
	};
	
	FILE *file;
	char dummy[32];
	for(int i = 0;i<nb_particle_spins;i++){
		file = fopen(table_names[i],"r");
		if(file == 0){
			printf("[read_asymp_fits] : ERROR COULD NOT OPEN FILE %s !\n",table_names[i]);
			fflush(stdout);
			return;
		};
		rewind(file);
		for(int j = 0;j<4;j++){
			fscanf(file,"%s",dummy);
		};
		for(int j = 0;j<nb_epsilon;j++){
			fscanf(file,"%s",dummy);
			for(int k = 0;k<3;k++){
				fscanf(file,"%lf",&(fits[i][j][k]));
			};
		};
		fclose(file);
	};
	for(int i = 0;i<nb_particle_spins;i++){
		free(table_names[i]);
	};
	return;
};

double nb_spin_0(double E,double *ms){
	// This function computes the number of spin 0 dof.
	
	double result = 0.;
	if(ADD_DM && spin_DM == 0. && E > m_DM){
		result += dof_DM; // ADD
	};
	if(E > ms[5]){ // Higgs contribution
		result += 1.; // Higgs boson (1 helicity state)
	};
	return result;
};

double nb_spin_1(double E,double *ms){
	// This function computes the number of spin 1 dof.
	
	double result = 2.; // photon contribution (2 helicity states)
	if(ADD_DM && spin_DM == 1. && E > m_DM){
		result += dof_DM; // ADD
	};
	if(E > ms[13]){ // W+- boson contribution
		result += 6.; // W+- bosons (2 particle species)*(3 helicity states)
	};
	if(E > ms[14]){ // Z0 boson contribution
		result += 3.; // Z0 boson (3 helicity states)
	};
	if(E > ms[4]){ // gluons contribution
		result += 16.; // gluons (8 particle species)*(2 helicity states)
	};
	return result;
};

double nb_spin_2(double E,double *ms){
	// This function computes the number of spin 2 dof.
	
	double result = grav*2.; // graviton (2 polarization states)
	if(ADD_DM && spin_DM == 2. && E > m_DM){
		result += dof_DM; // ADD
	};
	return result;
};

double nb_spin_05(double E,double *ms){
	// This function computes the number of spin 1/2 dof.
	
	double result = 6.; // neutrinos (3 particle species)*(1 helicity state)*(2 for antiparticles)
	if(ADD_DM && spin_DM == 0.5 && E > m_DM){
		result += dof_DM; // ADD
	};
	if(E > ms[0]){ // bottom quark contribution
		result += 12.; // bottom quark (2 helicity states)*(3 color states)*(2 for antiparticles)
	};
	if(E > ms[1]){ // charm quark contribution
		result += 12.; // charm quark (2 helicity states)*(3 color states)*(2 for antiparticles)
	};
	if(E > ms[2]){ // down quark contribution
		result += 12.; // down quark (2 helicity states)*(3 color states)*(2 for antiparticles)
	};
	if(E > ms[3]){ // electron contribution
		result += 4.; // electrons (2 helicity states)*(2 for antiparticles)
	};
	if(E > ms[6]){ // muon contribution
		result += 4.; // muon (2 helicity states)*(2 for antiparticles)
	};
	if(E > ms[9]){ // strange quark contribution
		result += 12.; // strange quark (2 helicity states)*(3 color states)*(2 for antiparticles)
	};
	if(E > ms[10]){ // tau contribution
		result += 4.; // tau (2 helicity states)*(2 for antiparticles)
	};
	if(E > ms[11]){ // top quark contribution
		result += 12.; // top quark (2 helicity states)*(3 color states)*(2 for antiparticles)
	};
	if(E > ms[12]){ // up quark contribution
		result += 12.; // up quark (2 helicity states)*(3 color states)*(2 for antiparticles)
	};
	return result;
};

double nb_spin_15(double E,double *ms){
	// This function computes the number of spin 3/2 dof.
	
	double result = 0.;
	if(ADD_DM && spin_DM == 1.5 && E > m_DM){
		result += dof_DM;
	};
	return result;
};

double sum_gammas_Kerr(double M,double a,double E,double ***gammas,int i,double *gamma_x,double ***fits,double *ms){
	// This function computes the sum over the dof of the Kerr emissivity.
	
	double x = 2.*E*M*G;
	double result = 0.;
	if(x <= 0.01){ // low-energy fit
		result += nb_spin_0(E,ms)*pow(10.,fits[0][i][0]*log10(2.*E*M*G) + fits[0][i][1])/(2.*pi); // spin 0 particles
		result += nb_spin_1(E,ms)*pow(10.,fits[1][i][0]*log10(2.*E*M*G) + fits[1][i][1])/(2.*pi); // spin 1 particles
		result += nb_spin_2(E,ms)*pow(10.,fits[2][i][0]*log10(2.*E*M*G) + fits[2][i][1])/(2.*pi); // spin 2 particles
		result += nb_spin_05(E,ms)*pow(10.,fits[3][i][0]*log10(2.*E*M*G) + fits[3][i][1])/(2.*pi); // spin 1/2 particles
		result += nb_spin_15(E,ms)*pow(10.,fits[4][i][0]*log10(2.*E*M*G) + fits[4][i][1])/(2.*pi); // spin 3/2 particles
	}
	else if(x >= 5.){ // high-energy fit
		result += nb_spin_0(E,ms)*pow(10.,fits[0][i][2]*(2.*M*E*G) + fits[0][i][3] + fits[0][i][4]*cos(fits[0][i][6]*2.*E*M*G) + fits[0][i][5]*sin(fits[0][i][6]*2.*E*M*G))/(2.*pi); // spin 0 particles
		result += nb_spin_1(E,ms)*pow(10.,fits[1][i][2]*(2.*M*E*G) + fits[1][i][3] + fits[1][i][4]*cos(fits[1][i][6]*2.*E*M*G) + fits[1][i][5]*sin(fits[1][i][6]*2.*E*M*G))/(2.*pi); // spin 1 particles
		result += nb_spin_2(E,ms)*pow(10.,fits[2][i][2]*(2.*M*E*G) + fits[2][i][3] + fits[2][i][4]*cos(fits[2][i][6]*2.*E*M*G) + fits[2][i][5]*sin(fits[2][i][6]*2.*E*M*G))/(2.*pi); // spin 2 particles
		result += nb_spin_05(E,ms)*pow(10.,fits[3][i][2]*(2.*M*E*G) + fits[3][i][3] + fits[3][i][4]*cos(fits[3][i][6]*2.*E*M*G) + fits[3][i][5]*sin(fits[3][i][6]*2.*E*M*G))/(2.*pi); // spin 1/2 particles
		result += nb_spin_15(E,ms)*pow(10.,fits[4][i][2]*(2.*M*E*G) + fits[4][i][3] + fits[4][i][4]*cos(fits[4][i][6]*2.*E*M*G) + fits[4][i][5]*sin(fits[4][i][6]*2.*E*M*G))/(2.*pi); // spin 3/2 particles
	}
	else{ // tabulated value
		int counter_x = 0;
		while(counter_x < nb_x - 1 && 2.*E*M*G >= gamma_x[counter_x]){
			counter_x++;
		};
		result += nb_spin_0(E,ms)*(gammas[0][i][counter_x-1] + (gammas[0][i][counter_x] - gammas[0][i][counter_x-1])/(gamma_x[counter_x] - gamma_x[counter_x-1])*(2.*E*M*G - gamma_x[counter_x-1]))/(2.*pi); // spin 0 particles
		result += nb_spin_1(E,ms)*(gammas[1][i][counter_x-1] + (gammas[1][i][counter_x] - gammas[1][i][counter_x-1])/(gamma_x[counter_x] - gamma_x[counter_x-1])*(2.*E*M*G - gamma_x[counter_x-1]))/(2.*pi); // spin 1 particles
		result += nb_spin_2(E,ms)*(gammas[2][i][counter_x-1] + (gammas[2][i][counter_x] - gammas[2][i][counter_x-1])/(gamma_x[counter_x] - gamma_x[counter_x-1])*(2.*E*M*G - gamma_x[counter_x-1]))/(2.*pi); // spin 2 particles
		result += nb_spin_05(E,ms)*(gammas[3][i][counter_x-1] + (gammas[3][i][counter_x] - gammas[3][i][counter_x-1])/(gamma_x[counter_x] - gamma_x[counter_x-1])*(2.*E*M*G - gamma_x[counter_x-1]))/(2.*pi); // spin 1/2 particles
		result += nb_spin_15(E,ms)*(gammas[4][i][counter_x-1] + (gammas[4][i][counter_x] - gammas[4][i][counter_x-1])/(gamma_x[counter_x] - gamma_x[counter_x-1])*(2.*E*M*G - gamma_x[counter_x-1]))/(2.*pi); // spin 3/2 particles
	};
	return result;
};

double sum_gammas_LQG(double M,double epsilon,double E,double ***gammas,int i,double *gamma_x,double ***fits,double *ms){
	// This function computes the sum over the dof of the Kerr emissivity.
	
	double x = 2.*E*M*G;
	double result = 0.;
	double temp = temp_LQG(M,epsilon);
	if(x <= 0.01){ // low-energy fit
		result += nb_spin_0(E,ms)*pow(10.,fits[0][i][0]*log10(2.*E*M*G) + fits[0][i][1])/(2.*pi); // spin 0 particles
		result += nb_spin_1(E,ms)*pow(10.,fits[1][i][0]*log10(2.*E*M*G) + fits[1][i][1])/(2.*pi); // spin 1 particles
		result += nb_spin_2(E,ms)*pow(10.,fits[2][i][0]*log10(2.*E*M*G) + fits[2][i][1])/(2.*pi); // spin 2 particles
		result += nb_spin_05(E,ms)*pow(10.,fits[3][i][0]*log10(2.*E*M*G) + fits[3][i][1])/(2.*pi); // spin 1/2 particles
	}
	else if(x >= 5.){ // high-energy fit -> this is the (tabulated) geometrical optics approximation
		result += nb_spin_0(E,ms)*(27./4.*pow(x,2.)/expbis(E/temp)*fits[0][i][2])/(2.*pi); // spin 0 particles
		result += nb_spin_1(E,ms)*(27./4.*pow(x,2.)/expbis(E/temp)*fits[1][i][2])/(2.*pi); // spin 1 particles
		result += nb_spin_2(E,ms)*(27./4.*pow(x,2.)/expbis(E/temp)*fits[2][i][2])/(2.*pi); // spin 2 particles
		result += nb_spin_05(E,ms)*(27./4.*pow(x,2.)/(exp(E/temp) + 1.)*fits[3][i][2])/(2.*pi); // spin 1/2 particles
	}
	else{ // tabulated value
		int counter_x = 0;
		while(counter_x < nb_x - 1 && 2.*E*M*G >= gamma_x[counter_x]){
			counter_x++;
		};
		result += nb_spin_0(E,ms)*(gammas[0][i][counter_x-1] + (gammas[0][i][counter_x] - gammas[0][i][counter_x-1])/(gamma_x[counter_x] - gamma_x[counter_x-1])*(2.*E*M*G - gamma_x[counter_x-1]))/(2.*pi); // spin 0 particles
		result += nb_spin_1(E,ms)*(gammas[1][i][counter_x-1] + (gammas[1][i][counter_x] - gammas[1][i][counter_x-1])/(gamma_x[counter_x] - gamma_x[counter_x-1])*(2.*E*M*G - gamma_x[counter_x-1]))/(2.*pi); // spin 1 particles
		result += nb_spin_2(E,ms)*(gammas[2][i][counter_x-1] + (gammas[2][i][counter_x] - gammas[2][i][counter_x-1])/(gamma_x[counter_x] - gamma_x[counter_x-1])*(2.*E*M*G - gamma_x[counter_x-1]))/(2.*pi); // spin 2 particles
		result += nb_spin_05(E,ms)*(gammas[3][i][counter_x-1] + (gammas[3][i][counter_x] - gammas[3][i][counter_x-1])/(gamma_x[counter_x] - gamma_x[counter_x-1])*(2.*E*M*G - gamma_x[counter_x-1]))/(2.*pi); // spin 1/2 particles
	};
	return result;
};

void compute_fM_Kerr(double **fM,double *masses,double ***gammas,double *gamma_a,double *gamma_x,double ***fits,double *ms){
	// This function computes the value of f(M,a*).
	
	double temp;
	double *energies;
	int nb_energies = 1000;
	energies = (double *)malloc(nb_energies*sizeof(double));
	for(int i = 0;i<nb_a;i++){
		printf("%i/%i ",i+1,nb_a);
		fflush(stdout);
		for(int j = 0;j<nb_masses;j++){
			temp = temp_Kerr(masses[j],gamma_a[i]);
			for(int k = 0;k<nb_energies;k++){
				energies[k] = pow(10.,log10(temp*1.e-5) + (log10(temp*1.e+5) - log10(temp*1.e-5))/(nb_energies - 1)*k);
			};
			for(int k = 0;k<nb_energies-1;k++){
				fM[j][i] += (energies[k+1] - energies[k])*energies[k]*pow(masses[j],2.)*sum_gammas_Kerr(masses[j],gamma_a[i],energies[k],gammas,i,gamma_x,fits,ms);
			};
		};
	};
	free(energies);
	return;
};

void compute_LQG(double **fM,double *masses,double ***gammas,double *gamma_epsilon,double *gamma_x,double ***fits,double *ms){
	// This function computes the value of f(M,epsilon).
	
	double temp;
	double *energies;
	int nb_energies = 1000;
	energies = (double *)malloc(nb_energies*sizeof(double));
	for(int i = 0;i<nb_epsilon;i++){
		printf("%i/%i ",i+1,nb_epsilon);
		fflush(stdout);
		for(int j = 0;j<nb_masses;j++){
			temp = temp_LQG(masses[j],gamma_epsilon[i]);
			for(int k = 0;k<nb_energies;k++){
				energies[k] = pow(10.,log10(temp*1.e-5) + (log10(temp*1.e+5) - log10(temp*1.e-5))/(nb_energies - 1)*k);
			};
			for(int k = 0;k<nb_energies-1;k++){
				fM[j][i] += (energies[k+1] - energies[k])*energies[k]*pow(masses[j],2.)*sum_gammas_LQG(masses[j],gamma_epsilon[i],energies[k],gammas,i,gamma_x,fits,ms);
			};
		};
	};
	free(energies);
	return;
};

void write_fM_Kerr(double **fM,double *masses,double *gamma_a,char *f_name){
	// This function writes the f(M,a*) tables.
	
	FILE *file = fopen(f_name,"w+"); // WARNING
	if(!file){
		printf("\n Could not open file %s\n",f_name);
		fflush(stdout);
		return;
	};
	rewind(file);
	fprintf(file,"%15s","mass/a");
	for(int i = 0;i<nb_a;i++){
		fprintf(file,"%15.5e",gamma_a[i]);
	};
	fprintf(file,"\n");
	for(int i = 0;i<nb_masses;i++){
		fprintf(file,"%15.5e",masses[i]);
		for(int j = 0;j<nb_a;j++){
			fprintf(file,"%15.5e",fM[i][j]);
		};
		fprintf(file,"\n");
	};
	fclose(file);
	return;
};

void write_LQG(double **fM,double *masses,double *gamma_epsilon,char *f_name){
	// This function writes the f(M,epsilon) tables.
	
	FILE *file = fopen(f_name,"w+"); // WARNING
	if(!file){
		printf("\n Could not open file %s\n",f_name);
		fflush(stdout);
		return;
	};
	rewind(file);
	fprintf(file,"%15s","mass/epsilon");
	for(int i = 0;i<nb_epsilon;i++){
		fprintf(file,"%15.5e",gamma_epsilon[i]);
	};
	fprintf(file,"\n");
	for(int i = 0;i<nb_masses;i++){
		fprintf(file,"%15.5e",masses[i]);
		for(int j = 0;j<nb_epsilon;j++){
			fprintf(file,"%15.5e",fM[i][j]);
		};
		fprintf(file,"\n");
	};
	fclose(file);
	return;
};

void read_gamma_tables_Kerr_2(double ***gammas,double *gamma_a,double *gamma_x){
	// This function reads the m*greybody factors tables.
	
	char **gamma_names = (char **)malloc(nb_particle_spins*sizeof(char *));
	for(int i = 0;i<nb_particle_spins;i++){
		gamma_names[i] = (char *)malloc(64*sizeof(char));
	};
	sprintf(gamma_names[0],"./spin_0_gM.txt");
	sprintf(gamma_names[1],"./spin_1_gM.txt");
	sprintf(gamma_names[2],"./spin_2_gM.txt");
	sprintf(gamma_names[3],"./spin_0.5_gM.txt");
	sprintf(gamma_names[4],"./spin_1.5_gM.txt");
	
	char dummy[32];
	FILE *table;
	for(int i = 0;i<nb_particle_spins;i++){ // WARNING the spin 3/2 table is dummy, it was only computed for the case a* = 0
		table = fopen(gamma_names[i],"r");
		rewind(table);
		if(!table){
			printf("\n[read_simple_tables] : ERROR COULD NOT FIND TABLE %s !\n",gamma_names[i]);
			fflush(stdout);
			return;
		};
		fscanf(table,"%s",dummy);
		for(int j = 0;j<nb_x;j++){
			fscanf(table,"%lf",&(gamma_x[j]));
		};
		for(int j = 0;j<nb_a;j++){
			fscanf(table,"%lf",&(gamma_a[j]));
			for(int k = 0;k<nb_x;k++){
				fscanf(table,"%lf",&(gammas[i][j][k]));
			};
		};
		fclose(table);
	};
	for(int i = 0;i<nb_particle_spins;i++){
		free(gamma_names[i]);
	};
	return;
};

void read_gM_fits(double ***fits){
	// This function reads the tables fitting m*asymptotic Kerr emissivities.
	
	char **table_names;
	table_names = (char **)malloc(nb_particle_spins*sizeof(char *));
	for(int i = 0;i<nb_particle_spins;i++){
		table_names[i] = (char *)malloc(64*sizeof(char));
	};
	sprintf(table_names[0],"./spin_0_fits_gM.txt");
	sprintf(table_names[1],"./spin_1_fits_gM.txt");
	sprintf(table_names[2],"./spin_2_fits_gM.txt");
	sprintf(table_names[3],"./spin_0.5_fits_gM.txt");
	sprintf(table_names[4],"./spin_1.5_fits_gM.txt");
	
	FILE *file;
	char dummy[32];
	for(int i = 0;i<nb_particle_spins;i++){ // WARNING the spin 3/2 table is dummy as it was only computed for the case a* = 0
		file = fopen(table_names[i],"r");
		if(file == 0){
			printf("[read_asymp_fits] : ERROR COULD NOT OPEN FILE %s !\n",table_names[i]);
			fflush(stdout);
			return;
		};
		rewind(file);
		for(int j = 0;j<8;j++){
			fscanf(file,"%s",dummy);
		};
		for(int j = 0;j<nb_a;j++){
			fscanf(file,"%s",dummy);
			for(int k = 0;k<7;k++){
				fscanf(file,"%lf",&(fits[i][j][k]));
			};
		};
		fclose(file);
	};
	for(int i = 0;i<nb_particle_spins;i++){
		free(table_names[i]);
	};
	return;
};

void compute_gM(double **gM,double *masses,double ***gammas,double *gamma_a,double *gamma_x,double ***fits,double *ms){
	// This function computes the value of g(M,a).
	
	double temp;
	double *energies;
	int nb_energies = 1000;
	energies = (double *)malloc(nb_energies*sizeof(double));
	for(int i = 0;i<nb_a;i++){
		printf("%i/%i ",i+1,nb_a);
		fflush(stdout);
		for(int j = 0;j<nb_masses;j++){
			temp = temp_Kerr(masses[j],gamma_a[i]);
			for(int k = 0;k<nb_energies;k++){
				energies[k] = pow(10.,log10(temp*1.e-5) + (log10(temp*1.e+5) - log10(temp*1.e-5))/(nb_energies - 1)*k);
			};
			for(int k = 0;k<nb_energies-1;k++){
				if(gamma_a[i] != 0.){
					gM[j][i] += (energies[k+1] - energies[k])*masses[j]/gamma_a[i]*sum_gammas_Kerr(masses[j],gamma_a[i],energies[k],gammas,i,gamma_x,fits,ms);
				};
			};
		};
	};
	for(int j = 0;j<nb_masses;j++){
		gM[j][0] = gM[j][1] - (gM[j][2] - gM[j][1])/(gamma_a[2] - gamma_a[1])*(gamma_a[1] - gamma_a[0]);
	};
	free(energies);
	return;
};

void write_gM(double **gM,double *masses,double *gamma_a,char *g_name){
	// This function writes the g(M,a) tables.
	
	FILE *file = fopen(g_name,"w+"); // WARNING
	if(!file){
		printf("\n Could not open file %s\n",g_name);
		fflush(stdout);
		return;
	};
	rewind(file);
	fprintf(file,"%15s","mass/a");
	for(int i = 0;i<nb_a;i++){
		fprintf(file,"%15.5e",gamma_a[i]);
	};
	fprintf(file,"\n");
	for(int i = 0;i<nb_masses;i++){
		fprintf(file,"%15.5e",masses[i]);
		for(int j = 0;j<nb_a;j++){
			fprintf(file,"%15.5e",gM[i][j]/G);
		};
		fprintf(file,"\n");
	};
	fclose(file);
	return;
};

/* MAIN FUNCTION */

int main(int argc, char **argv){
	
	int type;
	
	if(argc < 2){
		printf("\n [main] : ERROR this program needs at least one parameter, the type of tables");
		printf("\n END OF EXECUTION\n");
		fflush(stdout);
		exit(0);
	}
	else{
		sscanf(argv[1],"%i",&(type));
	};
	
	if(type == 1 && ADD_DM && spin_DM == 1.5){
		printf("\n [main] : ERROR the greybody factors for the spin 3/2 particle in polymerized metric is not available");
		printf("\n END OF EXECUTION\n");
		fflush(stdout);
		exit(0);
	};
		
	char *f_name = (char *)malloc(64*sizeof(char));
	char *g_name = (char *)malloc(64*sizeof(char));
	
	printf("\n Define the output file names:");
	printf("\n\t f file: ");
	fflush(stdout);
	scanf("%s",f_name);
	if(type == 0){
		printf("\n\t g file: ");
		fflush(stdout);
		scanf("%s",g_name);
	};
	
	double ms[nb_particles+1]; // masses in GeV
	ms[0] = 4.18; // bottom quark
	ms[1] = 1.27; // charm quark
	ms[2] = 4.7e-3; // down quark
	ms[3] = 0.5109989461e-3; // electron
	ms[4] = 200.e-3; // gluons, effective mass
	ms[5] = 125.03; // Higgs boson
	ms[6] = 105.6583745e-3; // muon
	ms[7] = 0.; // neutrinos
	ms[8] = 0.; // photon
	ms[9] = 96e-3; // strange quark
	ms[10] = 1.77686; // tau
	ms[11] = 173.21; // top quark
	ms[12] = 2.2e-3; // up quark
	ms[13] = 80.403; // W boson
	ms[14] = 91.1876; // Z boson
	ms[15] = 0.; // graviton
	
	double *gamma_a;
	double *gamma_epsilon;
	double *gamma_x;
	double ***gammas;
	double ***fits;
	double **fM;
	double *masses;
	masses = (double *)malloc(nb_masses*sizeof(double));
	for(int i = 0;i<nb_masses;i++){
		masses[i] = pow(10.,log10(Mp) + (log10(Mmax) - log10(Mp))/(nb_masses - 1)*i);
	};
	
	switch(type){
		
		case 0:{ // Kerr tables
			nb_particle_spins = 5; // spins 0, 1, 2, 1/2 and 3/2
			
			gamma_a = (double *)malloc(nb_a*sizeof(double));
			gamma_x = (double *)malloc(nb_x*sizeof(double));
			gammas = (double ***)malloc(nb_particle_spins*sizeof(double **));
			
			printf(" Reading gamma tables...");
			fflush(stdout);
			
			for(int i = 0;i<nb_particle_spins;i++){
				gammas[i] = (double **)malloc(nb_a*sizeof(double *));
				for(int j = 0;j<nb_a;j++){
					gammas[i][j] = (double *)malloc(nb_x*sizeof(double));
				};
			};
			
			read_gamma_tables_Kerr(gammas,gamma_a,gamma_x);
			
			printf("\t DONE");
			printf("\n Reading asymptotical fits...");
			fflush(stdout);
			
			fits = (double ***)malloc(nb_particle_spins*sizeof(double **));
			for(int i = 0;i<nb_particle_spins;i++){
				fits[i] = (double **)malloc(nb_a*sizeof(double *));
				for(int j = 0;j<nb_a;j++){
					fits[i][j] = (double *)malloc(7*sizeof(double));
				};
			};
			
			read_fM_fits_Kerr(fits);
			
			printf("\t DONE");
			printf("\n Computing f(M,a*) table...");
			fflush(stdout);
			
			fM = (double **)malloc(nb_masses*sizeof(double *));
			for(int i = 0;i<nb_masses;i++){
				fM[i] = (double *)malloc(nb_a*sizeof(double));
				for(int j = 0;j<nb_a;j++){
					fM[i][j] = 0.;
				};
			};
			
			compute_fM_Kerr(fM,masses,gammas,gamma_a,gamma_x,fits,ms);
	
			printf("\t DONE");
			printf("\n Writing f(M,a*) table...");
			
			write_fM_Kerr(fM,masses,gamma_a,f_name);
			
			for(int i = 0;i<nb_masses;i++){
				free(fM[i]);
			};
			
			printf("\t DONE");
			printf("\n Reading m*gamma tables...");
			fflush(stdout);
			
			read_gamma_tables_Kerr_2(gammas,gamma_a,gamma_x);
			
			printf("\t DONE");
			printf("\n Reading m*fits...");
			fflush(stdout);
			
			read_gM_fits(fits);
			
			printf("\t\t DONE");
			printf("\n Computing g(M,a*) table...");
			fflush(stdout);
			
			double **gM;
			gM = (double **)malloc(nb_masses*sizeof(double *));
			for(int i = 0;i<nb_masses;i++){
				gM[i] = (double *)malloc(nb_a*sizeof(double));
				for(int j = 0;j<nb_a;j++){
					gM[i][j] = 0.;
				};
			};
			
			compute_gM(gM,masses,gammas,gamma_a,gamma_x,fits,ms);
			
			printf("\t DONE");
			printf("\n Writing g(M,a) table...");
			fflush(stdout);
			
			write_gM(gM,masses,gamma_a,g_name);
			
			printf("\t DONE\n");
			fflush(stdout);
			
			break;
		};
		
		case 1:{ // polymerized tables
			nb_particle_spins = 4; // spins 0, 1, 2 and 1/2
			
			gamma_epsilon = (double *)malloc(nb_epsilon*sizeof(double));
			gamma_x = (double *)malloc(nb_x*sizeof(double));
			gammas = (double ***)malloc(nb_particle_spins*sizeof(double **));
			
			printf(" Reading gamma tables...");
			fflush(stdout);
			
			for(int i = 0;i<nb_particle_spins;i++){
				gammas[i] = (double **)malloc(nb_epsilon*sizeof(double *));
				for(int j = 0;j<nb_epsilon;j++){
					gammas[i][j] = (double *)malloc(nb_x*sizeof(double));
				};
			};
			
			read_gamma_tables_LQG(gammas,gamma_epsilon,gamma_x);
			
			printf("\t DONE");
			printf("\n Reading asymptotical fits...");
			fflush(stdout);
			
			fits = (double ***)malloc(nb_particle_spins*sizeof(double **));
			for(int i = 0;i<nb_particle_spins;i++){
				fits[i] = (double **)malloc(nb_epsilon*sizeof(double *));
				for(int j = 0;j<nb_epsilon;j++){
					fits[i][j] = (double *)malloc(3*sizeof(double));
				};
			};
			
			read_fits_LQG(fits);
			
			printf("\t DONE");
			printf("\n Computing f(M,epsilon) table...");
			fflush(stdout);
			
			fM = (double **)malloc(nb_masses*sizeof(double *));
			for(int i = 0;i<nb_masses;i++){
				fM[i] = (double *)malloc(nb_epsilon*sizeof(double));
				for(int j = 0;j<nb_epsilon;j++){
					fM[i][j] = 0.;
				};
			};
			
			compute_LQG(fM,masses,gammas,gamma_epsilon,gamma_x,fits,ms);
			
			printf("\t DONE");
			printf("\n Writing f(M,epsilon) table...");
			
			write_LQG(fM,masses,gamma_epsilon,f_name);
			
			for(int i = 0;i<nb_masses;i++){
				free(fM[i]);
			};
			
			printf("\t DONE");
			
			break;
		};
		
		default:{
			printf("\n [main] : ERROR wrong type!\n END OF EXECUTION");
			fflush(stdout);
			exit(0);
		};
	};
	
	printf("\n END OF EXECUTION\n");
	fflush(stdout);
	
	return 1;
};










