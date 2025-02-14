// Program that stacks the time-dependent spectra of BlackHawk results
// taking the redshift into account
// Authors: Alexandre Arbey, alexandre.arbey@ens-lyon.fr & Jérémy Auffinger, jeremy.auffinger.blackhawk@protonmail.com
// Last modification: 07 April 2024

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
//#include <direct.h>
#include <dirent.h>

int nb_params = 5;
int nb_fin_ener = 500;
double Emin = 1e-5; // in GeV
double Emax = 1e+19; // in GeV
double t_eq = 5e+4*365.*24.*3600.; // in s
double t_today = 13.8e+9*365.*24.*3600.; // in s
double t_CMB = 3.8e+5*365.*24.*3600.; // in s
int domination = 0; // 0: radiation, 1: matter, 2: alternate with default cosmology
int CMB = 0; // determines whether the stacking starts at CMB, 0: no it starts before, 1: yes
int today = 0; // determines whether the redshifting stops at today's time (1) or at BH total evaporation (0)

double redshift(double t,double t_evap,int today){
	switch(domination){
		case 0:{
			return pow(t_evap/t,1./2.) - 1.;
			break;
		}
		case 1:{
			return pow(t_evap/t,2./3.) - 1.;
			break;
		}
		case 2:{
			if(t < t_eq){ // before M-R equality: radiation domination
				if(today == 1){
					return pow(pow(t_today,4./3.)/(t*pow(t_eq,1./3.)),1./2.) - 1.;
				}
				else{
					return pow(pow(t_evap,4./3.)/(t*pow(t_eq,1./3.)),1./2.) - 1.;
				}
			}
			else{ // after M-R equality: matter domination
				if(today == 1){
					return pow(t_today/t,2./3.) - 1.;
				}
				else{
					return pow(t_evap/t,2./3.) - 1.;
				}
			}
			break;
		}
		default:{
			printf("[redshift] : ERROR wrong redshift choice!\n");
			printf("[main] : end of execution");
			fflush(stdout);
			exit(0);
		}
	}
}

int main(int argc, char **argv){
	
	printf("[main] : start of execution\n");
	fflush(stdout);
	
	char *path = (char *)malloc(256*sizeof(char));
	int nb_init_ener;
	
	if(argc < nb_params+1){
		printf("[main] : ERROR this program needs at least %i parameters!\n",nb_params);
		printf("[main] : end of execution");
		fflush(stdout);
		exit(0);
	}
	else{
		sscanf(argv[1],"%s",path);
		sscanf(argv[2],"%i",&nb_init_ener);
		sscanf(argv[3],"%i",&domination);
		sscanf(argv[4],"%i",&CMB);
		sscanf(argv[5],"%i",&today);
	}
	
	printf("[main] : recovering evolution data\n");
	fflush(stdout);
	char *name_life_evolutions = (char *)malloc(256*sizeof(char));
	sprintf(name_life_evolutions,"%slife_evolutions.txt",path);
	FILE *life_evolutions = fopen(name_life_evolutions,"r+");
	int nb_times;
	if(!life_evolutions){
		printf("[main] : ERROR file not found!\n");
		printf("\t%s\n",name_life_evolutions);
		printf("[main] : end of execution");
		fflush(stdout);
		exit(0);
	}
	rewind(life_evolutions);
	char *dumb = (char *)malloc(32*sizeof(char));
	for(int i=0;i<16;i++){
		fscanf(life_evolutions,"%s",dumb);
	}
	fscanf(life_evolutions,"%i",&nb_times);
	for(int i=0;i<3;i++){
		fscanf(life_evolutions,"%s",dumb);
	}
	double *times = (double *)malloc(nb_times*sizeof(double));
	double *masses = (double *)malloc(nb_times*sizeof(double));
	double *spins = (double *)malloc(nb_times*sizeof(double));
	double *dts = (double *)malloc(nb_times*sizeof(double));
	for(int i = 0;i<nb_times;i++){
		fscanf(life_evolutions,"%lf",&(times[i]));
		fscanf(life_evolutions,"%lf",&(masses[i]));
		fscanf(life_evolutions,"%lf",&(spins[i]));
	}
	fclose(life_evolutions);
	char *name_dts = (char *)malloc(256*sizeof(char));
	sprintf(name_dts,"%sdts.txt",path);
	FILE *dts_file = fopen(name_dts,"r+");
	if(!dts_file){
		printf("[main] : ERROR file not found!\n");
		printf("\t%s\n",name_dts);
		printf("[main] : end of execution");
		fflush(stdout);
		exit(0);
	}
	for(int i = 0;i<12;i++){
		fscanf(dts_file,"%s",dumb);
	}
	for(int i = 0;i<nb_times;i++){
		fscanf(dts_file,"%s",dumb);
		fscanf(dts_file,"%lf",&(dts[i]));
	}
	fclose(dts_file);
	
	printf("[main] : recovering spectrum data\n");
	fflush(stdout);
	char *name_spectrum = (char *)malloc(256*sizeof(char));
	sprintf(name_spectrum,"%sphoton_secondary_spectrum.txt",path); // ATTENTION
	FILE *spectrum_file = fopen(name_spectrum,"r+");
	if(!life_evolutions){
		printf("[main] : ERROR file not found!\n");
		printf("\t%s\n",name_spectrum);
		printf("[main] : end of execution");
		fflush(stdout);
		exit(0);
	}
	rewind(spectrum_file);
	double *init_ener = (double *)malloc(nb_init_ener*sizeof(double));
	double **init_spectrum = (double **)malloc(nb_times*sizeof(double *));
	for(int i = 0;i<nb_times;i++){
		init_spectrum[i] = (double *)malloc(nb_init_ener*sizeof(double));
	}
	for(int i = 0;i<9;i++){
		fscanf(spectrum_file,"%s",dumb);
	}
	for(int i = 0;i<nb_init_ener;i++){
		fscanf(spectrum_file,"%lf",&(init_ener[i]));
	}
	for(int i = 0;i<nb_times;i++){
		fscanf(spectrum_file,"%s",dumb);
		for(int j = 0;j<nb_init_ener;j++){
			fscanf(spectrum_file,"%lf",&(init_spectrum[i][j]));
		}
	}
	fclose(spectrum_file);
	
	printf("[main] : computing stacked spectrum\n");
	fflush(stdout);
	double *final_ener = (double *)malloc(nb_fin_ener*sizeof(double));
	for(int i = 0;i<nb_fin_ener;i++){
		final_ener[i] = pow(10.,log10(Emin) + (log10(Emax) - log10(Emin))/(nb_fin_ener - 1)*i);
	}
	
	double t_evap = times[nb_times-1];
	
	double **interm_spectrum = (double **)malloc(nb_times*sizeof(double *));
	for(int i = 0;i<nb_times;i++){
		interm_spectrum[i] = (double *)malloc(nb_fin_ener*sizeof(double));
		for(int j = 0;j<nb_fin_ener;j++){
			interm_spectrum[i][j] = 0.;
		}
	}
	int counter;
	int too_low;
	double z;
	for(int i = 0;i<nb_times;i++){
		if((CMB == 0 || (CMB == 1 && times[i] > t_CMB)) && (today == 0 || (today == 1 && times[i] < t_today))){
			z = redshift(times[i],t_evap,today);
			for(int j = 0;j<nb_fin_ener;j++){
				counter = 0;
				too_low = 1;
				while(counter < nb_init_ener && init_ener[counter]/(1.+z) <= final_ener[j]){
					too_low = 0;
					counter++;
				}
				if(too_low || counter == nb_init_ener){
					interm_spectrum[i][j] = 0.;
				}
				else{
					interm_spectrum[i][j] = fabs(init_spectrum[i][counter-1] + (init_spectrum[i][counter] - init_spectrum[i][counter-1])/(init_ener[counter] - init_ener[counter-1])*(final_ener[j] - init_ener[counter-1]/(1.+z))); // HERE
				}				
				printf("\n%.5e",init_spectrum[i][counter]); // HERE
				fflush(stdout);
			}
		}
	}
	
	double *final_spectrum = (double *)malloc(nb_fin_ener*sizeof(double));
	for(int i = 0;i<nb_fin_ener;i++){
		final_spectrum[i] = 0.;
	}
	for(int i = 0;i<nb_times-1;i++){
		z = redshift(times[i],t_evap,today);
		for(int j = 0;j<nb_fin_ener;j++){
			final_spectrum[j] += (1.+z)*interm_spectrum[i][j]*dts[i];
		}
	}
	
	FILE *results;
	char *results_file = (char *)malloc(64*sizeof(char));
	if(domination == 0){
		sprintf(results_file,"results_RD.txt");
	}
	else if(domination == 1){
		sprintf(results_file,"results_MD.txt");
	}
	else{
		sprintf(results_file,"results_alternate.txt");
	}
	results = fopen(results_file,"w+");
	if(!results){
		printf("[main] : ERROR could not open file %s\n",results_file);
		printf("[main] : end of execution");
		fflush(stdout);
		exit(0);
	}
	rewind(results);
	fprintf(results,"%15s%15s","energy","spectrum");
	for(int i = 0;i<nb_fin_ener;i++){
		fprintf(results,"\n%15.5e%15.5e",final_ener[i],final_spectrum[i]);
	}
	fclose(results);
	
	printf("[main] : end of execution");
	fflush(stdout);
	
	return 1;
}
