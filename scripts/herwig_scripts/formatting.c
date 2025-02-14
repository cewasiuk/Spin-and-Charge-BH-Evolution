// Script to format the HERWIG hadronization tables

#include <stdbool.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <dirent.h>

int nb_initial_particles = 14;
int nb_final_particles = 11;
int nb_final_energies = 100;
double Emin_init = 25.0;
double Emax_init = 100000.0;
double Emin_fin = 1.0e-6;
double Emax_fin = 100000.0;
int nb_initial_energies = 100;

void read_yoda(char *path,double **result){
	FILE *yoda = fopen(path,"r");
	if(!yoda){ // in this case the file LEP.yoda doesn't exist, so the computation was not made, so it seems the energy was too low or smth so we return all zeros
		return;
	};
	rewind(yoda);
	char dummy[256];
	for(int i = 0;i<97;i++){
		fgets(dummy,256,yoda);
		// printf("%s",dummy);
	};
	double dumb;
	for(int i = 0;i<nb_final_particles;i++){
		for(int j = 0;j<16;j++){
			fgets(dummy,256,yoda);
			// printf("%s",dummy);
		};
		for(int j = 0;j<nb_final_energies;j++){
			for(int k = 0;k<6;k++){
				fscanf(yoda,"%lf",&dumb);
			};
			fscanf(yoda,"%lf",&result[i][j]);
			// printf("%f\n",result[i][j]);
		};
		fgets(dummy,256,yoda);
	};
	for(int i = 0;i<7;i++){
		fgets(dummy,256,yoda);
		// printf("%s",dummy);
	};
	int nb_events;
	fscanf(yoda,"%lf",&dumb);
	fscanf(yoda,"%lf",&dumb);
	fscanf(yoda,"%i",&nb_events);
	for(int i = 0;i<nb_final_particles;i++){
		for(int j = 0;j<nb_final_energies;j++){
			if(nb_events != 0){
				result[i][j] = result[i][j]/nb_events;
			}
			else{
				result[i][j] = 0.;
			};
		};
	};
	fclose(yoda);
	return;
};

void write_tables(double ****tables,double *initial_energies,double *final_energies){
	char **table_names;
	table_names = (char **)malloc(nb_final_particles*sizeof(char *));
	for(int i = 0;i<nb_final_particles;i++){
		table_names[i] = (char *)malloc(20*sizeof(char));
	};
	sprintf(table_names[0],"electron.txt");
	sprintf(table_names[1],"nue.txt");
	sprintf(table_names[2],"muon.txt");
	sprintf(table_names[3],"K0L.txt");
	sprintf(table_names[4],"numu.txt");
	sprintf(table_names[5],"nutau.txt");
	sprintf(table_names[6],"pipm.txt");
	sprintf(table_names[7],"neutron.txt");
	sprintf(table_names[8],"photon.txt");
	sprintf(table_names[9],"proton.txt");
	sprintf(table_names[10],"Kpm.txt");
	
	for(int i = 0;i<nb_final_particles;i++){
		FILE *out = fopen(table_names[i],"w");
		fprintf(out,"%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n","initial_energy","final_energy","photon","gluons","higgs","Wpm","Z0","electron","muon","tau","up","down","charm","strange","top","bottom");
		for(int j = 0;j<nb_initial_energies;j++){
			for(int k = 0;k<nb_final_energies;k++){
				fprintf(out,"%15.5e%15.5e",initial_energies[j],final_energies[k]);
				for(int l = 0;l<nb_initial_particles;l++){
					if(l==2){ // in this case it is e+ e- -> Higgs, no need to divide by two
						fprintf(out,"%15.5e",tables[i][j][k][l]);
					}
					else{
						fprintf(out,"%15.5e",tables[i][j][k][l]/2.);
					};
				};
				fprintf(out,"\n");
			};
		};
		fclose(out);
	};
};

int main(){
	double *initial_energies = (double *)malloc(nb_initial_energies*sizeof(double));
	for(int i = 0;i<nb_initial_energies;i++){
		initial_energies[i] = pow(10.0,log10(Emin_init) + (log10(Emax_init) - log10(Emin_init))/(nb_initial_energies-1)*i);
	};
	double *final_energies = (double *)malloc(nb_final_energies*sizeof(double));
	for(int i = 0;i<nb_final_energies;i++){
		final_energies[i] = pow(10.0,log10(Emin_fin) + (log10(Emax_fin) - log10(Emin_fin))/(nb_final_energies-1)*i);
	};
	double ****result = (double ****)malloc(nb_initial_particles*sizeof(double ***));
	for(int i = 0;i<nb_initial_particles;i++){
		result[i] = (double ***)malloc(nb_initial_energies*sizeof(double **));
		for(int j = 0;j<nb_initial_energies;j++){
			result[i][j] = (double **)malloc(nb_final_particles*sizeof(double *));
			for(int k = 0;k<nb_final_particles;k++){
				result[i][j][k] = (double *)malloc(nb_final_energies*sizeof(double));
			};
		};
	};
	
	char **folder_names = (char **)malloc(nb_initial_particles*sizeof(char *));
	for(int i = 0;i<nb_initial_particles;i++){
		folder_names[i] = (char *)malloc(20*sizeof(char));
	};
	sprintf(folder_names[0],"photon");
	sprintf(folder_names[1],"gluons");
	sprintf(folder_names[2],"higgs");
	sprintf(folder_names[3],"Wpm");
	sprintf(folder_names[4],"Z0");
	sprintf(folder_names[5],"electron");
	sprintf(folder_names[6],"muon");
	sprintf(folder_names[7],"tau");
	sprintf(folder_names[8],"up");
	sprintf(folder_names[9],"down");
	sprintf(folder_names[10],"charm");
	sprintf(folder_names[11],"strange");
	sprintf(folder_names[12],"top");
	sprintf(folder_names[13],"bottom");

	char result_name[100];
	
	// Reading the computed tables
	
	printf("[main] : Reading computed tables...\n");

	for(int i = 0;i<nb_initial_particles;i++){
		printf("%i :",i);
		for(int j = 0;j<nb_initial_energies;j++){
			printf("%i\t",j);
			sprintf(result_name,"./%s/%.5eGeV/LEP.yoda",folder_names[i],initial_energies[j]);
			read_yoda(result_name,result[i][j]);
		};
		printf("\n");
	};

	printf("[main] : Reading complete\n");

	// Format conversion

	printf("[main] : Formating tables...\n");

	double ****tables = (double ****)malloc(nb_final_particles*sizeof(double ***));
	for(int i = 0;i<nb_final_particles;i++){
		tables[i] = (double ***)malloc(nb_initial_energies*sizeof(double **));
		for(int j = 0;j<nb_initial_energies;j++){
			tables[i][j] = (double **)malloc(nb_final_energies*sizeof(double *));
			for(int k = 0;k<nb_final_energies;k++){
				tables[i][j][k] = (double *)malloc(nb_initial_particles*sizeof(double));
				for(int l = 0;l<nb_initial_particles;l++){
					tables[i][j][k][l] = result[l][j][i][k]; // format conversion
				};
			};
		};
	};
	
	printf("[main] : Formating complete\n");

	// Writing table files

	printf("[main] : Writing into files 'particles.txt'...\n");

	write_tables(tables,initial_energies,final_energies);

	printf("[main] : Writing complete\n");
	
	return 1;
};
