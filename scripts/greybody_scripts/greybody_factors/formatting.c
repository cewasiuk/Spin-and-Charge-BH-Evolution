// Script to format simple greybody tables

#include <stdbool.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <dirent.h>

int nb_a = 50;
int nb_x = 200;
double xmin = 0.01;
double xmax = 5.;

void read_tables(char ***table_names,double ***gammas,double ***mgammas){
	FILE *table;
	for(int i = 0;i<4;i++){
		table = fopen(table_names[i][0],"r"); //  a =/= 0 tables
		if(!table){
			printf("[read_table] : Unable to read table %s\n",table_names[i][0]);
			fflush(stdout);
		};
		rewind(table);
		for(int j = 0;j<nb_a;j++){
			for(int k = 0;k<nb_x;k++){
				fscanf(table,"%lf",&(gammas[i][j][k]));
			};
		};
		fclose(table);
	};
	for(int i = 0;i<4;i++){
		table = fopen(table_names[i][1],"r"); //  a =/= 0 tables
		if(!table){
			printf("[read_table] : Unable to read table %s\n",table_names[i][1]);
			fflush(stdout);
		};
		rewind(table);
		for(int j = 0;j<nb_a;j++){
			for(int k = 0;k<nb_x;k++){
				fscanf(table,"%lf",&(mgammas[i][j][k]));
			};
		};
		fclose(table);
	};
	return;
};

void write_tables(char ***new_tables,double ***gammas,double ***mgammas,double *a,double *x){
	FILE *table;
	
	for(int i = 0;i<4;i++){
		table = fopen(new_tables[i][0],"w+");
		if(!table){
			printf("[write_table] : Unable to open file %s\n",new_tables[i][0]);
			fflush(stdout);
			return;
		};
		rewind(table);
		fprintf(table,"%15s","a/x");
		for(int j = 0;j<nb_x;j++){
			fprintf(table,"%15.5e",x[j]);
		};
		fprintf(table,"\n");
		for(int j = 0;j<nb_a;j++){
			fprintf(table,"%15.5e",a[j]);
			for(int k = 0;k<nb_x;k++){
				fprintf(table,"%15.5e",gammas[i][j][k]);
			};
			fprintf(table,"\n");
		};
		fclose(table);
	};
	for(int i = 0;i<4;i++){
		table = fopen(new_tables[i][1],"w+");
		if(!table){
			printf("[write_table] : Unable to open file %s\n",new_tables[i][1]);
			fflush(stdout);
			return;
		};
		rewind(table);
		fprintf(table,"%15s","a/x");
		for(int j = 0;j<nb_x;j++){
			fprintf(table,"%15.5e",x[j]);
		};
		fprintf(table,"\n");
		for(int j = 0;j<nb_a;j++){
			fprintf(table,"%15.5e",a[j]);
			for(int k = 0;k<nb_x;k++){
				fprintf(table,"%15.5e",mgammas[i][j][k]);
			};
			fprintf(table,"\n");
		};
		fclose(table);
	};
	
	free(a);
	free(x);
	
	return;
};

void read_fits(char ***fits_tables,double ****fits){
	FILE *table;
	for(int i = 0;i<4;i++){
		table = fopen(fits_tables[i][0],"r"); // fM fits
		if(!table){
			printf("[read_fits] : Unable to read table %s\n",fits_tables[i][0]);
			fflush(stdout);
			return;
		};
		rewind(table);
		for(int j = 0;j<nb_a;j++){
			for(int k = 0;k<7;k++){
				fscanf(table,"%lf",&(fits[i][0][j][k]));
			};
		};
		fclose(table);
	};
	for(int i = 0;i<4;i++){
		table = fopen(fits_tables[i][1],"r"); // gM fits
		if(!table){
			printf("[read_fits] : Unable to read table %s\n",fits_tables[i][1]);
			fflush(stdout);
			return;
		};
		rewind(table);
		for(int j = 0;j<nb_a;j++){
			for(int k = 0;k<7;k++){
				fscanf(table,"%lf",&(fits[i][1][j][k]));
			};
		};
		fclose(table);
	};
	return;
};

void write_fits(char ***new_fits,double ****fits,double *a){
	FILE *table;
	for(int i = 0;i<4;i++){
		table = fopen(new_fits[i][0],"w+");
		if(!table){
			printf("[write_fits] : Unable to open file %s\n",new_fits[i][0]);
			fflush(stdout);
			return;
		};
		rewind(table);
		fprintf(table,"%15s%15s%15s%15s%15s%15s%15s%15s\n","a/fits","a1","b1","a2","b2","c2","d2","e2");
		for(int j = 0;j<nb_a;j++){
			fprintf(table,"%15.5e",a[j]);
			for(int k = 0;k<7;k++){
				fprintf(table,"%15.5e",fits[i][0][j][k]);
			};
			fprintf(table,"\n");
		};
		fclose(table);
	};
	for(int i = 0;i<4;i++){
		table = fopen(new_fits[i][1],"w+");
		if(!table){
			printf("[write_fits] : Unable to open file %s\n",new_fits[i][1]);
			fflush(stdout);
			return;
		};
		rewind(table);
		fprintf(table,"%15s%15s%15s%15s%15s%15s%15s%15s\n","a/fits","a1","b1","a2","b2","c2","d2","e2");
		for(int j = 0;j<nb_a;j++){
			fprintf(table,"%15.5e",a[j]);
			for(int k = 0;k<7;k++){
				if(i == 0 && (k == 0 || k == 1) && j < 8){
					fprintf(table,"%15.5e",fits[i][1][8][k]);
				}
				else{
					fprintf(table,"%15.5e",fits[i][1][j][k]);
				};
			};
			fprintf(table,"\n");
		};
		fclose(table);
	};
	return;
};

int main(){
	
	char ***table_names;
	table_names = (char ***)malloc(4*sizeof(char **));
	for(int i = 0;i<4;i++){
		table_names[i] = (char **)malloc(2*sizeof(char *));
		for(int j = 0;j<2;j++){
			table_names[i][j] = (char *)malloc(32*sizeof(char));
		};
	};
	
	sprintf(table_names[0][0],"test_0_fM.txt");
	sprintf(table_names[1][0],"test_1_fM.txt");
	sprintf(table_names[2][0],"test_2_fM.txt");
	sprintf(table_names[3][0],"test_0.5_fM.txt");	
	sprintf(table_names[0][1],"test_0_gM.txt");
	sprintf(table_names[1][1],"test_1_gM.txt");
	sprintf(table_names[2][1],"test_2_gM.txt");
	sprintf(table_names[3][1],"test_0.5_gM.txt");
	
	
	double ***gammas;
	gammas = (double ***)malloc(4*sizeof(double **));
	for(int i = 0;i<4;i++){
		gammas[i] = (double **)malloc(nb_a*sizeof(double *));
		for(int j = 0;j<nb_a;j++){
			gammas[i][j] = (double *)malloc(nb_x*sizeof(double));
		};
	};
	double ***mgammas;
	mgammas = (double ***)malloc(4*sizeof(double **));
	for(int i = 0;i<4;i++){
		mgammas[i] = (double **)malloc(nb_a*sizeof(double *));
		for(int j = 0;j<nb_a;j++){
			mgammas[i][j] = (double *)malloc(nb_x*sizeof(double));
		};
	};
	
	printf("\nReading tables...\n");
	fflush(stdout);
	
	read_tables(table_names,gammas,mgammas);
	
	for(int i = 0;i<4;i++){
		for(int j = 0;j<2;j++){
			free(table_names[i][j]);
		};
	};
	
	printf("Writing new tables...\n");
	fflush(stdout);
	
	char ***new_tables;
	new_tables = (char ***)malloc(4*sizeof(char **));
	for(int i = 0;i<4;i++){
		new_tables[i] = (char **)malloc(2*sizeof(char *));
		for(int j = 0;j<2;j++){
			new_tables[i][j] = (char *)malloc(32*sizeof(char));
		};
	};
	
	sprintf(new_tables[0][0],"spin_0_fM.txt");
	sprintf(new_tables[1][0],"spin_1_fM.txt");
	sprintf(new_tables[2][0],"spin_2_fM.txt");
	sprintf(new_tables[3][0],"spin_0.5_fM.txt");
	sprintf(new_tables[0][1],"spin_0_gM.txt");
	sprintf(new_tables[1][1],"spin_1_gM.txt");
	sprintf(new_tables[2][1],"spin_2_gM.txt");
	sprintf(new_tables[3][1],"spin_0.5_gM.txt");
	
	double *x = (double *)malloc(nb_x*sizeof(double));
	for(int i = 0;i<nb_x;i++){
		x[i] = pow(10.,log10(0.01) + (log10(5.) - log10(0.01))/(nb_x - 1)*i);
	};
	
	double *a = (double *)malloc(nb_a*sizeof(double));
	a[0] = 0.;
	for(int i = 30;i<nb_a;i++){
		a[50+30-i-1] = 1. - pow(10.,log10(0.0001) + (log10(0.1) - log10(0.0001))/(nb_a-1 - 30)*(i - 30));
	};
	for(int i = 1;i<9;i++){
		a[i] = pow(10.,log10(0.001) + (log10(0.1) - log10(0.001))/7.*(i - 1));
	};
	for(int i = 8;i<31;i++){
		a[i] = 0.1 + (0.9 - 0.1)/22.*(i - 8);
	};
	
	write_tables(new_tables,gammas,mgammas,a,x);
	
	/*
	for(int i = 0;i<4;i++){
		free(new_tables[i]);
		for(int j = 0;j<nb_a;j++){
			free(gammas[i][j]);
			free(mgammas[i][j]);
		};
	};
	free(x);
	*/
	
	printf("Reading fits...\n");
	fflush(stdout);
	
	char ***fits_tables;
	fits_tables = (char ***)malloc(4*sizeof(char **));
	for(int i = 0;i<4;i++){
		fits_tables[i] = (char **)malloc(2*sizeof(char *));
		for(int j = 0;j<2;j++){	
			fits_tables[i][j] = (char *)malloc(32*sizeof(char));
		};
	};
	
	sprintf(fits_tables[0][0],"fits0_fM.txt");
	sprintf(fits_tables[1][0],"fits1_fM.txt");
	sprintf(fits_tables[2][0],"fits2_fM.txt");
	sprintf(fits_tables[3][0],"fits0.5_fM.txt");
	sprintf(fits_tables[0][1],"fits0_gM.txt");
	sprintf(fits_tables[1][1],"fits1_gM.txt");
	sprintf(fits_tables[2][1],"fits2_gM.txt");
	sprintf(fits_tables[3][1],"fits0.5_gM.txt");
	
	double ****fits;
	fits = (double ****)malloc(4*sizeof(double ***));
	for(int i = 0;i<4;i++){
		fits[i] = (double ***)malloc(2*sizeof(double **));
		for(int j = 0;j<2;j++){
			fits[i][j] = (double **)malloc(nb_a*sizeof(double *));
			for(int k = 0;k<nb_a;k++){
				fits[i][j][k] = (double *)malloc(7*sizeof(double));
			};
		};
	};
	
	read_fits(fits_tables,fits);
	
	printf("Writing new fits...\n");
	fflush(stdout);
	
	char ***new_fits;
	new_fits = (char ***)malloc(4*sizeof(char **));
	for(int i = 0;i<4;i++){
		new_fits[i] = (char **)malloc(2*sizeof(char *));
		for(int j = 0;j<2;j++){
			new_fits[i][j] = (char *)malloc(32*sizeof(char));
		};
	};
	
	sprintf(new_fits[0][0],"spin_0_fits_fM.txt");
	sprintf(new_fits[1][0],"spin_1_fits_fM.txt");
	sprintf(new_fits[2][0],"spin_2_fits_fM.txt");
	sprintf(new_fits[3][0],"spin_0.5_fits_fM.txt");
	sprintf(new_fits[0][1],"spin_0_fits_gM.txt");
	sprintf(new_fits[1][1],"spin_1_fits_gM.txt");
	sprintf(new_fits[2][1],"spin_2_fits_gM.txt");
	sprintf(new_fits[3][1],"spin_0.5_fits_gM.txt");
	
	a = (double *)malloc(nb_a*sizeof(double));
	a[0] = 0.;
	for(int i = 30;i<nb_a;i++){
		a[50+30-i-1] = 1. - pow(10.,log10(0.0001) + (log10(0.1) - log10(0.0001))/(nb_a-1 - 30)*(i - 30));
	};
	for(int i = 1;i<9;i++){
		a[i] = pow(10.,log10(0.001) + (log10(0.1) - log10(0.001))/7.*(i - 1));
	};
	for(int i = 8;i<31;i++){
		a[i] = 0.1 + (0.9 - 0.1)/22.*(i - 8);
	};
	
	write_fits(new_fits,fits,a);
	
	for(int i = 0;i<4;i++){
		free(new_fits[i]);
		for(int j = 0;j<2;j++){
			for(int k = 0;k<nb_a;k++){
				free(fits[i][j][k]);
			};
		};
	};
	
	return 1;
};















