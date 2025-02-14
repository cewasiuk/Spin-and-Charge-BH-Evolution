// Script to compute multiple HERWIG hadronization tables

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <stdbool.h>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <unistd.h>
using namespace std;

int nb_events = 100000;
double Emin = 25.;
double Emax = 100000.;
int nb_init_en = 100;
int nb_fin_en = 100;

int main(){

	//system("rm */ -rf");
	
	double *initial_energies = new double [nb_init_en];
	for(int i = 0;i<nb_init_en;i++){
		initial_energies[i] = pow(10.,log10(Emin) + (log10(Emax) - log10(Emin))/(nb_init_en-1)*i);
	};
	
	// creating the output folders
	
	cout << "[main] : Creation of the output folders..." << endl;
	
	char *command = new char [200];
	char **folders = new char *[nb_init_en];
	for(int i = 0;i<nb_init_en;i++){
		folders[i] = new char [20];
		sprintf(folders[i],"%.5eGeV",initial_energies[i]);
		sprintf(command,"mkdir %s",folders[i]);
		system(command);
	};
	
	cout << "[main] : Folders ...GeV created" << endl;
	
	// defining input files
	
	cout << "[main] : Defining input files..." << endl;
	
	ifstream in("LEP.in",ios::in);
	if(!in){
		cout << "[main] : ERROR cannot open LEP.in file" << endl;
		return 0;
	};
	
	char *new_input_file_path = new char [40];
	string str;
	int cond;
	
	for(int i = 0;i<nb_init_en;i++){
		
		ifstream in("LEP.in",ios::in);
		if(!in){
			cout << "[main] : ERROR cannot open LEP.in file" << endl;
			return 0;
		};
	
		sprintf(new_input_file_path,"./%s/LEP.in",folders[i]);
		
		ofstream out(new_input_file_path,ios::out);
		while(in){
			cond = 0;
			getline(in,str);
			for(int j = 0;j<80;j++){
				if(str[j] == '%'){
					str[j] = ' ';
					cond = 1;
				};
			};
			if(cond){ // in this case we are on the right line : we add the right energy
				out << str << initial_energies[i] << endl;
			}
			else{
				out << str << endl;
			};
		};
		out.close();
		in.close();
	};
	
	cout << "[main] : Input files LEP.in created" << endl;
	
	// running the computation
	
	cout << "[main] : Computing the events..." << endl;
	
	for(int i = 0;i<nb_init_en;i++){
		chdir(folders[i]);
		sprintf(command,"Herwig read LEP.in");
		system(command);
		sprintf(command,"Herwig run LEP.run -N %i",nb_events);
		system(command);
		chdir("..");
	};
	
	cout << "[main] : Computation done" << endl;
	
	return 1;
};
