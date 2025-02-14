// Program that launches BlackHawk with a list of parameters files
// that it also creates.
// Author : Jérémy Auffinger, j.auffinger@ipnl.in2p3.fr
// Last modification: 4 January 2021

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
//#include <direct.h>
#include <dirent.h>
#include <string.h>
#include <unistd.h>

#define pi 3.141592653

int nb_params = 0; // number of parameters this program needs

struct param{ // this structure will contain all the parameters of the run
	char destination_folder[200]; // destination folder of the output
	int full_output; // 0: reduced output, 1: extensive output
	int interpolation_method; // 0: we perform linear interpolations, 1: we perform logarithmic interpolations
	
	int metric; // 0: Kerr BHs, 1: polymerized BHs, 2: charged BHs, 3: higher-dimensional BHs
	
	int BH_number; // number of initial black hole masses
	double Mmin; // minimum initial mass of black holes
	double Mmax; // maximum initial mass of black holes
	int param_number; // number of black holes parameters
	double amin; // minimum initial black hole spin
	double amax; // maximum initial black hole spin
	double Qmin; // minimum initial black hole charge
	double Qmax; // maximum initial black hole charge
	double epsilon_LQG; // the epsilon parameter for the polymerized metric (epsilon < 0.794)
	double a0_LQG; // the a0 parameter for the polymerized metric (a0 = 0. or 0.11)
	double n; // the number of additional spatial dimensions
	
	int spectrum_choice; // chooses the initial black holes mass distribution. 0: Dirac, 1: lognormal distribution (mass), 11: lognormal distribution (number), 2: power-law distribution, 
							// 3: critical collapse distribution, 4: peak-theory distribution, 5:uniform, -1: User's distribution
	int spectrum_choice_param; // chooses the initial black holes spin distribution. 0: Dirac, 1: uniform, 2: gaussian.
	double amplitude_lognormal; // amplitude of the distribution
	double amplitude_lognormal2; // amplitude of the distribution
	double stand_dev_lognormal; // standard deviation of the distribution
	double crit_mass_lognormal; // mean of the distribution
	
	double amplitude_powerlaw; // amplitude of the distribution
	double eqstate_powerlaw; // equation of state parameter P = w*rho
	
	double amplitude_critical_collapse; // amplitude of the distribution
	double crit_mass_critical_collapse; // mean of the distribution
	
	double amplitude_uniform; // amplitude of the distribution
	
	double stand_dev_param_gaussian; // standard deviation of the parameter distribution
	double mean_param_gaussian; // mean of the parameter distribution
	
	char table[32]; // name of the User's mass/spin distribution file
	
	int tmin_manual; // 0: automatic tmin for Dirac distribution, 1: manual tmin
	double tmin; // time of black hole formation
	int limit; // maximum iteration limit when computing the life evolution of black holes
	int BH_remnant; // 0: total evaporation after Planck mass, 1: stop the evaporation at M_relic
	double M_remnant; // mass of the BH remnant
	
	int E_number; // number of initial energies
	double Emin; // minimum initial energy
	double Emax; // maximum initial energy
	
	int grav; // 0: no gravitons emitted, 1: gravitons emitted
	int add_DM; // 0: no DM, 1: one primary DM particle
	double m_DM; // DM mass in GeV
	double spin_DM; // spin of DM among 0., 1., 2., 0.5, 1.5
	double dof_DM; // number of DM degrees of freedom
	
	int primary_only; // if set to 1, the code will skip the hadronization part, if set to 0, the secondary spectra are computed
	int hadronization_choice; // choice of hadronization tables (0: PYTHIA, 1: HERWIG, 2: PYTHIA NEW, 3: HAZMA)
};

void create_param_file(char *path,struct param *parameters){
	/*char dir_name[32];
	sprintf(dir_name,"parameters");
	DIR *directory = opendir(dir_name);*/
	char *command = (char *)malloc(256*sizeof(char));
	sprintf(command,"cp ./parameters.txt %s",path);
	system(command);
	FILE *param_file;
	param_file = fopen(path,"r+");
	if(!param_file){
		printf("\n\t [create_param_file] : ERROR could not open file %s",path);
		fflush(stdout);
		exit(0);
	};
	rewind(param_file);
	int counter = 1;
	char ch;
	while((ch = fgetc(param_file)) != EOF && counter < 45){
		if(ch == '&'){
			switch(counter){
				case 1:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%s #",parameters->destination_folder);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 2:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%i #",parameters->full_output);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 3:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%i #",parameters->interpolation_method);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 4:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%i #",parameters->metric);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 5:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%i #",parameters->BH_number);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 6:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->Mmin);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 7:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->Mmax);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 8:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%i #",parameters->param_number);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 9:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->amin);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 10:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->amax);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 11:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->Qmin);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 12:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->Qmax);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 13:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->epsilon_LQG);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 14:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->a0_LQG);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 15:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->n);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 16:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%i #",parameters->spectrum_choice);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 17:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%i #",parameters->spectrum_choice_param);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 18:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->amplitude_lognormal);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 19:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->amplitude_lognormal2);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 20:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->stand_dev_lognormal);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 21:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->crit_mass_lognormal);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 22:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->amplitude_powerlaw);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 23:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->eqstate_powerlaw);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 24:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->amplitude_critical_collapse);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 25:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->crit_mass_critical_collapse);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 26:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->amplitude_uniform);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 27:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->stand_dev_param_gaussian);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 28:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->mean_param_gaussian);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 29:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%s #",parameters->table);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 30:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%i #",parameters->tmin_manual);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 31:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->tmin);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 32:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%i #",parameters->limit);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 33:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%i #",parameters->BH_remnant);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 34:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->M_remnant);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 35:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%i #",parameters->E_number);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 36:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->Emin);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 37:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->Emax);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 38:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%i #",parameters->grav);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 39:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%i #",parameters->add_DM);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 40:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->m_DM);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 41:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->spin_DM);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 42:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%.5e #",parameters->dof_DM);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 43:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%i #",parameters->primary_only);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				case 44:{
					fseek(param_file,-1,SEEK_CUR);
					fprintf(param_file,"%i #",parameters->hadronization_choice);
					fseek(param_file,0,SEEK_CUR);
					counter++;
					break;
				};
				default:{
					printf("\n\t [create_param_file] : ERROR in the counter value: %i",counter);
					fflush(stdout);
					exit(0);
					break;
				};
			};
		};
	};
	fclose(param_file);
	
	return;
};

void fill_parameters_default(struct param *parameters,char *run_name){
	// this routine fills a parameters structure with default parameters
	
	sprintf(parameters->destination_folder,run_name); // destination folder of the output
	parameters->full_output = 0; // 0: reduced output, 1: extensive output
	parameters->interpolation_method = 0; // 0: we perform linear interpolations, 1: we perform logarithmic interpolations
	
	parameters->metric = 0; // 0: Kerr BHs, 1: polymerized BHs, 2: charged BHs, 3: higher-dimensional BHs
	
	parameters->BH_number = 1; // number of initial black hole masses
	parameters->Mmin = 1e+15; // minimum initial mass of black holes
	parameters->Mmax = 1e+16; // maximum initial mass of black holes
	parameters->param_number = 1; // number of black holes parameters
	parameters->amin = 0.; // minimum initial black hole spin
	parameters->amax = 0.9999; // maximum initial black hole spin
	parameters->Qmin = 0.; // minimum initial black hole charge
	parameters->Qmax = 0.9; // maximum initial black hole charge
	parameters->epsilon_LQG = 0.; // the epsilon parameter for the polymerized metric (epsilon < 0.794)
	parameters->a0_LQG = 0.; // the a0 parameter for the polymerized metric (a0 = 0. or 0.11)
	parameters->n = 0; // the number of additional spatial dimensions
	
	parameters->spectrum_choice = 0; // chooses the initial black holes mass distribution. 0: Dirac, 1: lognormal distribution (mass), 11: lognormal distribution (number), 2: power-law distribution, 
							// 3: critical collapse distribution, 4: peak-theory distribution, 5:uniform, -1: User's distribution
	parameters->spectrum_choice_param = 0; // chooses the initial black holes spin distribution. 0: Dirac, 1: uniform, 2: gaussian.
	parameters->amplitude_lognormal = 1.; // amplitude of the distribution
	parameters->amplitude_lognormal2 = 1.; // amplitude of the distribution
	parameters->stand_dev_lognormal = 1.; // standard deviation of the distribution
	parameters->crit_mass_lognormal = 1.; // mean of the distribution
	
	parameters->amplitude_powerlaw = 1.; // amplitude of the distribution
	parameters->eqstate_powerlaw = 0.3333; // equation of state parameter P = w*rho
	
	parameters->amplitude_critical_collapse = 1.; // amplitude of the distribution
	parameters->crit_mass_critical_collapse = 1.; // mean of the distribution
	
	parameters->amplitude_uniform = 1.; // amplitude of the distribution
	
	parameters->stand_dev_param_gaussian = 1.; // standard deviation of the parameter distribution
	parameters->mean_param_gaussian = 1.; // mean of the parameter distribution
	
	sprintf(parameters->table,"table.txt"); // name of the User's mass/spin distribution file
	
	parameters->tmin_manual = 0; // 0: automatic tmin for Dirac distribution, 1: manual tmin
	parameters->tmin = 1e-30; // time of black hole formation
	parameters->limit = 5000; // maximum iteration limit when computing the life evolution of black holes
	parameters->BH_remnant = 0; // 0: total evaporation after Planck mass, 1: stop the evaporation at M_relic
	parameters->M_remnant = 1e-4; // mass of the BH remnant
	
	parameters->E_number = 1000; // number of initial energies
	parameters->Emin = 5; // minimum initial energy
	parameters->Emax = 1e+5; // maximum initial energy
	
	parameters->grav = 0; // 0: no gravitons emitted, 1: gravitons emitted
	parameters->add_DM = 0; // 0: no DM, 1: one primary DM particle
	parameters->m_DM = 0.; // DM mass in GeV
	parameters->spin_DM = 0.; // spin of DM among 0., 1., 2., 0.5, 1.5
	parameters->dof_DM = 1.; // number of DM degrees of freedom
	
	parameters->primary_only = 0; // if set to 1, the code will skip the hadronization part, if set to 0, the secondary spectra are computed
	parameters->hadronization_choice = 0; // choice of hadronization tables (0: PYTHIA, 1: HERWIG, 2: PYTHIA NEW, 3: HAZMA)
	
	return;
}

void fill_parameters_variable(struct param *parameters,int i,double parammin,double parammax,int nb_runs,int choice_variable_param,int scan_type){
	// this routine modifies some variable parameter in the parameters structure
	
	double temp;
	if(scan_type == 0 || choice_variable_param == 5){ // linear scan
		temp = parammin + (parammax - parammin)/(nb_runs - 1)*i;
	}
	else{
		temp = pow(10.,log10(parammin) + (log10(parammax) - log10(parammin))/(nb_runs - 1)*i);
	}
	switch(choice_variable_param){
		case 1:{
			parameters->Mmin = temp;
			break;
		}
		case 2:{
			parameters->amin = temp;
			break;
		}
		case 3:{
			parameters->Qmin = temp;
			break;
		}
		case 4:{
			parameters->epsilon_LQG = temp;
			break;
		}
		case 5:{
			parameters->n = temp;
			break;
		}
		default:{
			printf("\n\t [fill_parameters_variable] : ERROR wrong choice of variable parameter !");
			fflush(stdout);
			exit(0);
			break;
		}
	}
	
	return;
}

void fill_parameters_hadronization(struct param *parameters,int hadronization_choice){
	// this routine fills the hadronization parameter and adapts the primary energies
	
	parameters->hadronization_choice = hadronization_choice-1;
	switch(hadronization_choice){
		case 1:case 3:{ // PYTHIA tables, PYTHIA tables new // HERE extrapolation !!!
			parameters->Emin = 1e-6;
			parameters->Emax = 1e+5;
			break;
		}
		case 2:{ // HERWIG tables
			parameters->Emin = 25.;
			parameters->Emax = 1e+5;
			break;
		}
		case 4:{ // Hazma tables
			parameters->Emin = 1e-6;
			parameters->Emax = 5.;
			break;
		}
		default:{
			printf("\n\t\t [fill_parameters_hadronization] : ERROR wrong hadronization choice %i !",hadronization_choice);
			fflush(stdout);
			exit(0);
			break;
		}
	}
	
	return;
}

void create_runs_file(char **runs,int nb_runs,char *runs_name){
	// this routine creates a file "runs.txt" to be read by constraints.x
	
	char *file_name = (char *)malloc(64*sizeof(char));
	sprintf(file_name,runs_name);
	FILE *file = fopen(file_name,"w+");
	if(!file){
		printf("\n\t [create_runs_file] : ERROR could not open file %s !",file_name);
		fflush(stdout);
		exit(0);
	}
	rewind(file);
	fprintf(file,"nb_runs = %i\n\n",nb_runs);
	for(int i = 0;i<nb_runs;i++){
		fprintf(file,"%s\n",runs[i]);
	}
	fclose(file);
	
	return;
}

void BH_launcher(char **param_names,int nb_runs,int sessions){
	// this routine launches BlackHawk_inst and BlackHawk_tot on the parameters files created before
	
	chdir("../../..");
	char command[256];
	if(sessions == 0 || sessions == 2){
		printf("\n\t launching runs with BlackHawk_inst...");
		fflush(stdout);
		for(int i = 0;i<nb_runs;i++){
			printf(" %i",i+1);
			fflush(stdout);
			// warning remove "./" in case of cygwin
			sprintf(command,"./BlackHawk_inst.x ../Isatis/BH_launcher/%s > nohup_inst_%i.txt &",param_names[i],i+1);
			system(command);
			fflush(stdout);
		}
		printf(" DONE");
		fflush(stdout);
	}
	if(sessions == 1 || sessions == 2){
		printf("\n\t launching runs with BlackHawk_tot...");
		fflush(stdout);
		for(int i = 0;i<nb_runs;i++){
			printf(" %i",i+1);
			fflush(stdout);
			// warning remove "./" in case of cygwin
			sprintf(command,"./BlackHawk_tot.x ../Isatis/BH_launcher/%s > nohup_tot_%i.txt &",param_names[i],i+1);
			system(command);
			fflush(stdout);
		}
		printf(" DONE");
		fflush(stdout);
	}
	printf("\n");
	fflush(stdout);
	chdir("./scripts/Isatis/BH_launcher/");
	return;
};

int main(int argc, char **argv){
	
	printf("\n\t # Comprehensive BlackHawk launcher device #\n\n");
	printf("[main] : Start of execution\n");
	fflush(stdout);
	
	if(argc < nb_params+1){
		printf("[main] : ERROR this program needs at least %i parameters\n[main] : End of execution",nb_params);
		fflush(stdout);
		exit(0);
	};
	
	char *runs_name = (char *)malloc(128*sizeof(char));
	printf("\t > What is the name of the runs file? Answer: ");
	fflush(stdout);
	scanf("%s",runs_name);
	
	int nb_runs;
	printf("\t > How many runs do you want to launch? Answer: ");
	fflush(stdout);
	scanf("%i",&(nb_runs));

	if(nb_runs < 1){
		printf("[main] : ERROR wrong number of runs %i",nb_runs);
		fflush(stdout);
		exit(0);
	}
	
	int choice_variable_param;
	int scan_type;
	double parammin,parammax;
	if(nb_runs > 1){
		printf("\t > What parameter do you want to vary?");
		printf("\n\t\t 1 for BH mass");
		printf("\n\t\t 2 for BH spin");
		printf("\n\t\t 3 for BH charge");
		printf("\n\t\t 4 for BH LQG parameter");
		printf("\n\t\t 5 for number of dimensions");
		printf("\n\t   Answer: ");
		fflush(stdout);
		scanf("%i",&(choice_variable_param));
		if(choice_variable_param > 0 && choice_variable_param < 5){
			printf("\t > Do you want a linear (type 0) or logarithmic (type 1) scan? Answer: ");
			fflush(stdout);
			scanf("%i",&(scan_type));
		}
		switch(choice_variable_param){
			case 1:{
				printf("\t > What is the minimum mass (in g)? Answer: ");
				fflush(stdout);
				scanf("%lf",&(parammin));
				printf("\t > What is the maximum mass (in g)? Answer: ");
				fflush(stdout);
				scanf("%lf",&(parammax));
				break;
			}
			case 2:{
				printf("\t > What is the minimum spin (dimensionless)? Answer: ");
				fflush(stdout);
				scanf("%lf",&(parammin));
				printf("\t > What is the maximum spin (dimensionless)? Answer: ");
				fflush(stdout);
				scanf("%lf",&(parammax));
				break;
			}
			case 3:{
				printf("\t > What is the minimum charge (dimensionless)? Answer: ");
				fflush(stdout);
				scanf("%lf",&(parammin));
				printf("\t > What is the maximum charge (dimensionless)? Answer: ");
				fflush(stdout);
				scanf("%lf",&(parammax));
				break;
			}
			case 4:{
				printf("\t > What is the minimum epsilon_0 (dimensionless)? Answer: ");
				fflush(stdout);
				scanf("%lf",&(parammin));
				printf("\t > What is the maximum epsilon_0 (dimensionless)? Answer: ");
				fflush(stdout);
				scanf("%lf",&(parammax));
				break;
			}
			case 5:{
				printf("\t > What is the minimum number of dimensions (dimensionless)? Answer: ");
				fflush(stdout);
				scanf("%lf",&(parammin));
				printf("\t > What is the maximum number of dimensions (dimensionless)? Answer: ");
				fflush(stdout);
				scanf("%lf",&(parammax));
				break;
			}
			default:{
				printf("[main] : ERROR wrong choice of variable parameter!");
				fflush(stdout);
				exit(0);
				break;
			}
		}
	}
	
	int hadronization_choice;
	printf("\t > What hadronization tables do you want to use?");
	printf("\n\t\t 1 for PYTHIA (BBN) tables");
	printf("\n\t\t 2 for HERWIG (BBN) tables");
	printf("\n\t\t 3 for PYTHIA (today) tables");
	printf("\n\t\t 4 for Hazma (today) tables");
	printf("\n\t   Answer: ");
	fflush(stdout);
	scanf("%i",&(hadronization_choice));
	
	int sessions;
	printf("\t > What instance do you want to launch?");
	printf("\n\t\t 1 for BlackHawk_inst");
	printf("\n\t\t 2 for BlackHawk_tot");
	printf("\n\t\t 3 for both");
	printf("\n\t   Answer: ");
	fflush(stdout);
	scanf("%i",&(sessions));
	
	char basic_name[64];
	printf("\t > Enter the basic run name: ");
	fflush(stdout);
	scanf("%s",basic_name);
	
	printf("[main] : Creating parameters files...");
	fflush(stdout);
	
	char **runs = (char **)malloc(nb_runs*sizeof(char *));
	for(int i = 0;i<nb_runs;i++){
		runs[i] = (char *)malloc(256*sizeof(char));
	};
	
	for(int i = 0;i<nb_runs;i++){
		sprintf(runs[i],"%s_%i",basic_name,i+1);
	}

	char **param_names = (char **)malloc(nb_runs*sizeof(char *));
	for(int i = 0;i<nb_runs;i++){
		param_names[i] = (char *)malloc(256*sizeof(char));
		sprintf(param_names[i],"%s.txt",runs[i]);
	};
	
	struct param *parameters = (struct param *)malloc(nb_runs*sizeof(struct param));
	
	for(int i = 0;i<nb_runs;i++){
		fill_parameters_default(&(parameters[i]),runs[i]); // we fill using the default parameters
		fill_parameters_hadronization(&(parameters[i]),hadronization_choice); // we modify the hadronization parameters
		if(nb_runs > 1){
			fill_parameters_variable(&(parameters[i]),i,parammin,parammax,nb_runs,choice_variable_param,scan_type); // for scan of parameters, we modify the parameter
		}
		create_param_file(param_names[i],&(parameters[i])); // we create the parameters files
		printf(" %i",i+1);
		fflush(stdout);
	}
	
	printf(" DONE\n");
	printf("[main] : Creating runs file...");
	fflush(stdout);
	
	create_runs_file(runs,nb_runs,runs_name);
	
	printf(" DONE\n");
	
	printf("[main] : Launching BlackHawk...");
	fflush(stdout);
	
	BH_launcher(param_names,nb_runs,sessions-1);
	
	// warning remove "\n" in case of cygwin
	printf("[main] : End of execution\n");
	fflush(stdout);
	
	return 1;
};











