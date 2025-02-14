// This is the source file where the methods computing
// the time evolution of black holes are implemented.
// Last modification: 07 April 2024
// Authors: Jérémy Auffinger jeremy.auffinger.blackhawk@protonmail.com & Alexandre Arbey alexandre.arbey@ens-lyon.fr

#include "include.h"

double rplus_Kerr(double M,double a){
	// This function computes the external Kerr radius
	// of a BH of mass M and spin parameter a.
	
	return M*(1. + sqrt(1. - pow(a,2.)));
}

double temp_Kerr(double M,double a){
	// This function computes the Hawking temperature of a Kerr
	// PBH of mass M and spin parameter a.
	
	return 1./(2.*pi*G)*(rplus_Kerr(M,a) - M)/(pow(rplus_Kerr(M,a),2.) + pow(a*M,2.));
}

double P_LQG(double epsilon){
	// This function computes the polymerization fonction P(epsilon).
	
	return (sqrt(1. + pow(epsilon,2.)) - 1.)/(sqrt(1. + pow(epsilon,2.)) + 1.);
};

double m_LQG(double M,double epsilon){
	// This function computes the modified mass m in polymeriezd metrics.
	
	return M/pow(1. + P_LQG(epsilon),2.);
};

double temp_LQG(double M,double epsilon,double a0){
	// This function computes the Hawking temperature of a polymerized BH
	// with parameters epsilon and a0.
	
	return 4.*pow(m_LQG(M,epsilon),3.)*(1-pow(P_LQG(epsilon),2.))/(32.*pi*pow(m_LQG(M,epsilon),4.) + 2.*pi*pow(a0,2.))/G;
};

double rplus_charged(double M,double Q){
	// This function computes the external radius of a charged BH.
	
	return M*(1. + sqrt(1-pow(Q,2.)));
};

double rminus_charged(double M,double Q){
	// This function computes the external radius of a charged BH.
	
	return M*(1-sqrt(1-pow(Q,2.)));
};

double temp_charged(double M, double Q){
	// This function computes the Hawking temperature of a charged BH with dimensionless charge Q.
	
	return (rplus_charged(M,Q) - rminus_charged(M,Q))/(4.*pi*pow(rplus_charged(M,Q),2.))/G;
};

double rH_higher(double M,double n,double M_star){
	// This function computes the radius of a higher dimensional BH with n extra dimensions,
	// mass M and modified Planck mass M_star.
	
	double gamma;
	if(n == 0.){
		gamma = sqrt(pi)/2.; // Gamma(3/2)
	}
	else if(n == 1.){
		gamma = 1.; // Gamma(2)
	}
	else if(n == 2.){
		gamma = 3.*sqrt(pi)/4.; // Gamma(5/2)
	}
	else if(n == 3.){
		gamma = 2.; // Gamma(3)
	}
	else if(n == 4.){
		gamma = 15./8.*sqrt(pi); // Gamma(7/2)
	}
	else if(n == 5.){
		gamma = 6.; // Gamma(4)
	}
	else if(n == 6.){
		gamma = 105./16.*sqrt(pi); // Gamma(9/2)
	}
	else{
		printf("\n\t [rH_higher] : ERROR WRONG NUMBER OF EXTRA DIMENSIONS !\n");
		fflush(stdout);
		exit(0);
	}
	
	return 1./(sqrt(pi)*M_star)*pow(M/M_star*8.*gamma/(n+2.),1./(n+1.));
};

double temp_higher(double M,double n,double M_star){
	// This function computes the Hawking temperature of a higher dimensional BH
	// with n extra dimensions and a modified Planck mass M_star.
	
	return (n+1.)/(4.*pi*rH_higher(M,n,M_star))/G;
};

void read_fM_table(double **fM_table,double *fM_masses,double *fM_param,struct param *parameters){
	// This function reads the f(M,param) factor table used to compute the evolution
	// of BH, stored in the folder /fM_tables.
	
	FILE *table_file;
	char name[64];
	if(parameters->metric == 0){
		if(!parameters->add_DM){
			if(parameters->grav){
				sprintf(name,"./src/tables/fM_tables/fM.txt");
				//sprintf(name,"./src/tables/fM_tables/fM_marco.txt"); // HERE
			}
			else{
				sprintf(name,"./src/tables/fM_tables/fM_nograv.txt");
			}
		}
		else{
			if(parameters->grav){
				if(parameters->spin_DM == 0.){
					sprintf(name,"./src/tables/fM_tables/fM_add0.txt");
				}
				else if(parameters->spin_DM == 1.){
					sprintf(name,"./src/tables/fM_tables/fM_add1.txt");
				}
				else if(parameters->spin_DM == 2.){
					sprintf(name,"./src/tables/fM_tables/fM_add2.txt");
				}
				else if(parameters->spin_DM == 0.5){
					sprintf(name,"./src/tables/fM_tables/fM_add05.txt");
				}
				else{
					sprintf(name,"./src/tables/fM_tables/fM_add15.txt");
				}
			}
			else{
				if(parameters->spin_DM == 0.){
					sprintf(name,"./src/tables/fM_tables/fM_add0_nograv.txt");
				}
				else if(parameters->spin_DM == 1.){
					sprintf(name,"./src/tables/fM_tables/fM_add1_nograv.txt");
				}
				else if(parameters->spin_DM == 2.){
					sprintf(name,"./src/tables/fM_tables/fM_add2_nograv.txt");
				}
				else if(parameters->spin_DM == 0.5){
					sprintf(name,"./src/tables/fM_tables/fM_add05_nograv.txt");
				}
				else{
					sprintf(name,"./src/tables/fM_tables/fM_add15_nograv.txt");
				}
			}
		}
	}
	else if(parameters->metric == 1){
		if(parameters->a0_LQG == 0.){
			if(parameters->grav){
				sprintf(name,"./src/tables/fM_tables/fM_LQG.txt");
			}
			else{
				sprintf(name,"./src/tables/fM_tables/fM_LQG_nograv.txt");
			}
		}
		else{ // by default the a0 != 0 case corresponds to the standard a0 = 0.11
			if(parameters->grav){
				sprintf(name,"./src/tables/fM_tables/fM_LQG_a0.txt");
			}
			else{
				sprintf(name,"./src/tables/fM_tables/fM_LQG_a0_nograv.txt");
			}
		}
	}
	table_file = fopen(name,"r");
	if(!table_file){
		printf("\n\t [read_fM_table] : ERROR COULD NOT READ TABLE '%s' !\n",name);
		fflush(stdout);
		exit(0);
	}
	rewind(table_file);
	if(parameters->full_output){
		printf("\n\n\t reading:\t%s\n",name);
		fflush(stdout);
	}
	char dummy[20];
	fscanf(table_file,"%s",dummy);
	for(int i = 0;i<parameters->nb_fM_param;i++){
		fscanf(table_file,"%lf",&(fM_param[i]));
	}
	for(int i = 0;i<parameters->nb_fM_masses;i++){
		fscanf(table_file,"%lf",&(fM_masses[i]));
		for(int j = 0;j<parameters->nb_fM_param;j++){
			fscanf(table_file,"%lf",&(fM_table[i][j]));
		}
	}
	fclose(table_file);
	return;
}

void read_gM_table(double **gM_table,double *fM_masses,double *fM_param,struct param *parameters){
	// This function reads the g(M,param) factor table used to compute the evolution
	// of BH, stored in the folder /fM_tables.
	
	FILE *table_file;
	char name[64];
	if(!parameters->add_DM){
		if(parameters->grav){
			sprintf(name,"./src/tables/fM_tables/gM.txt");
		}
		else{
			sprintf(name,"./src/tables/fM_tables/gM_nograv.txt");
		}
	}
	else{
		if(parameters->grav){
			if(parameters->spin_DM == 0.){
				sprintf(name,"./src/tables/fM_tables/gM_add0.txt");
			}
			else if(parameters->spin_DM == 1.){
				sprintf(name,"./src/tables/fM_tables/gM_add1.txt");
			}
			else if(parameters->spin_DM == 2.){
				sprintf(name,"./src/tables/fM_tables/gM_add2.txt");
			}
			else if(parameters->spin_DM == 0.5){
				sprintf(name,"./src/tables/fM_tables/gM_add05.txt");
			}
			else{
				sprintf(name,"./src/tables/fM_tables/gM_add15.txt");
			}
		}
		else{
			if(parameters->spin_DM == 0.){
				sprintf(name,"./src/tables/fM_tables/gM_add0_nograv.txt");
			}
			else if(parameters->spin_DM == 1.){
				sprintf(name,"./src/tables/fM_tables/gM_add1_nograv.txt");
			}
			else if(parameters->spin_DM == 2.){
				sprintf(name,"./src/tables/fM_tables/gM_add2_nograv.txt");
			}
			else if(parameters->spin_DM == 0.5){
				sprintf(name,"./src/tables/fM_tables/gM_add05_nograv.txt");
			}
			else{
				sprintf(name,"./src/tables/fM_tables/gM_add15_nograv.txt");
			}
		}
	}
	table_file = fopen(name,"r");
	if(!table_file){
		printf("\n\t [read_gM_table] : ERROR COULD NOT READ TABLE '%s' !\n",name);
		fflush(stdout);
		exit(0);
	}
	rewind(table_file);
	if(parameters->full_output){
		printf("\t\t\t%s\n\n",name);
		fflush(stdout);
	}
	char dummy[20];
	fscanf(table_file,"%s",dummy);
	for(int i = 0;i<parameters->nb_fM_param;i++){
		fscanf(table_file,"%lf",&(fM_param[i]));
	}
	for(int i = 0;i<parameters->nb_fM_masses;i++){
		fscanf(table_file,"%lf",&(fM_masses[i]));
		for(int j = 0;j<parameters->nb_fM_param;j++){
			fscanf(table_file,"%lf",&(gM_table[i][j]));
		}
	}
	fclose(table_file);
	return;
}

double loss_rate_M(double M,double param,double **fM_table,double *fM_masses,double *fM_param,int counter_M,int counter_param,struct param *parameters){
	// This function computes the BH mass loss rate dM/dt by using the f table.
	
	double fM;
	if(parameters->interpolation_method == 0){ // linear interpolation
		fM = fM_table[counter_M][counter_param] + (fM_table[counter_M+1][counter_param] - fM_table[counter_M][counter_param])/(fM_masses[counter_M+1] - fM_masses[counter_M])*(M - fM_masses[counter_M]) + (fM_table[counter_M][counter_param+1] - fM_table[counter_M][counter_param])/(fM_param[counter_param+1] - fM_param[counter_param])*(param - fM_param[counter_param]);
	}
	else if(parameters->interpolation_method == 1){ // logarithmic interpolation
		fM = pow(10.,log10(fM_table[counter_M][counter_param]) + (log10(fM_table[counter_M+1][counter_param]) - log10(fM_table[counter_M][counter_param]))/(log10(fM_masses[counter_M+1]) - log10(fM_masses[counter_M]))*(log10(M) - log10(fM_masses[counter_M]))) + (fM_table[counter_M][counter_param+1] - fM_table[counter_M][counter_param])/(fM_param[counter_param+1] - fM_param[counter_param])*(param - fM_param[counter_param]);
	}
	else{
		printf("\n\t [loss_rate_M] : ERROR WRONG INTERPOLATION METHOD !\n");
		fflush(stdout);
		exit(1);
	}
	return -fM/pow(M,2.);
}

double loss_rate_param(double M,double param,double **fM_table,double **gM_table,double *fM_masses,double *fM_param,int counter_M,int counter_param,struct param *parameters){
	// This function computes the BH spin loss rate da/dt by using the f(M,a) and g(M,a) tables.
	
	if(parameters->metric == 0){
		if(param == 0.){ // no more evolution of the spin
			return 0.;
		}
		else{
			double fM;
			double gM;
			if(parameters->interpolation_method == 0 || counter_param == 0){ // linear interpolation
				fM = fM_table[counter_M][counter_param] + (fM_table[counter_M+1][counter_param] - fM_table[counter_M][counter_param])/(fM_masses[counter_M+1] - fM_masses[counter_M])*(M - fM_masses[counter_M]) + (fM_table[counter_M][counter_param+1] - fM_table[counter_M][counter_param])/(fM_param[counter_param+1] - fM_param[counter_param])*(param - fM_param[counter_param]);
				gM = gM_table[counter_M][counter_param] + (gM_table[counter_M+1][counter_param] - gM_table[counter_M][counter_param])/(fM_masses[counter_M+1] - fM_masses[counter_M])*(M - fM_masses[counter_M]) + (gM_table[counter_M][counter_param+1] - gM_table[counter_M][counter_param])/(fM_param[counter_param+1] - fM_param[counter_param])*(param - fM_param[counter_param]);
			}
			else if(parameters->interpolation_method == 1){ // logarithmic interpolation
				fM = pow(10.,log10(fM_table[counter_M][counter_param]) + (log10(fM_table[counter_M+1][counter_param]) - log10(fM_table[counter_M][counter_param]))/(log10(fM_masses[counter_M+1]) - log10(fM_masses[counter_M]))*(log10(M) - log10(fM_masses[counter_M]))) + (fM_table[counter_M][counter_param+1] - fM_table[counter_M][counter_param])/(fM_param[counter_param+1] - fM_param[counter_param])*(param - fM_param[counter_param]);
				gM = pow(10.,log10(gM_table[counter_M][counter_param]) + (log10(gM_table[counter_M+1][counter_param]) - log10(gM_table[counter_M][counter_param]))/(log10(fM_masses[counter_M+1]) - log10(fM_masses[counter_M]))*(log10(M) - log10(fM_masses[counter_M]))) + (gM_table[counter_M][counter_param+1] - gM_table[counter_M][counter_param])/(fM_param[counter_param+1] - fM_param[counter_param])*(param - fM_param[counter_param]);
			}
			else{
				printf("\n\t [loss_rate_param] : ERROR WRONG INTERPOLATION METHOD !\n");
				fflush(stdout);
				exit(0);
			}
			return -gM*param/pow(M,3.) + 2.*fM*param/pow(M,3.);
		}
	}
	else if(parameters->metric == 1){
		return 0.; // the epsilon parameter is constant
	}
}

void evolution_times(double *init_masses,double *init_params,double **evol_times,double **fM_table,double **gM_table,double *fM_masses,double *fM_param,
	struct param *parameters){
	// This function computes the lifetimes of each BH in the spec_table grid, in order to
	// sort them and predict the correct integration order for life_evolution(). The limit
	// parameter should not be reached.
	
	double life_mass; // running BH mass
	double life_param; // running BH spin
	double dt; // running time step
	int counter_M; // running position in f(M,a) / g(M,a) tables
	int counter_param; // running position in f(M,a) / g(M,a) tables
	double estimate_M; // intermediate evaluations
	double estimate_param; // intermediate evaluations
	int counter; // number of time steps for each mass, should not reach the limit parameter
	double evol_M; // relative BH mass evolution
	double evol_param; // relative BH spin evolution
	
	double limiting_mass;
	if(parameters->BH_remnant){
		limiting_mass = parameters->M_remnant;
	}
	else{
		limiting_mass = 1.1*Mp;
	}
	
	for(int i = 0;i<parameters->BH_number;i++){
		for(int j = 0;j<parameters->param_number;j++){
			life_mass = init_masses[i]; // initialization of the mass
			life_param = init_params[j]; // initialization of the spin
			dt = parameters->tmin; // initialization of the timestep
			evol_times[i][j] = parameters->tmin; // initialization of the lifetime
			counter_M = parameters->nb_fM_masses-1; // the BH mass is decreasing
			counter_param = parameters->nb_fM_param-1; // the BH parameter is decreasing (spin)
			counter = 0;
			while(counter < parameters->limit-1 && life_mass > limiting_mass){ // recursion limit, evaporation until Planck mass
				counter_M = parameters->nb_fM_masses-1;
				counter_param = parameters->nb_fM_param-1;
				while(counter_M > 0 && life_mass < fM_masses[counter_M]){ // we find the closest tabulated value
					counter_M--;
				}
				while(counter_param > 0 && life_param < fM_param[counter_param]){ // we find the closest tabulated value
					counter_param--;
				}
				
				estimate_M = life_mass + dt*loss_rate_M(life_mass,life_param,fM_table,fM_masses,fM_param,counter_M,counter_param,parameters); // Euler estimation of the next mass
				estimate_param = life_param + dt*loss_rate_param(life_mass,life_param,fM_table,gM_table,fM_masses,fM_param,counter_M,counter_param,parameters); // Euler estimation of the next spin
				evol_M = estimate_M/life_mass;
				if(life_param != 0.){
					evol_param = estimate_param/life_param;
				}
				else{
					evol_param = 0.;
				}
				
				if(evol_M < 0.9 || (parameters->metric == 0 && (evol_param != 0. && ((evol_param < 1. && evol_param < 0.9) || (evol_param > 1. && evol_param > 1.1))))){ // one of the variations is too high, we diminish the time step
					dt = dt/2.;
				}
				else if(evol_M > 0.99 && (parameters-> metric != 0 || evol_param == 0. || ((evol_param <= 1. && evol_param > 0.99) || (evol_param > 1. && evol_param < 1.01))) && dt/evol_times[i][j] < 1.){ // all the variations are too low and the time step is not too high compared to the time scale, we increase the time step
					dt = dt*2.;
				}
				else{ // all the variations are fine, we record the evolution and move forward
					evol_times[i][j] = evol_times[i][j] + dt;
					life_mass = estimate_M;
					if(estimate_param < 1.e-3){ // we stop evolving the parameter and approximate it to 0
						life_param = 0.;
					}
					else{
						life_param = estimate_param;
					}
				}
				
				counter++;
			}
			if(counter >= parameters->limit){
				printf("\n\t [evolution_times] : ERROR ITERATION LIMIT REACHED !\n");
				fflush(stdout);
				exit(0);
			}
		}
	}
	
	return;
}

void sort_lifetimes(double **evol_times,double *sorted_times,int **rank,struct param *parameters){
	// This function sorts the lifetime of BHs to use the simple algorithm of
	// life_evolutions().
	
	for(int i = 0;i<parameters->BH_number;i++){
		for(int j = 0;j<parameters->param_number;j++){
			sorted_times[parameters->param_number*i+j] = evol_times[i][j];
		}
	}
	sort_fusion(sorted_times,parameters->BH_number*parameters->param_number);
	int counter_rank;
	for(int i = 0;i<parameters->BH_number;i++){
		for(int j = 0;j<parameters->param_number;j++){
			counter_rank = 0;
			while(counter_rank < parameters->BH_number*parameters->param_number-1 && fabs(evol_times[i][j] - sorted_times[counter_rank])/sorted_times[counter_rank] > machine_precision){
				counter_rank++;
			}
			rank[i][j] = counter_rank;
		}
	}
	
	return;
}

void life_evolution(double ***life_masses,double ***life_params,double *life_times,double *dts,double *init_masses,
	double *init_params,int **rank,double **fM_table,double **gM_table,double *fM_masses,double *fM_param,struct param *parameters){
	// This function computes the evolution of Kerr BHs from their initial mass to the Planck mass.
	// It uses the functions loss_rate_M_Kerr and loss_rate_a_Kerr to compute an estimation of M(t+dt) and a(t+dt) through a Euler
	// method, and adapts the time step to keep relative variations between 0.001 and 0.1.
	// The intermediate masses, spins, times and time intervals are stored in the arrays life_masses[][][],
	// life_params[][][], life_times[] and dts[] respectively. A recursion limit is set by the limit parameter
	// but should not be reached.
	
	double dt = parameters->tmin; // running time step
	int current_ind = 0; // running index of computation
	int **counters_M = (int **)malloc(parameters->BH_number*sizeof(int *)); // positions in the f / g tables
	int **counters_param = (int **)malloc(parameters->BH_number*sizeof(int *)); // positions in the f / g tables
	double **estimates_M = (double **)malloc(parameters->BH_number*sizeof(double *)); // Euler estimates
	double **estimates_param = (double **)malloc(parameters->BH_number*sizeof(double *)); // Euler estimates
	life_times[0] = parameters->tmin;
	dts[0] = parameters->tmin;
	for(int i = 0;i<parameters->BH_number;i++){
		estimates_M[i] = (double *)malloc(parameters->param_number*sizeof(double));
		estimates_param[i] = (double *)malloc(parameters->param_number*sizeof(double));
		counters_M[i] = (int *)malloc(parameters->param_number*sizeof(int));
		counters_param[i] = (int *)malloc(parameters->param_number*sizeof(int));
		for(int j = 0;j<parameters->param_number;j++){
			life_masses[i][j][0] = init_masses[i]; // initialization of the masses
			life_params[i][j][0] = init_params[j]; // initialization of the spins
			counters_M[i][j] = parameters->nb_fM_masses-1; // masses are decreasing
			counters_param[i][j] = parameters->nb_fM_param-1; // spins are decreasing
		}
	}
	double ratio_M; // contains the relative evolution of the BH mass
	double ratio_param; // contains the relative evolution of the BH parameter
	int counter = 0; // contains the current total number of tried steps
	int interesting_BH = 0; // contains the running index in the sorted lifetimes
	int interesting_i; // contains the position of the interesting BH mass
	int interesting_j; // contains the position of the interesting BH spin
	int output_BH = 0; // decides whether we print the computation advancement
	parameters->nb_fin_times = 0;
	double limiting_mass;
	if(parameters->BH_remnant){
		limiting_mass = parameters->M_remnant;
	}
	else{
		limiting_mass = 1.1*Mp;
	}
	if(parameters->full_output){
		printf("\n\n\t computing:\t");
		fflush(stdout);
	}
	while(counter<parameters->limit-1 && interesting_BH < parameters->BH_number*parameters->param_number){ // limit of recursion
#if defined(_OPENMP)
#pragma omp parallel for
#endif
		for(int i = 0;i<parameters->BH_number;i++){
			for(int j = 0;j<parameters->param_number;j++){
				if(life_masses[i][j][current_ind] > limiting_mass){ // We continue to evolve this BH
					counters_M[i][j] = parameters->nb_fM_masses-1;
					counters_param[i][j] = parameters->nb_fM_param-1;
					while(counters_M[i][j] > 0 && life_masses[i][j][current_ind] < fM_masses[counters_M[i][j]]){ // we find the closest tabulated value
						counters_M[i][j]--;
					}
					while(counters_param[i][j] > 0 && life_params[i][j][current_ind] < fM_param[counters_param[i][j]]){ // we find the closest tabulated value
						counters_param[i][j]--;
					}
					estimates_M[i][j] = life_masses[i][j][current_ind] + dt*loss_rate_M(life_masses[i][j][current_ind],life_params[i][j][current_ind],fM_table,fM_masses,fM_param,counters_M[i][j],counters_param[i][j],parameters); // Euler estimation of the next mass
					estimates_param[i][j] = life_params[i][j][current_ind] + dt*loss_rate_param(life_masses[i][j][current_ind],life_params[i][j][current_ind],fM_table,gM_table,fM_masses,fM_param,counters_M[i][j],counters_param[i][j],parameters); // Euler estimation of the next parameter
				}
				else{ // end of the evolution for this BH
					if(rank[i][j] >= interesting_BH){
						interesting_BH++;
					}
					if(parameters->BH_remnant){
						estimates_M[i][j] = limiting_mass;
					}
					else{
						estimates_M[i][j] = 0.;
					}
					estimates_param[i][j] = 0.;
				}
			}
		}
		if(interesting_BH < parameters->BH_number*parameters->param_number){
			for(int i = 0;i<parameters->BH_number;i++){
				for(int j = 0;j<parameters->param_number;j++){
					if(rank[i][j] == interesting_BH){
						interesting_i = i;
						interesting_j = j;
					}
				}
			}
			ratio_M = estimates_M[interesting_i][interesting_j]/life_masses[interesting_i][interesting_j][current_ind]; // the relative mass variation
			if(life_params[interesting_i][interesting_j][current_ind] != 0.){
				ratio_param = estimates_param[interesting_i][interesting_j]/life_params[interesting_i][interesting_j][current_ind]; // the relative parameter variation
			}
			else{
				ratio_param = 0.;
			}
			if(ratio_M < 0.9 || (parameters->metric == 0 && (ratio_param != 0. && ((ratio_param < 1. && ratio_param < 0.9) || (ratio_param > 1. && ratio_param > 1.1))))){ // one of the variations is too high, we diminish the time step
				dt = dt/2.;
			}
			else if(ratio_M > 0.99 && (parameters-> metric != 0 || ratio_param == 0. || ((ratio_param <= 1. && ratio_param > 0.99) || (ratio_param > 1. && ratio_param < 1.01))) && dt/life_times[current_ind] < 1.){ // all the variations are too low and the time step is not too high compared to the time scale, we increase the time step
				dt = dt*2.;
			}
			else{ // all the variations are fine, we record the evolution and move forward
				current_ind++;
				life_times[current_ind] = life_times[current_ind-1] + dt;
				dts[current_ind] = dt;
				parameters->nb_fin_times = current_ind;
				for(int i = 0;i<parameters->BH_number;i++){
					for(int j = 0;j<parameters->param_number;j++){
						if(estimates_M[i][j] > Mp){
							life_masses[i][j][current_ind] = estimates_M[i][j];
							if(estimates_param[i][j] < 1.e-3){ // limit spin under which we consider that it is 0
								life_params[i][j][current_ind] = 0.;
							}
							else{
								life_params[i][j][current_ind] = estimates_param[i][j];
							}
						}
						else{
							life_masses[i][j][current_ind] = 0.;
							life_params[i][j][current_ind] = 0.;
						}
					}
				}
			}
			counter++;
			if(parameters->full_output && output_BH != interesting_BH){
				output_BH = interesting_BH;
				if(parameters->BH_number*parameters->param_number > 10 && output_BH%(parameters->BH_number*parameters->param_number/10) == 0){ // bug corrected : if there is less than 10 BHs we don't seek to show the current status (too fast anyway)
					printf("%i/%i ",output_BH,parameters->BH_number*parameters->param_number);
					fflush(stdout);
				}
			}
		}
	}
	if(parameters->full_output){
		printf("%i/%i\n\n",parameters->BH_number*parameters->param_number,parameters->BH_number*parameters->param_number);
		fflush(stdout);
	}
	free2D_int(counters_M,parameters->BH_number);
	free2D_int(counters_param,parameters->BH_number);
	free2D_double(estimates_M,parameters->BH_number);
	free2D_double(estimates_param,parameters->BH_number);
	if(counter == parameters->limit-1){
		printf("\n\t [life_evolution] : ERROR ITERATION LIMIT REACHED !\n");
		fflush(stdout);
	}
	return;
}

void write_life_evolutions(double ***life_masses,double ***life_params,double *life_times,double *dts,struct param *parameters){
	// This function writes the number of time steps and the masses and spins as
	// functions of time computed by the function life_evolution for
	// all initial BH masses and spins in the file 'life_evolutions.txt'.
	
	char file_name[500];
	sprintf(file_name,"./results/%s/life_evolutions.txt",parameters->destination_folder);
	FILE *file = fopen(file_name,"w+");
	if(!file){
		printf("\n\t [write_life_evolutions] : ERROR COULD NOT OPEN FILE '%s' !\n",file_name);
		fflush(stdout);
		exit(0);
	}
	rewind(file);
	fprintf(file,"Evolution of the BH masses and spins as functions of time.\n");
	fprintf(file,"Total number of time iterations: %i\n\n",parameters->nb_fin_times);
	
	for(int i = 0;i<parameters->BH_number;i++){
		fprintf(file,"%15s","t");
		for(int j = 0;j<parameters->param_number;j++){
			if(parameters->metric == 0){
				fprintf(file,"%5s%15s%15s","","M","a");
			}
			else if(parameters->metric == 1){
				fprintf(file,"%5s%15s%15s","","M","epsilon");
			}
		}
		fprintf(file,"\n");
		for(int k = 0;k<parameters->nb_fin_times;k++){
			fprintf(file,"%15.5e",life_times[k]/time_conversion);
			for(int j = 0;j<parameters->param_number;j++){
				fprintf(file,"%5s%15.5e%15.5e","",life_masses[i][j][k]/mass_conversion,life_params[i][j][k]);
			}
			fprintf(file,"\n");
		}
		fprintf(file,"\n\n");
	}
	fclose(file);
	
	sprintf(file_name,"./results/%s/dts.txt",parameters->destination_folder);
	file = fopen(file_name,"w+");
	if(!file){
		printf("\n\t [write_life_evolutions] : ERROR COULD NOT OPEN FILE '%s' !\n",file_name);
		fflush(stdout);
		exit(0);
	}
	rewind(file);
	fprintf(file,"Evolution of the integration timestep as a function of time.\n");
	fprintf(file,"%15s%15s\n","t","dt");
	for(int i = 0;i<parameters->nb_fin_times;i++){
		fprintf(file,"%15.5e%15.5e\n",life_times[i]/time_conversion,dts[i]/time_conversion);
	}
	fclose(file);
	
	return;
}
