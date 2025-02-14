// Program that computes PBH constraints from a set of BlackHawk results.
// Author: Jérémy Auffinger, j.auffinger@ipnl.in2p3.fr
// Last modification: 4 January 2021

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <dirent.h>
#include <string.h>
#include <unistd.h>

#define tab_len 200 // this is the default length of the tabulated values; if a segmentation fault occurs, just expand it

// conversion factors
#define kpc_to_cm 		3.085678e+21 // 1 kpc = XXX cm
#define mass_conversion 5.60958884e+23 // 1 g = XXX GeV
#define leng_conversion 5.06773058e+13 // 1 cm = XXX GeV

#define Mp		1.221e+19 // Planck mass in GeV				
#define G 		6.67430e-11 // SI units
#define c		2.99792458e+10 // CGS units cm.s^-1
#define pi		3.141592
#define GF		1.1663787e-5 // in GeV^-2
#define sin2W	0.213 // Weinberg angle squared (PDG)

// masses of the Standard Model particles (PDG 2017) in GeV
#define m_electron  0.5109989461e-3
#define m_neutrino  0. // simplifying hypothesis, considering that even the heaviest primordial black hole has a temperature >> m_neutrinos
#define m_muon  	105.6583745e-3
#define m_tau  		1.77686
#define m_up 		2.2e-3
#define m_down  	4.7e-3
#define m_strange  	96e-3
#define m_charm  	1.27
#define m_top  		173.21
#define m_bottom  	4.18
#define m_higgs  	125.03
#define m_photon   	0.
#define m_gluon   	0.2 // effective mass
#define m_wpm  		80.403
#define m_z0  		91.1876
#define m_graviton  0.

#define m_pi0		134.9766e-3
#define m_pipm  	139.57018e-3
#define m_K0L   	497.7e-3
#define m_Kpm   	493.7e-3
#define m_proton    938.272e-3
#define m_neutron   939.5654e-3

int nb_params = 2; // number of parameters needed for one run


void free1D_double(double *array){
	// This function frees the memory allocated to a 1D double array.
	
	if(array!=NULL) free(array);
	array = NULL;
	return;
}

void free2D_double(double **array,int l_1stD){
	// This function frees the memory allocated to a 2D double array
	// of first dimension l_1stD.
	
	if(array==NULL) return;

	for(int i = 0;i<l_1stD;i++){
		if(array[i]!=NULL) free(array[i]);
		array[i] = NULL;
	}
	if(array!=NULL) free(array);
	array = NULL;
	return;
}

double gaussian(double x,double mean, double standdev){
	// this routine computes the gaussian function of mean "mean" and relative standard deviation "standdev" at point x 
	
	return 1./(sqrt(2.*pi)*standdev*mean)*exp(-pow(x - mean,2.)/(2.*pow(standdev*mean,2.)));
}

struct param{ // this structure will contain all the parameters of a BlackHawk run
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
	double M_star; // modified Planck mass in extra dimension models
	int spectrum_choice; // chooses the initial black holes mass distribution. 0: Dirac, 1: lognormal distribution (mass), 11: lognormal distribution (number), 2: power-law distribution, 
							// 3: critical collapse distribution, 4: peak-theory distribution, 5:uniform, -1: User's distribution
	int spectrum_choice_param; // chooses the initial black holes spin distribution. 0: Dirac, 1: uniform, 2: gaussian.
	double amplitude; // amplitude of the distribution
	double stand_dev; // standard deviation of the mass distribution
	double crit_mass; // characteristic mass
	double eqstate; // equation of state parameter P = w*rho
	double stand_dev_param; // standard deviation of the parameter distribution
	double mean_param; // mean of the parameter distribution
	char table[32]; // name of the User's mass/spin distribution file
	
	int tmin_manual; // 0: automatic tmin for Dirac distribution, 1: manual tmin
	double tmin; // time of black hole formation
	int limit; // maximum iteration limit when computing the life evolution of black holes
	int nb_fin_times; // contains the number of integration times
	int BH_remnant; // 0: total evaporation after Planck mass, 1: stop the evaporation at M_relic
	double M_remnant; // mass of the BH remnant
	
	int E_number; // number of initial energies
	double Emin; // minimum initial energy
	double Emax; // maximum initial energy
	int particle_number; // number of particle species to be emitted
	int grav; // 0: no gravitons emitted, 1: gravitons emitted
	int add_DM; // 0: no DM, 1: one primary DM particle
	double m_DM; // DM mass in GeV
	double spin_DM; // spin of DM among 0., 1., 2., 0.5, 1.5
	double dof_DM; // number of DM degrees of freedom
	
	int primary_only; // if set to 1, the code will skip the hadronization part, if set to 0, the secondary spectra are computed
	int hadronization_choice; // choice of hadronization tables (0: PYTHIA, 1: HERWIG, 2: PYTHIA NEW, 3: HAZMA)
		
	double Mmin_fM; // minimum BH mass in the f / g tables
	double Mmax_fM; // maximum BH mass in the f / g tables
	int nb_fM_masses; // number of black hole masses in the f / g tables
	int nb_fM_param; // number of black hole spins in the f(M,a) / g(M,a) tables, or BH espilon in the f(M,epsilon) tables
	double param_min; // minimum BH secondary parameter in the f,g and gamma tables
	double param_max; // maximum BH secondary parameter in the f,g and gamma tables
	int nb_gamma_param; // number of particle spins in the gamma(E,M,a) tables, or BH epsilon in the gamma(E,M,epsilon,a0) tables
	int nb_gamma_x; // number of the BH masses in the gamma tables
	int nb_gamma_spins; // number of spins (0, 1, 2, 1/2, [3/2]) tabulated
	int nb_gamma_fits; // number of fitting coefficients in the gamma fits tables
	double Emin_hadro; // minimum initial energy in the hadronization tables
	double Emax_hadro; // maximum initial energy in the hadronization tables
	int nb_init_en; // number of initial energies in the hadronization tables
	int nb_fin_en; // number of final energies in the hadronization tables
	int nb_init_part; // number of initial particles in the hadronization tables
	int nb_fin_part; // number of final particles in the hadronization tables
};

struct Isatis{ // this structure contains the parameters used in the statistical constraint method
	char output[64]; // name of the output file
	char path[256]; // path to the BlackHawk version
	int sessions; // 0: instantaneous only, 1: time-dependent only, 2: both
	double local_DM_GeV; // in GeV/cm^3 https://arxiv.org/abs/1212.3670
	double local_DM; // in g/cm^3 
	double global_DM; // in g/cm^3 (Planck results 2018 Omega_DM*h**2 = 0.12 and h = 0.674)
	double r_0; // distance Sun-GC in Kpc->cm
	int profile; // density profile of the galaxy, 0: generalized NFW, 1: Einasto
	double rho_c; // characteristic halo density in g/cm^3
	double r_c; // characteristic halo radius in kpc
	double gamma; // density profile inner slope
	double t_eq; // matter-radiation equality in s
	double t_today; // age of the universe in s
	double t_CMB; // CMB last scattering surface in s
	int domination; // redshift history, 0: standard, 1: radiation domination, 2: matter domination
	int energy_method; // statistical sub-method used: 1: full energy range, 2: binned energy range, see arXiv:XXXX.XXXXX
	int background_type; // 0: gal+egal arXiv:2010.04797, 1: gal+egal arXiv:2110.03333, 2: gal arXiv:1101.1381 and egal arXiv:1906.04750
	int confidence_level; // number of error bars above the data
	double S_to_N; // signal to noise ratio for prospective instruments
};

struct type_1{ // this structure contains all the relevant parameters for type 2 constraints (prospective X/gamma-ray detectors)
	char name[64]; // name of the constraint
	int type; // 0: galactic, 1: extragalactic
	int nb_ener; // number of energies in the flux(E) file
	int format; // number of columns in the flux(E) file
	double l_min; // minimum latitude
	double l_max; // maximum latitude
	double b_min; // minimum longitude
	double b_max; // maximum longitude
};

struct type_2{ // this structure contains all the relevant parameters for type 2 constraints (prospective X/gamma-ray detectors)
	char name[64]; // name of the constraint
	int nb_ener; // number of energies in the Aeff(E) file
	double l_min; // minimum latitude
	double l_max; // maximum latitude
	double b_min; // minimum longitude
	double b_max; // maximum longitude
	int nb_bckg_gal; // number of galactic background energies
	int nb_bckg_egal; // number of extragalactic background energies
	double T_obs; // time of observation
};

void read_Isatis(struct Isatis *isatis,char *name){
	// This function reads the Isatis parameters in the file name
	// and fills the Isatis structure.
	
	FILE *param_file;
	param_file = fopen(name,"r");
	if(!param_file){
		printf("\n\t\t [read_Isatis] : ERROR COULD NOT OPEN FILE '%s' !\n",name);
		fflush(stdout);
		exit(0);
	}
	rewind(param_file);
	
	char dummy[64]; // contains irrelevant characters of the parameters file (2 for each parameter, 1st is its name and 2nd is =)

	while(EOF != fscanf(param_file,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(param_file,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		if(!strcasecmp(dummy,"output"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%s",isatis->output);
		}
		else if(!strcasecmp(dummy,"path"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%s",isatis->path);
		}
		else if(!strcasecmp(dummy,"sessions"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(isatis->sessions));
		}
		else if(!strcasecmp(dummy,"local_DM_GeV"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(isatis->local_DM_GeV));
			isatis->local_DM = isatis->local_DM_GeV/mass_conversion; // in g.cm^-3
		}
		else if(!strcasecmp(dummy,"global_DM"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(isatis->global_DM));
		}
		else if(!strcasecmp(dummy,"r_0"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(isatis->r_0));
			isatis->r_0 = isatis->r_0*kpc_to_cm; // conversion from kpc to cm
		}
		else if(!strcasecmp(dummy,"profile"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(isatis->profile));
		}
		else if(!strcasecmp(dummy,"rho_c_halo"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(isatis->rho_c)); // already in g/cm^3
		}
		else if(!strcasecmp(dummy,"r_c_halo"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(isatis->r_c));
			isatis->r_c = isatis->r_c*kpc_to_cm; // conversion from kpc to cm
		}
		else if(!strcasecmp(dummy,"gamma_halo"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(isatis->gamma));
		}
		else if(!strcasecmp(dummy,"t_eq"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(isatis->t_eq));
		}
		else if(!strcasecmp(dummy,"t_CMB"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(isatis->t_CMB));
		}
		else if(!strcasecmp(dummy,"t_today"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(isatis->t_today));
		}
		else if(!strcasecmp(dummy,"domination"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(isatis->domination));
		}
		else if(!strcasecmp(dummy,"energy_method"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(isatis->energy_method));
		}
		else if(!strcasecmp(dummy,"background_type"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(isatis->background_type));
		}
		else if(!strcasecmp(dummy,"confidence_level"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(isatis->confidence_level));
		}
		else if(!strcasecmp(dummy,"signal_to_noise_ratio"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(isatis->S_to_N));
		}
	}
	
	return;
}

void read_params(struct param *parameters,char *name){
	// This function reads the parameters file of the run and fills
	// the parameters structure.
	
	FILE *param_file;
	param_file = fopen(name,"r");
	if(!param_file){
		printf("\n\t\t [read_params] : ERROR COULD NOT OPEN FILE '%s' !\n",name);
		fflush(stdout);
		exit(0);
	}
	rewind(param_file);
	
	char dummy[64]; // contains irrelevant characters of the parameters file (2 for each parameter, 1st is its name and 2nd is =)
	
	double Emin_hadro_PYTHIA,Emax_hadro_PYTHIA;
	int nb_fin_en_PYTHIA,nb_init_en_PYTHIA,nb_init_part_PYTHIA,nb_fin_part_PYTHIA;
	double Emin_hadro_HERWIG,Emax_hadro_HERWIG;
	int nb_fin_en_HERWIG,nb_init_en_HERWIG,nb_init_part_HERWIG,nb_fin_part_HERWIG;
	double Emin_hadro_PYTHIA_new,Emax_hadro_PYTHIA_new;
	int nb_fin_en_PYTHIA_new,nb_init_en_PYTHIA_new,nb_init_part_PYTHIA_new,nb_fin_part_PYTHIA_new;
	
	Emin_hadro_PYTHIA=Emax_hadro_PYTHIA=Emin_hadro_HERWIG=Emax_hadro_HERWIG=Emin_hadro_PYTHIA_new=Emax_hadro_PYTHIA_new=0.;
	nb_fin_en_PYTHIA=nb_init_en_PYTHIA=nb_init_part_PYTHIA=nb_fin_part_PYTHIA=nb_fin_en_HERWIG=nb_init_en_HERWIG=nb_init_part_HERWIG=nb_fin_part_HERWIG=nb_fin_en_PYTHIA_new=nb_init_en_PYTHIA_new=nb_init_part_PYTHIA_new=nb_fin_part_PYTHIA_new=0;	
	
	double amplitude_lognormal,amplitude_lognormal2,stand_dev_lognormal,crit_mass_lognormal,amplitude_powerlaw,eqstate_powerlaw,amplitude_critical_collapse,crit_mass_critical_collapse,amplitude_uniform;
	double stand_dev_param_gaussian,mean_param_gaussian;
	
	amplitude_lognormal=amplitude_lognormal2=stand_dev_lognormal=crit_mass_lognormal=amplitude_powerlaw=eqstate_powerlaw=amplitude_critical_collapse=crit_mass_critical_collapse=amplitude_uniform=0.;
	stand_dev_param_gaussian=mean_param_gaussian=0.;
	
	while(EOF != fscanf(param_file,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(param_file,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		if(!strcasecmp(dummy,"destination_folder"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%s",parameters->destination_folder);
		}
		else if(!strcasecmp(dummy,"full_output"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->full_output));
		}
		else if(!strcasecmp(dummy,"interpolation_method"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->interpolation_method));
		}
		else if(!strcasecmp(dummy,"metric"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->metric));
		}
		else if(!strcasecmp(dummy,"BH_number"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->BH_number));
		}
		else if(!strcasecmp(dummy,"Mmin"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(parameters->Mmin));
		}
		else if(!strcasecmp(dummy,"Mmax"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(parameters->Mmax));
		}
		else if(!strcasecmp(dummy,"param_number"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->param_number));
		}
		else if(!strcasecmp(dummy,"amin"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(parameters->amin));
		}
		else if(!strcasecmp(dummy,"amax"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(parameters->amax));
		}
		else if(!strcasecmp(dummy,"Qmin"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(parameters->Qmin));
		}
		else if(!strcasecmp(dummy,"Qmax"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(parameters->Qmax));
		}
		else if(!strcasecmp(dummy,"epsilon_LQG"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(parameters->epsilon_LQG));
		}
		else if(!strcasecmp(dummy,"a0_LQG"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(parameters->a0_LQG));
		}
		else if(!strcasecmp(dummy,"n"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(parameters->n));
		}
		else if(!strcasecmp(dummy,"spectrum_choice"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->spectrum_choice));
		}
		else if(!strcasecmp(dummy,"spectrum_choice_param"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->spectrum_choice_param));
		}
		else if(!strcasecmp(dummy,"amplitude_lognormal"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(amplitude_lognormal));

		}
		else if(!strcasecmp(dummy,"amplitude_lognormal2"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(amplitude_lognormal2));

		}
		else if(!strcasecmp(dummy,"stand_dev_lognormal"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(stand_dev_lognormal));	
		}
		else if(!strcasecmp(dummy,"crit_mass_lognormal"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(crit_mass_lognormal));
			
		}
		else if(!strcasecmp(dummy,"amplitude_powerlaw"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(amplitude_powerlaw));

		}
		else if(!strcasecmp(dummy,"eqstate_powerlaw"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(eqstate_powerlaw));
		}
		else if(!strcasecmp(dummy,"amplitude_critical_collapse"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(amplitude_critical_collapse));
		}
		else if(!strcasecmp(dummy,"crit_mass_critical_collapse"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&crit_mass_critical_collapse);
		}
		else if(!strcasecmp(dummy,"amplitude_uniform"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&amplitude_uniform);
		}
		else if(!strcasecmp(dummy,"stand_dev_param_gaussian"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&stand_dev_param_gaussian);
		}
		else if(!strcasecmp(dummy,"mean_param_gaussian"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&mean_param_gaussian);
		}
		else if(!strcasecmp(dummy,"table"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%s",parameters->table);
		}
		else if(!strcasecmp(dummy,"tmin_manual"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->tmin_manual));
		}
		else if(!strcasecmp(dummy,"tmin"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(parameters->tmin));
			if(!(parameters->tmin_manual)){
				double gamma_BH = 0.2; // fraction of particle horizon mass that collapses into a BH during radiation-dominated phase (see Arbey et al. 2012.09867)
				parameters->tmin = parameters->Mmin*G/(gamma_BH*pow(c,3.));
			}
		}
		else if(!strcasecmp(dummy,"limit"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->limit));
			parameters->limit = parameters->limit*parameters->BH_number*parameters->param_number;
		}
		else if(!strcasecmp(dummy,"BH_remnant"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->BH_remnant));
		}
		else if(!strcasecmp(dummy,"M_remnant"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(parameters->M_remnant));
		}
		else if(!strcasecmp(dummy,"E_number"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->E_number));
		}
		else if(!strcasecmp(dummy,"Emin"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(parameters->Emin));
		}
		else if(!strcasecmp(dummy,"Emax"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(parameters->Emax));
		}
		else if(!strcasecmp(dummy,"grav"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->grav));
		}
		else if(!strcasecmp(dummy,"add_DM"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->add_DM));
		}
		else if(!strcasecmp(dummy,"m_DM"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(parameters->m_DM));
		}
		else if(!strcasecmp(dummy,"spin_DM"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(parameters->spin_DM));
		}
		else if(!strcasecmp(dummy,"dof_DM"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(parameters->dof_DM));
		}
		else if(!strcasecmp(dummy,"primary_only"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->primary_only));
		}
		else if(!strcasecmp(dummy,"hadronization_choice"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->hadronization_choice));
		}
	}
	fclose(param_file);
	
	if(parameters->metric == 1){
		parameters->param_number = 1; // epsilon is the same for all BHs
		parameters->spectrum_choice_param = 0;
	}
	
	if(parameters->metric == 3){
		parameters->param_number = 1; // n is the same for all BHs
		parameters->spectrum_choice_param = 0;
		parameters->M_star = 1.;
	}
	
	if(parameters->spectrum_choice==1)
	{
		parameters->amplitude=amplitude_lognormal;
		parameters->stand_dev=stand_dev_lognormal;
		parameters->crit_mass=crit_mass_lognormal;
	}
	if(parameters->spectrum_choice==11)
	{
		parameters->amplitude=amplitude_lognormal2;
		parameters->stand_dev=stand_dev_lognormal;
		parameters->crit_mass=crit_mass_lognormal;
	}
	else if(parameters->spectrum_choice==2)
	{
		parameters->amplitude=amplitude_powerlaw;
		parameters->eqstate=eqstate_powerlaw;
		
		double gamma = 2.*parameters->eqstate/(1.+parameters->eqstate);	
	}
	else if(parameters->spectrum_choice==3)
	{
		parameters->amplitude=amplitude_critical_collapse;
		parameters->crit_mass=crit_mass_critical_collapse;
	}
	
	if(parameters->spectrum_choice_param == 2){
		parameters->stand_dev_param = stand_dev_param_gaussian;
		parameters->mean_param = mean_param_gaussian;
	}
				
	if(parameters->spectrum_choice == 0 && parameters->BH_number > 1)
	{
		parameters->BH_number = 1;
	}
	if(!(parameters->spectrum_choice == -1) && (parameters->metric == 0 || parameters->metric == 2) && parameters->spectrum_choice_param == 0 && parameters->param_number > 1)
	{
		parameters->param_number = 1;
	}
	
	if(parameters->hadronization_choice == 3){
		if(parameters->add_DM == 0){
			parameters->add_DM = 1;
			parameters->m_DM = m_pi0; // we have m_pi0 < m_pipm and we put back the constraint of E>m_pipm manually for the FSR 
			parameters->spin_DM = 0.;
			parameters->dof_DM = 1.;
		}
	}
	
	parameters->nb_init_part = 15+parameters->grav+parameters->add_DM;
	
	switch(parameters->hadronization_choice){
		case 0:{ // PYTHIA tables
			parameters->nb_fin_en = 500;
			parameters->nb_fin_part = 11;
			break;
		}
		case 1:{ // HERWIG tables
			parameters->nb_fin_en = 100;
			parameters->nb_fin_part = 11;
			break;
		}
		case 2:{ // PYTHIA new tables
			parameters->nb_fin_en = 500;
			parameters->nb_fin_part = 6;
			break;
		}
		case 3:{ // Hazma tables
			parameters->nb_fin_en = 500;
			parameters->nb_fin_part = 2;
			break;
		}
		default:{
			exit(0);
			break;
		}
	}
	
	return;
}

void read_data(char *path,double *BH_masses,double *BH_params,double *BH_spectrum,double *init_ener,double *fin_ener,double *inst_primary_spectra,double *inst_secondary_spectra,double *dts,struct param *parameters,struct Isatis *isatis){
	// This routine reads the spectra data in the BlackHawk results
	// and the time intervals.
	
	char dummy[64];
	
	printf("BH spectrum...");
	fflush(stdout);
	
	char *spectrum_name = (char *)malloc(256*sizeof(char));
	sprintf(spectrum_name,"%s/results/%s/BH_spectrum.txt",isatis->path,path);
	FILE *spectrum_file = fopen(spectrum_name,"r+");
	if(!spectrum_file){
		printf("\n\t\t  [read_data] : ERROR could not open file %s !");
		fflush(stdout);
		exit(0);
	}
	rewind(spectrum_file);
	
	for(int i = 0;i<14;i++){ // spectrum file header
		fscanf(spectrum_file,"%s",dummy);
	}
	for(int i = 0;i<parameters->param_number;i++){
		fscanf(spectrum_file,"%lf",&(BH_params[i]));
	}
	for(int i = 0;i<parameters->BH_number;i++){
		fscanf(spectrum_file,"%lf",&(BH_masses[i]));
		for(int j = 0;j<parameters->param_number;j++){
			fscanf(spectrum_file,"%lf",&(BH_spectrum[i*parameters->param_number+j]));
		}
	}
	
	if(isatis->sessions == 0 || isatis->sessions == 2){
		
		printf("instantaneous spectra...");
		fflush(stdout);
		
		char *primary_name = (char *)malloc(256*sizeof(char));
		char *secondary_name = (char *)malloc(256*sizeof(char));
		sprintf(primary_name,"%s/results/%s/instantaneous_primary_spectra.txt",isatis->path,path);
		sprintf(secondary_name,"%s/results/%s/instantaneous_secondary_spectra.txt",isatis->path,path);
		
		FILE *primary_file = fopen(primary_name,"r+");
		if(!primary_file){
			printf("\n\t\t  [read_data] : ERROR could not open file %s !",primary_file);
			fflush(stdout);
			exit(0);
		}
		rewind(primary_file);
		for(int i = 0;i<8+parameters->nb_init_part;i++){
			fscanf(primary_file,"%s",dummy);
		}
		for(int i = 0;i<parameters->E_number;i++){
			fscanf(primary_file,"%lf",&(init_ener[i]));
			for(int j = 0;j<parameters->nb_init_part;j++){
				fscanf(primary_file,"%lf",&(inst_primary_spectra[i*parameters->nb_init_part+j]));
			}
		}
		fclose(primary_file);
		
		FILE *secondary_file = fopen(secondary_name,"r+");
		if(!secondary_file){
			printf("\n\t\t [read_data] : ERROR could not open file %s !",secondary_file);
			fflush(stdout);
			exit(0);
		}
		rewind(secondary_file);
		for(int i = 0;i<8+parameters->nb_fin_part;i++){ // file header
			fscanf(secondary_file,"%s",dummy);
		}
		for(int i = 0;i<parameters->nb_fin_en;i++){
			fscanf(secondary_file,"%lf",&(fin_ener[i]));
			for(int j = 0;j<parameters->nb_fin_part;j++){
				fscanf(secondary_file,"%lf",&(inst_secondary_spectra[i*parameters->nb_fin_part+j]));
			}
		}
		fclose(secondary_file);
	}
	if(isatis->sessions == 1 || isatis->sessions == 2){
		
		printf("time-dependent spectra...");
		fflush(stdout);
		
		char *dts_name = (char *)malloc(256*sizeof(char));
		sprintf(dts_name,"%s/results/%s/dts.txt",isatis->path,path);
		FILE *dts_file = fopen(dts_name,"r+");
		if(!dts_file){
			printf("\n\t\t [read_data] : ERROR could not open file %s !",dts_name);
			fflush(stdout);
			exit(0);
		}
		rewind(dts_file);
		for(int i = 0;i<12;i++){
			fscanf(dts_file,"%s",dummy);
		}
		for(int i = 0;i<parameters->nb_fin_times;i++){
			fscanf(dts_file,"%lf",&(dts[i*2 + 0]));
			fscanf(dts_file,"%lf",&(dts[i*2 + 1]));
		}
		fclose(dts_file);
	}
	
	printf(" ");
	fflush(stdout);
	return;
}

double local_amplitude_DM(double *BH_masses,double *BH_params,double *BH_spectrum,struct param *parameters,struct Isatis *isatis){
	// This routine computes the local amplitude of the BH distribution (for local normalization)
	
	double result = 0.;
	for(int i = 0;i<parameters->BH_number;i++){
		for(int j = 0;j<parameters->param_number;j++){
			result += BH_masses[i]/isatis->local_DM*BH_spectrum[i*parameters->param_number + j]; // this is dimensionless
		}
	}
	
	return result;
}

double total_amplitude_DM(double *BH_masses,double *BH_params,double *BH_spectrum,struct param *parameters){
	// this routine computes the total amplitude of the BH distribution (for galactic normalization)
	
	double result = 0.;
	for(int i = 0;i<parameters->BH_number;i++){
		for(int j = 0;j<parameters->param_number;j++){
			result += BH_masses[i]*BH_spectrum[i*parameters->param_number + j]; // in g^-1.cm^3
		}
	}
	
	return result;
}

double global_amplitude_DM(double *BH_masses,double *BH_params,double *BH_spectrum,struct param *parameters,struct Isatis *isatis){
	// this routine computes the cosmological amplitude of the PBH distribution (for cosmological normalization)
	
	double result = 0.;
	for(int i = 0;i<parameters->BH_number;i++){
		for(int j = 0;j<parameters->param_number;j++){
			result += BH_masses[i]/isatis->global_DM*BH_spectrum[i*parameters->param_number + j]; // in cm^3
		}
	}
	
	return result;
}

double galactic_profile(double r,struct Isatis *isatis){
	// this routine computes the density of DM depending on the galactic coordinates in various profiles
	
	double rho;
	switch(isatis->profile){
		case 0:{ // generalized Navarro-Frenk-White profile (arXiv:9611107 and arXiv:1906.06133)
			// "convenient"parameters of arXiv:1102.4340: rho_c = 0.0125 M_sol/pc^3 = 8.5e-25 g/cm^3, r_c = 17 kpc, r_0 = 8.5 kpc, gamma = 1
			rho = isatis->rho_c*pow(isatis->r_c/r,isatis->gamma)*pow(1 + r/isatis->r_c,isatis->gamma-3);
			break;
		}
		case 1:{ // Einasto profile (original Einasto paper and arXiv:1906.06133)
			rho = isatis->rho_c*exp(-2./isatis->gamma*(pow(r/isatis->r_c,isatis->gamma) - 1.));
			break;
		}
		default:{
			printf("\n\t\t [galactic_profile] : ERROR wrong galactic density profile %i !",isatis->profile);
			fflush(stdout);
			exit(0);
			break;
		}
	}
	return rho;
}

double J_D(double l_min,double l_max,double b_min,double b_max,struct Isatis *isatis){
	// this routine computes the density factor J_D of DM as in eq. (5) of arXiv:2110.03333
	
	int nb_angles = 100;
	int nb_radii = 100;
	double r_max = 200.*kpc_to_cm; // maximum radius of the Galactic halo in cm
	double s_max;
	double result = 0.;
	double *l = (double *)malloc(nb_angles*sizeof(double));
	double *b = (double *)malloc(nb_angles*sizeof(double));
	double *s = (double *)malloc(nb_radii*sizeof(double));
	double metric;
	double r;
	for(int i = 0;i<nb_angles;i++){
		l[i] = l_min + (l_max - l_min)/(nb_angles - 1)*i;
		b[i] = b_min + (b_max - b_min)/(nb_angles - 1)*i;
	}
	for(int i = 0;i<nb_angles-1;i++){ // integral over l
		for(int j = 0;j<nb_angles-1;j++){ // integral over b
			s_max = isatis->r_0*cos(l[i])*cos(b[j]) + sqrt(pow(r_max,2.) - pow(isatis->r_0,2.)*(1. - pow(cos(l[i])*cos(b[j]),2.)));
			for(int k = 0;k<nb_radii;k++){
				s[k] = 0. + (s_max - 0.)/(nb_radii - 1)*k;
			}
			for(int k = 0;k<nb_radii-1;k++){ // integral over s(r(l,b))
				metric = fabs(cos(b[i]))*(l[i+1] - l[i])*(b[j+1] - b[j])*(s[k+1] - s[k]);
				r = sqrt(pow(isatis->r_0,2.) + pow(s[k],2.) - 2.*isatis->r_0*s[k]*cos(l[i])*cos(b[j]));
				result += metric*galactic_profile(r,isatis);
			}
		}
	}
	double Delta = 0.; // the field of view of the telescope
	for(int i = 0;i<nb_angles-1;i++){
		for(int j = 0;j<nb_angles-1;j++){
			Delta += fabs(cos(b[i]))*(l[i+1] - l[i])*(b[j+1] - b[j]);
		}
	}
	
	return result/Delta; // in g.cm^(-2)
}

void galactic_contribution_type1(double *galactic,double *inst_secondary_spectra,double A,struct type_1 *type1,struct param *parameters,struct Isatis *isatis){
	// this routine computes the galactic contribution from some spectrum to the local flux
	
	double density_factor = J_D(type1->l_min,type1->l_max,type1->b_min,type1->b_max,isatis);
	for(int i = 0;i<parameters->nb_fin_en;i++){ // following arXiv:2110.03333
		galactic[i] = 1./(4.*pi*A)*density_factor*inst_secondary_spectra[i*parameters->nb_fin_part + 0]; // 0: photons
	}
	
	return;
}

void galactic_contribution_type2(double *galactic,double *inst_secondary_spectra,double A,struct type_2 *type2,struct param *parameters,struct Isatis *isatis){
	// this routine computes the galactic contribution from some spectrum to the local flux
	
	double density_factor = J_D(type2->l_min,type2->l_max,type2->b_min,type2->b_max,isatis);
	for(int i = 0;i<parameters->nb_fin_en;i++){ // following arXiv:2110.03333
		galactic[i] = 1./(4.*pi*A)*density_factor*inst_secondary_spectra[i*parameters->nb_fin_part + 0]; // 0: photons
	}
	
	return;
}

double redshift(double t,double t_evap,double t_eq,double t_today,int domination,int today){
	switch(domination){
		case 0:{ // standard redshift history
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
		case 1:{ // radiation domination
			return pow(t_evap/t,1./2.) - 1.;
			break;
		}
		case 2:{ // matter domination
			return pow(t_evap/t,2./3.) - 1.;
			break;
		}
		default:{
			printf("[redshift] : ERROR wrong redshift history choice %i !\n",domination);
			printf("[main] : end of execution");
			fflush(stdout);
			exit(0);
		}
	}
}

void extragalactic_contribution(double *extragalactic,double *stacked_ener,int CMB,int today,int nb_ener,double *fin_ener,double *dts,char *path,double A,struct param *parameters,struct Isatis *isatis){
	// this routine computes the extragalactic contribution of spectrum flux stacked from some initial time until today
	// it reads the secondary spectrum line by line to avoid extensive storage in RAM
	
	if(isatis->sessions == 0){
		printf("\n\t\t [extragalactic_contribution] : WARNING the extragalactic contribution will not be computed because sessions = %i !",isatis->sessions);
		fflush(stdout);
	}
	else{
		char dummy[64];
		char *spectrum_name = (char *)malloc(256*sizeof(char));
		sprintf(spectrum_name,"%s/results/%s/photon_secondary_spectrum.txt",isatis->path,path); // secondary photons
		FILE *spectrum_file = fopen(spectrum_name,"r+");
		if(!(spectrum_file)){
			printf("\n\t\t [extragalactic_contribution] : ERROR could not open file %s !",spectrum_name);
			fflush(stdout);
			exit(0);
		}
		rewind(spectrum_file);
		
		for(int i = 0;i<9+parameters->nb_fin_en;i++){
			fscanf(spectrum_file,"%s",dummy);
		}
		
		double t_evap = dts[(parameters->nb_fin_times-1)*2 + 0];
		
		double *init_spectrum = (double *)malloc(parameters->nb_fin_en*sizeof(double));
		double *interm_spectrum = (double *)malloc(nb_ener*sizeof(double));
		int counter;
		int too_low;
		double z;
		double t_today = isatis->t_today;
		double t_eq = isatis->t_eq;
		double t_CMB = isatis->t_CMB;
		int domination = isatis->domination;
		for(int i = 0;i<parameters->nb_fin_times;i++){
			fscanf(spectrum_file,"%s",dummy);
			for(int j = 0;j<parameters->nb_fin_en;j++){ // we read the time line of interest in the BlackHawk spectrum
				fscanf(spectrum_file,"%lf",&(init_spectrum[j]));
			}
			if((CMB == 0 || (CMB == 1 && dts[i*2 + 0] > t_CMB)) && (today == 0 || (today == 1 && dts[i*2 + 0] < t_today))){ // if this time line is of interest, we compute the redshifted spectrum and ad it to the stacked spectrum
				z = redshift(dts[i*2 + 0],t_evap,t_eq,t_today,domination,today);
				for(int j = 0;j<nb_ener;j++){
					counter = 0;
					too_low = 1;
					while(counter < parameters->nb_fin_en && fin_ener[counter]/(1.+z) <= stacked_ener[j]){
						too_low = 0;
						counter++;
					}
					if(too_low || counter == parameters->nb_fin_en){
						interm_spectrum[j] = 0.;
					}
					else{
						interm_spectrum[j] = fabs(init_spectrum[counter-1] + (init_spectrum[counter] - init_spectrum[counter-1])/(fin_ener[counter] - fin_ener[counter-1])*(stacked_ener[j] - fin_ener[counter-1]/(1.+z)));
					}
					extragalactic[j] += (1.+z)*interm_spectrum[j]*dts[i*2 + 1];
				}
			}
		}
		for(int i = 0;i<nb_ener;i++){
			extragalactic[i] = extragalactic[i]*c/A/(4.*pi);
		}
		
		free1D_double(init_spectrum);
		free1D_double(interm_spectrum);
	}
	
	return;
}

void background_type2(double *bckg_gal_ener,double *background_gal,double *bckg_egal_ener,double *background_egal,struct type_2 *type2,struct Isatis *isatis){
	// this routine computes the galactic and extragalactic backgrounds for type2 constraints
	
	char background_name[256];
	FILE *background_file;
	char dummy[64];
	switch(isatis->background_type){
		case 0:{ // see arXiv:1703.02546
			sprintf(background_name,"./constraints/photons/background_gal_egal_1703.02546.txt"); // in this case we set the "total background" in the galactic part and keep extragalactic to 0
			background_file = fopen(background_name,"r+");
			if(!background_file){
				printf("\n\t\t\t[constraint_type2] : ERROR could not open file %s !",background_name);
				fflush(stdout);
				exit(0);
			}
			rewind(background_file);
			if(strncasecmp("#",dummy,1)){
				while((EOF!=fscanf(background_file,"%c",dummy)) && (strncasecmp("\n",dummy,1))){
					//printf("%s",dummy);
					//fflush(stdout);
				}
			}
			fscanf(background_file,"%s",dummy);
			fscanf(background_file,"%s",dummy);
			fscanf(background_file,"%i",&(type2->nb_bckg_gal));
			type2->nb_bckg_egal = type2->nb_bckg_gal;
			fscanf(background_file,"%s",dummy);
			fscanf(background_file,"%s",dummy);
			
			for(int i = 0;i<type2->nb_bckg_gal;i++){
				fscanf(background_file,"%lf",&(bckg_gal_ener[i]));
				bckg_egal_ener[i] = bckg_gal_ener[i];
				fscanf(background_file,"%lf",&(background_gal[i]));
				background_gal[i] = background_gal[i]/pow(bckg_gal_ener[i],2.); // E².flux in GeV^2.cm^-2.s^-1.sr^-1.GeV^-1 to flux in cm^-2.s^-1.sr^-1.GeV^-1
				background_egal[i] = 0.;
			}
			fclose(background_file);
			break;
		}
		case 1:{ // see arXiv:2110.03333 
			double E0 = 35.6966e-6; // conversion from keV to GeV
			double n1 = 1.4199;
			double n2 = 2.8956;
			double A_eg = 64.2e+3; // conversion from MeV^-1 to GeV^-1
			sprintf(background_name,"./constraints/photons/background_gal_egal_2110.03333.txt");
			background_file = fopen(background_name,"r+");
			if(!background_file){
				printf("\n\t\t\t[constraint_type2] : ERROR could not open file %s !",background_name);
				fflush(stdout);
				exit(0);
			}
			rewind(background_file);
			if(strncasecmp("#",dummy,1)){
				while((EOF!=fscanf(background_file,"%c",dummy)) && (strncasecmp("\n",dummy,1))){
					//printf("%s",dummy);
					//fflush(stdout);
				}
			}
			fscanf(background_file,"%s",dummy);
			fscanf(background_file,"%s",dummy);
			fscanf(background_file,"%i",&(type2->nb_bckg_gal));
			type2->nb_bckg_egal = type2->nb_bckg_gal;
			fscanf(background_file,"%s",dummy);
			fscanf(background_file,"%s",dummy);
			for(int i = 0;i<type2->nb_bckg_gal;i++){
				fscanf(background_file,"%lf",&(bckg_gal_ener[i]));
				bckg_egal_ener[i] = bckg_gal_ener[i];
				fscanf(background_file,"%lf",&(background_gal[i]));
				background_egal[i] = A_eg/(pow(bckg_egal_ener[i]/E0,n1) + pow(bckg_egal_ener[i]/E0,n2)); // see eq. (8) of arXiv:2110.03333
			}
			fclose(background_file);
			break;
		}
		case 2:{ // see arXiv:1101.1381 for the galactic part and arXiv:1906.04750 for the extragalactic one (direct data from various instruments)
			sprintf(background_name,"./constraints/photons/background_gal_1101.1381.txt");
			background_file = fopen(background_name,"r+");
			if(!background_file){
				printf("\n\t\t\t[constraint_type2] : ERROR could not open file %s !",background_name);
				fflush(stdout);
				exit(0);
			}
			rewind(background_file);
			if(strncasecmp("#",dummy,1)){
				while((EOF!=fscanf(background_file,"%c",dummy)) && (strncasecmp("\n",dummy,1))){
					//printf("%s",dummy);
					//fflush(stdout);
				}
			}
			fscanf(background_file,"%s",dummy);
			fscanf(background_file,"%s",dummy);
			fscanf(background_file,"%i",&(type2->nb_bckg_gal));
			fscanf(background_file,"%s",dummy);
			fscanf(background_file,"%s",dummy);
			for(int i = 0;i<type2->nb_bckg_gal;i++){
				fscanf(background_file,"%lf",&(bckg_gal_ener[i]));
				fscanf(background_file,"%lf",&(background_gal[i]));
				background_gal[i] = background_gal[i]/pow(bckg_gal_ener[i],2.); // conversion from E².flux to flux in GeV^-1.cm^-2.s^-1.sr^-1
			}
			fclose(background_file);
			
			sprintf(background_name,"./constraints/photons/background_egal_1906.04750.txt");
			background_file = fopen(background_name,"r+");
			if(!background_file){
				printf("\n\t\t\t[constraint_type2] : ERROR could not open file %s !",background_name);
				fflush(stdout);
				exit(0);
			}
			rewind(background_file);
			if(strncasecmp("#",dummy,1)){
				while((EOF!=fscanf(background_file,"%c",dummy)) && (strncasecmp("\n",dummy,1))){
					//printf("%s",dummy);
					//fflush(stdout);
				}
			}
			fscanf(background_file,"%s",dummy);
			fscanf(background_file,"%s",dummy);
			fscanf(background_file,"%i",&(type2->nb_bckg_egal));
			fscanf(background_file,"%s",dummy);
			fscanf(background_file,"%s",dummy);
			for(int i = 0;i<type2->nb_bckg_egal;i++){
				fscanf(background_file,"%lf",&(bckg_egal_ener[i]));
				fscanf(background_file,"%lf",&(background_egal[i]));
			}
			fclose(background_file);
			break;
		}
		default:{
			printf("\n\t\t [background_type2] : ERROR wrong background type %i !",isatis->background_type);
			fflush(stdout);
			break;
		}
	}
	return;
}

void refined_flux(double *refined_flux,double *refined_ener,int nb_refined_en,double *flux,double *energies,int nb_ener){
	// this routine computes an energy refined version of the fluxes in order to correctly account for the energy resolution
	
	int counter = 0;
	for(int i = 0;i<nb_refined_en;i++){
		while(counter < nb_ener && energies[counter] < refined_ener[i]){
			counter++;
		}
		if(counter > 0 && counter < nb_ener && flux[counter-1] != 0. && flux[counter-1] != 0.){
			refined_flux[i] = pow(10.,log10(flux[counter-1]) + (log10(flux[counter]) - log10(flux[counter-1]))/(log10(energies[counter]) - log10(energies[counter-1]))*(log10(refined_ener[i]) - log10(energies[counter-1])));
		}
		else{
			refined_flux[i] = 0.;
		}
	}
	
	return;
}

double energy_resolution_type2(double E,struct type_2 *type2){
	// this routine computes the energy resolution relevant for a type2 constraint
	// be careful that E is given in GeV and epsilon is the relative dimensionless Delta(E)/E
	
	double epsilon;
	if(!strcmp(type2->name,"XGIS-THESEUS_2110.03333")){ // XGIS-THESEUS
		if(E < 0.15e-3){ // in GeV (see table 1 of arXiv:2110.03333)
			epsilon = 8.5e-2;
		}
		else{
			epsilon = 2.5e-2;
		}
		epsilon = 0.424*epsilon;
		return epsilon;
	}
	if(!strcmp(type2->name,"AdEPT_2010.04797")){ // AdEPT
		epsilon = 15e-2; // 15% energy resolution at 70 MeV see Table 1 of arXiv:1703.02546 or arXiv:1311.2059 (warning it is FWHM so it has to be divided by 2)
		return epsilon;
	}
	if(!strcmp(type2->name,"GRAMS_2010.04797")){ // GRAMS
		double delta_s = 0.01*E/sqrt(E/2.5e-3);
		double delta_e = 5e-6; // in GeV
		epsilon = sqrt(pow(delta_s/E,2.) + pow(delta_e/E,2.)); // see Section 3.3 of arXiv:1901.03430
		return epsilon;
	}
	if(!strcmp(type2->name,"eASTROGRAM_2010.04797") || !strcmp(type2->name,"AS-ASTROGRAM_2104.06168")){ // eASTROGRAM see arXiv:1611.02232 fig. 19, arXiv:1711.01265 fig. 1.3.1 or AS-ASTROGRAM see https://doi.org/10.22323/1.358.0579 fig. 3
		if(E < 4.96424e-4){
			epsilon = pow(10.,log10(6.93748e-6) + (log10(8.79440e-6) - log10(6.93748e-6))/(log10(4.96424e-4) - log10(3.01483e-4))*(log10(E) - log10(3.01483e-4)))/E;
		}
		else if(E < 9.89445e-4){
			epsilon = pow(10.,log10(8.79440e-6) + (log10(1.26766e-5) - log10(8.79440e-6))/(log10(9.89445e-4) - log10(4.96424e-4))*(log10(E) - log10(4.96424e-4)))/E;
		}
		else if(E < 2.01441e-3){
			epsilon = pow(10.,log10(1.26766e-5) + (log10(1.51446e-5) - log10(1.26766e-5))/(log10(2.01441e-3) - log10(9.89445e-4))*(log10(E) - log10(9.89445e-4)))/E;
		}
		else if(E < 4.96424e-3){
			epsilon = pow(10.,log10(1.51446e-5) + (log10(2.93637e-5) - log10(1.51446e-5))/(log10(4.96424e-3) - log10(2.01441e-3))*(log10(E) - log10(2.01441e-3)))/E;
		}
		else if(E < 1.00000e-2){
			epsilon = pow(10.,log10(2.93637e-5) + (log10(4.86066e-5) - log10(2.93637e-5))/(log10(1.00000e-2) - log10(4.96424e-3))*(log10(E) - log10(4.96424e-3)))/E;
		}
		else{
			epsilon = 0.30;
		}
		return epsilon;
	}
	if(!strcmp(type2->name,"GECCO_2010.04797")){ // GECCO see arXiv:2101.10370 Section 2
		epsilon = 1e-2;
		return epsilon;
	}
	if(!strcmp(type2->name,"PANGU_2010.04797")){ // PANGU see https://doi.org/10.22323/1.236.0964
		epsilon = 40e-2;
		return epsilon;
	}
	if(!strcmp(type2->name,"AMEGO_2010.04797")){ // AMEGO see Fig. 17 of https://asd.gsfc.nasa.gov/amego/files/AMEGO_Decadal_RFI.pdf see also arXiv:1907.07558 and arXiv:2101.03105 (also effective area here)
		if(E < 1e-2){ // in GeV
			epsilon = 1e-2/2.; // 1% FWHM
		}
		else{
			epsilon = 30e-2/2.; // 30% FWHM
		}
		return epsilon;
	}
	if(!strcmp(type2->name,"MAST_1902.01491")){ // MAST see fig. 2 of arXiv:1902.01491
		epsilon = 12.5e-2; // 12.5% mean on the 10 MeV - 1 GeV energy range
		return epsilon;
	}
	
	printf("\n\t\t [energy_resolution_type2] : ERROR wrong type2 constraint %s !",type2->name);
	fflush(stdout);
	exit(0);
}

double number_of_particles_type2(double *flux,double *energies,int nb_ener,double *Aeff,struct type_2 *type2,struct param *parameters){
	// this routine computes the number of particles observed by a telescope (see eq. (10) of arXiv:2110.03333)
	
	double result = 0.;
	double metric;
	double epsilon; // the energy resolution Delta(E)/E
	double R; // the spectral resolution function given in eq. (11) of arXiv:2110.03333
	for(int i = 0;i<type2->nb_ener-1;i++){
		for(int j = 0;j<nb_ener-1;j++){
			metric = (Aeff[(i+1)*2+0] - Aeff[i*2+0])*(energies[j+1] - energies[j]);
			epsilon = energy_resolution_type2(energies[j],type2);
			R = gaussian(Aeff[i*2+0],energies[j],epsilon); // the window function assures that the energy ranges are the same for all fluxes
			result += type2->T_obs*metric*Aeff[i*2+1]*R*flux[j];
		}
	}
	int nb_angles = 100;
	double *l = (double *)malloc(nb_angles*sizeof(double));
	double *b = (double *)malloc(nb_angles*sizeof(double));
	for(int i = 0;i<nb_angles;i++){
		l[i] = type2->l_min + (type2->l_max - type2->l_min)/(nb_angles - 1)*i;
		b[i] = type2->b_min + (type2->b_max - type2->b_min)/(nb_angles - 1)*i;
	}
	double Delta = 0.; // the field of view of the telescope
	for(int i = 0;i<nb_angles-1;i++){
		for(int j = 0;j<nb_angles-1;j++){
			Delta += fabs(cos(b[i]))*(l[i+1] - l[i])*(b[j+1] - b[j]);
		}
	}
	
	return result*Delta;
}

void binned_number_of_particles_type2(double *binned_flux,double *flux,double *energies,int nb_ener,double *Aeff,struct type_2 *type2,struct param *parameters){
	// this routine computes the (energy) binned number of particles taking the energies of Aeff as reference
	
	int nb_angles = 100;
	double *l = (double *)malloc(nb_angles*sizeof(double));
	double *b = (double *)malloc(nb_angles*sizeof(double));
	for(int i = 0;i<nb_angles;i++){
		l[i] = type2->l_min + (type2->l_max - type2->l_min)/(nb_angles - 1)*i;
		b[i] = type2->b_min + (type2->b_max - type2->b_min)/(nb_angles - 1)*i;
	}
	double Delta = 0.; // the field of view of the telescope
	for(int i = 0;i<nb_angles-1;i++){
		for(int j = 0;j<nb_angles-1;j++){
			Delta += fabs(cos(b[i]))*(l[i+1] - l[i])*(b[j+1] - b[j]);
		}
	}
	
	double result = 0.;
	double metric;
	double epsilon; // the energy resolution Delta(E)/E
	double R; // the spectral resolution function given in eq. (11) of arXiv:2110.03333
	for(int i = 0;i<type2->nb_ener-1;i++){
		binned_flux[i] = 0.;
		for(int j = 0;j<nb_ener-1;j++){
			metric = (Aeff[(i+1)*2+0] - Aeff[i*2+0])*(energies[j+1] - energies[j]);
			epsilon = energy_resolution_type2(energies[j],type2);
			R = gaussian(Aeff[i*2+0],energies[j],epsilon); // the window function assures that the energy ranges are the same for all fluxes
			binned_flux[i] += type2->T_obs*metric*Aeff[i*2+1]*R*flux[j]*Delta;
		}
	}
	
	return;
}

void binned_flux_type1(double *binned_flux,double *flux_refined,double *refined_ener,int nb_refined_en,double *flux,struct type_1 *type1){
	// this routine gives the flux integrated over the energy bins of the instruments
	
	int counter;
	for(int i = 0;i<type1->nb_ener;i++){
		counter = 0;
		binned_flux[i] = 0.;
		while(counter < nb_refined_en && refined_ener[counter] < flux[i*type1->format+0] - flux[i*type1->format+1]){ // finding the lower bound of the energy bin E - DeltaE-
			counter++;
		}
		if(counter > 0 && counter+1 < nb_refined_en){
			while(counter < nb_refined_en && refined_ener[counter] < flux[i*type1->format+0] + flux[i*type1->format+2]){ // spanning the energy bin until the upper bound E + DeltaE+
				binned_flux[i] += (refined_ener[counter+1] - refined_ener[counter])*(flux_refined[counter+1] + flux_refined[counter])/2.;
				counter++;
			}
		}
	}
	
	return;
}

void fill_type1(struct type_1 *type1){
	// this routine fills the type_1 structure depending on the instrument
	
	if(!strcmp(type1->name,"COMPTEL_1502.06116")){ // COMPTEL EGRB of arXiv:1502.06116
		type1->type = 1; // extragalactic
		return;
	}
	if(!strcmp(type1->name,"Fermi-LAT_1410.3696")){ // Fermi-LAT EGRB of arXiv:1502.06116 (by default column A for the constraint)
		type1->type = 1; // extragalactic
		return;
	}
	if(!strcmp(type1->name,"EGRET_0405441")){ // EGRET EGRB of arXiv:0405441
		type1->type = 1; // extragalactic
		return;
	}
	if(!strcmp(type1->name,"HEAO_balloon_9903492")){ // HEAO+balloon EGRB of arXiv:9903492
		type1->type = 1; // extragalactic
		return;
	}
	if(!strcmp(type1->name,"COMPTEL_1107.0200")){ // COMPTEL GC of arXiv:1107.0200
		type1->type = 0; // galactic
		type1->l_min = -30./180.*pi; // -30°
		type1->l_max = +30./180.*pi; // +30°
		type1->b_min = -15./180.*pi; // -15°
		type1->b_max = +15./180.*pi; // +15°
		return;
	}
	if(!strcmp(type1->name,"Fermi-LAT_1101.1381")){ // Fermi-LAT GC of arXiv:1101.1381
		type1->type = 0; // galactic
		type1->l_min = -30./180.*pi; // -30°
		type1->l_max = +30./180.*pi; // +30°
		type1->b_min = -10./180.*pi; // -10°
		type1->b_max = +10./180.*pi; // +10°
		return;
	}
	if(!strcmp(type1->name,"EGRET_9811211")){ // EGRET GC of arXiv:9811211
		type1->type = 0; // galactic
		type1->l_min = -30./180.*pi; // -30°
		type1->l_max = +30./180.*pi; // +30°
		type1->b_min = -5./180.*pi; // -5°
		type1->b_max = +5./180.*pi; // +5°
		return;
	}
	if(!strcmp(type1->name,"INTEGRAL_1107.0200")){ // INTEGRAL GC of arXiv:1107.0200
		type1->type = 0; // galactic
		type1->l_min = -30./180.*pi; // -30°
		type1->l_max = +30./180.*pi; // +30°
		type1->b_min = -15./180.*pi; // -15°
		type1->b_max = +15./180.*pi; // +15°
		return;
	}
	
	return;
}

void fill_type2(struct type_2 *type2){
	// this routine fills the type_2 structure depending on the instrument
	
	if(!strcmp(type2->name,"XGIS-THESEUS_2110.03333")){ // XGIS-THESEUS of arXiv:2110.03333
		type2->l_min = -5.*2.*pi/360.; // -5°
		type2->l_max = +5.*2.*pi/360.; // +5°
		type2->b_min = -5.*2.*pi/360.; // -5°
		type2->b_max = +5.*2.*pi/360.; // +5°
		type2->T_obs = 1e+8; // in s
		return;
	}
	if(!strcmp(type2->name,"AdEPT_2010.04797")){ // AdEPT of arXiv:2010.04797
		type2->l_min = -5.*2.*pi/360.; // -5°
		type2->l_max = +5.*2.*pi/360.; // +5°
		type2->b_min = -5.*2.*pi/360.; // -5°
		type2->b_max = +5.*2.*pi/360.; // +5°
		type2->T_obs = 1e+8; // in s
		return;
	}
	if(!strcmp(type2->name,"GRAMS_2010.04797")){ // GRAMS of arXiv:2010.04797
		type2->l_min = -5.*2.*pi/360.; // -5°
		type2->l_max = +5.*2.*pi/360.; // +5°
		type2->b_min = -5.*2.*pi/360.; // -5°
		type2->b_max = +5.*2.*pi/360.; // +5°
		type2->T_obs = 1e+8; // in s
		return;
	}
	if(!strcmp(type2->name,"eASTROGRAM_2010.04797")){ // eASTROGRAM of arXiv:2010.04797
		type2->l_min = -5.*2.*pi/360.; // -5°
		type2->l_max = +5.*2.*pi/360.; // +5°
		type2->b_min = -5.*2.*pi/360.; // -5°
		type2->b_max = +5.*2.*pi/360.; // +5°
		type2->T_obs = 1e+8; // in s
		return;
	}
	if(!strcmp(type2->name,"GECCO_2010.04797")){ // GECCO of arXiv:2101.10370, arXiv:2010.04797
		type2->l_min = -5.*2.*pi/360.; // -5°
		type2->l_max = +5.*2.*pi/360.; // +5°
		type2->b_min = -5.*2.*pi/360.; // -5°
		type2->b_max = +5.*2.*pi/360.; // +5°
		type2->T_obs = 1e+8; // in s
		return;
	}
	if(!strcmp(type2->name,"PANGU_2010.04797")){ // PANGU of arXiv:2010.04797
		type2->l_min = -5.*2.*pi/360.; // -5°
		type2->l_max = +5.*2.*pi/360.; // +5°
		type2->b_min = -5.*2.*pi/360.; // -5°
		type2->b_max = +5.*2.*pi/360.; // +5°
		type2->T_obs = 1e+8; // in s
		return;
	}
	if(!strcmp(type2->name,"AS-ASTROGRAM_2104.06168")){ // AS-ASTROGRAM of arXiv:2104.06168
		type2->l_min = -5.*2.*pi/360.; // -5°
		type2->l_max = +5.*2.*pi/360.; // +5°
		type2->b_min = -5.*2.*pi/360.; // -5°
		type2->b_max = +5.*2.*pi/360.; // +5°
		type2->T_obs = 1e+8; // in s
		return;
	}
	if(!strcmp(type2->name,"AMEGO_2010.04797")){ // AMEGO of arXiv:2010.04797
		type2->l_min = -5.*2.*pi/360.; // -5°
		type2->l_max = +5.*2.*pi/360.; // +5°
		type2->b_min = -5.*2.*pi/360.; // -5°
		type2->b_max = +5.*2.*pi/360.; // +5°
		type2->T_obs = 1e+8; // in s
		return;
	}
	if(!strcmp(type2->name,"MAST_1902.01491")){ // MAST of arXiv:2010.04797
		type2->l_min = -5.*2.*pi/360.; // -5°
		type2->l_max = +5.*2.*pi/360.; // +5°
		type2->b_min = -5.*2.*pi/360.; // -5°
		type2->b_max = +5.*2.*pi/360.; // +5°
		type2->T_obs = 1e+8; // in s
		return;
	}
	
	// if the constraint is not in the list
	printf("\n\t\t  [fill_type2] : ERROR wrong type2 constraint %s !",type2->name);
	fflush(stdout);
	exit(0);
}

double constraint_type2(double *fin_ener,double *inst_secondary_spectra,double *BH_masses,double *BH_params,double *BH_spectrum,double *dts,char *path,struct type_2 *type2,struct param *parameters,struct Isatis *isatis){
	// this routine calculates the constraint for type2 constraints: galactic and extragalactic constraint with background for prospective telescopes

	double result = -1.; // the final constraint f_PBH for this particular run
	
	printf("reading Aeff file...");
	fflush(stdout);
	
	char *Aeff_name = (char *)malloc(128*sizeof(char));
	sprintf(Aeff_name,"./constraints/photons/Aeff_%s.txt",type2->name);
	FILE *Aeff_file = fopen(Aeff_name,"r+");
	if(!Aeff_file){
		printf("\n\t\t\t[constraint_type2] : ERROR could not open file %s !",Aeff_name);
		fflush(stdout);
		exit(0);
	}
	rewind(Aeff_file);
	
	char dummy[64];
	if(strncasecmp("#",dummy,1)){
		while((EOF!=fscanf(Aeff_file,"%c",dummy)) && (strncasecmp("\n",dummy,1))){
			//printf("%s",dummy);
			//fflush(stdout);
		}
	}
	
	fscanf(Aeff_file,"%s",dummy);
	fscanf(Aeff_file,"%s",dummy);
	fscanf(Aeff_file,"%i",&(type2->nb_ener));
	
	// constraint file is of the type "energy / A_eff(E)"	
	double *Aeff = (double *)malloc(type2->nb_ener*2*sizeof(double *));
	
	for(int i = 0;i<2;i++){ // header of file
		fscanf(Aeff_file,"%s",dummy);
	}
	for(int i = 0;i<type2->nb_ener;i++){
		for(int j = 0;j<2;j++){
			fscanf(Aeff_file,"%lf",&(Aeff[i*2+j]));
		}
	}
	fclose(Aeff_file);
	free(Aeff_name);
	
	printf("applying constraint...");
	fflush(stdout);
	
	// normalization of the PBH distribution for galactic constraint
	double A = total_amplitude_DM(BH_masses,BH_params,BH_spectrum,parameters); // obtain the amplitude factor for extended mass/param distributions
	
	// compute the galactic contribution
	
	double *galactic = (double *)malloc(parameters->nb_fin_en*sizeof(double));
	
	galactic_contribution_type2(galactic,inst_secondary_spectra,A,type2,parameters,isatis);
	
	// compute the extra-galactic contribution
	int nb_stacked_en = 500;
	double Emin_stacked = 1e-6; // in GeV
	double Emax_stacked = 1e+19; // in GeV
	double *stacked_ener = (double *)malloc(nb_stacked_en*sizeof(double));
	double *extragalactic = (double *)malloc(nb_stacked_en*sizeof(double));
	for(int i = 0;i<nb_stacked_en;i++){
		stacked_ener[i] = pow(10.,log10(Emin_stacked) + (log10(Emax_stacked) - log10(Emin_stacked))/(nb_stacked_en - 1)*i);
		extragalactic[i] = 0.;
	}
	int today = 1; // we stop the redshifting at todays time
	int CMB = 1; // we start stacking after CMB
	
	// normalization of the PBH density function
	A = global_amplitude_DM(BH_masses,BH_params,BH_spectrum,parameters,isatis); // the amplitude factor for PBH=DM cosmologicaly
	
	extragalactic_contribution(extragalactic,stacked_ener,CMB,today,nb_stacked_en,fin_ener,dts,path,A,parameters,isatis);
	
	// compute the background
	double *bckg_gal_ener = (double *)malloc(tab_len*sizeof(double));
	double *background_gal = (double *)malloc(tab_len*sizeof(double));
	double *bckg_egal_ener = (double *)malloc(tab_len*sizeof(double));
	double *background_egal = (double *)malloc(tab_len*sizeof(double));
	
	background_type2(bckg_gal_ener,background_gal,bckg_egal_ener,background_egal,type2,isatis);
	
	// compute refined spectra in the correct energy range
	int nb_refined_en = 10000;
	double *refined_ener = (double *)malloc(nb_refined_en*sizeof(double)); // refined energy range in the range [E_min/10,E_max*10] around the Aeff table energy range
	for(int i = 0;i<nb_refined_en;i++){
		refined_ener[i] = pow(10.,log10(Aeff[0*2+0]/10.) + (log10(Aeff[(type2->nb_ener-1)*2+0]*10.) - log10(Aeff[0*2+0]/10.))/(nb_refined_en - 1)*i);
	}
	double *galactic_refined = (double *)malloc(nb_refined_en*sizeof(double));
	double *extragalactic_refined = (double *)malloc(nb_refined_en*sizeof(double));
	double *background_gal_refined = (double *)malloc(nb_refined_en*sizeof(double));
	double *background_egal_refined = (double *)malloc(nb_refined_en*sizeof(double));
	
	refined_flux(galactic_refined,refined_ener,nb_refined_en,galactic,fin_ener,parameters->nb_fin_en);
	refined_flux(extragalactic_refined,refined_ener,nb_refined_en,extragalactic,stacked_ener,nb_stacked_en);
	refined_flux(background_gal_refined,refined_ener,nb_refined_en,background_gal,bckg_gal_ener,type2->nb_bckg_gal);
	refined_flux(background_egal_refined,refined_ener,nb_refined_en,background_egal,bckg_egal_ener,type2->nb_bckg_egal);
	
	switch(isatis->energy_method){
		case 0:{ // in this case we compare the integrated number of photons over the full energy range with some SNR
			// compute the number of photons for signal and background
			double N_signal_gal = number_of_particles_type2(galactic_refined,refined_ener,nb_refined_en,Aeff,type2,parameters);
			double N_signal_egal = number_of_particles_type2(extragalactic_refined,refined_ener,nb_refined_en,Aeff,type2,parameters);
			double N_background_gal = number_of_particles_type2(background_gal_refined,refined_ener,nb_refined_en,Aeff,type2,parameters);
			double N_background_egal = number_of_particles_type2(background_egal_refined,refined_ener,nb_refined_en,Aeff,type2,parameters);
			
			// deduce the constraint
			result = 1./((N_signal_gal + N_signal_egal)/sqrt(N_background_gal + N_background_egal)/isatis->S_to_N);
			break;
		}
		case 1:{ // in this case we compare the number of photons by energy bins with some SNR
			// compute the binned number of photons for signal and background
			double *binned_signal_gal = (double *)malloc((type2->nb_ener-1)*sizeof(double));
			double *binned_signal_egal = (double *)malloc((type2->nb_ener-1)*sizeof(double));
			double *binned_background_gal = (double *)malloc((type2->nb_ener-1)*sizeof(double));
			double *binned_background_egal = (double *)malloc((type2->nb_ener-1)*sizeof(double));
			
			binned_number_of_particles_type2(binned_signal_gal,galactic_refined,refined_ener,nb_refined_en,Aeff,type2,parameters);
			binned_number_of_particles_type2(binned_signal_egal,extragalactic_refined,refined_ener,nb_refined_en,Aeff,type2,parameters);
			binned_number_of_particles_type2(binned_background_gal,background_gal_refined,refined_ener,nb_refined_en,Aeff,type2,parameters);
			binned_number_of_particles_type2(binned_background_egal,background_egal_refined,refined_ener,nb_refined_en,Aeff,type2,parameters);
			
			for(int i = 0;i<type2->nb_ener-1;i++){
				if((Aeff[i*2+0] > bckg_gal_ener[0] && Aeff[i*2+0] > bckg_egal_ener[0] && Aeff[i*2+0] < bckg_gal_ener[type2->nb_bckg_gal-1] && Aeff[i*2+0] < bckg_egal_ener[type2->nb_bckg_egal-1]) && (binned_signal_gal[i] + binned_signal_egal[i] != 0.) && (result == -1. || (result > 0. && 1./((binned_signal_gal[i] + binned_signal_egal[i])/sqrt(binned_background_gal[i] + binned_background_egal[i])/isatis->S_to_N) < result))){
					result = 1./((binned_signal_gal[i] + binned_signal_egal[i])/sqrt(binned_background_gal[i] + binned_background_egal[i])/isatis->S_to_N);
				}
			}
			
			free1D_double(binned_signal_gal);
			free1D_double(binned_signal_egal);
			free1D_double(binned_background_gal);
			free1D_double(binned_background_egal);
			break;
		}
		default:{
			printf("\n\t\t [constraints_type2] : ERROR wrong energy method %i !",isatis->energy_method);
			fflush(stdout);
			exit(0);
			break;
		}
	}
	
	free1D_double(Aeff);
	free1D_double(bckg_gal_ener);
	free1D_double(bckg_egal_ener);
	free1D_double(background_gal);
	free1D_double(background_egal);
	free1D_double(galactic);
	free1D_double(extragalactic);
	free1D_double(stacked_ener);
	free1D_double(galactic_refined);
	
	printf(" DONE",result);
	fflush(stdout);
	
	return result;
}

double constraint_type1(double *fin_ener,double *inst_secondary_spectra,double *BH_masses,double *BH_params,double *BH_spectrum,double *dts,char *path,struct type_1 *type1,struct param *parameters,struct Isatis *isatis){
	// this routine computes the constraints for existing instruments
	
	double result = -1.; // the final constraint f_PBH for this particular run
	
	printf("reading flux file...");
	fflush(stdout);
	
	char *flux_name = (char *)malloc(128*sizeof(char));
	sprintf(flux_name,"./constraints/photons/flux_%s.txt",type1->name);
	FILE *flux_file = fopen(flux_name,"r+");
	if(!flux_file){
		printf("\n\t\t\t[constraint_type1] : ERROR could not open file %s !",flux_name);
		fflush(stdout);
		exit(0);
	}
	rewind(flux_file);
	
	char dummy[64];
	if(strncasecmp("#",dummy,1)){
		while((EOF!=fscanf(flux_file,"%c",dummy)) && (strncasecmp("\n",dummy,1))){
			//printf("%s",dummy);
			//fflush(stdout);
		}
	}
	
	fscanf(flux_file,"%s",dummy);
	fscanf(flux_file,"%s",dummy);
	fscanf(flux_file,"%i",&(type1->nb_ener));
	fscanf(flux_file,"%s",dummy);
	fscanf(flux_file,"%s",dummy);
	fscanf(flux_file,"%i",&(type1->format));
	
	// flux file is of the type "energy(GeV) | DeltaE- | DeltaE+ | flux(GeV^-1.cm^-2.s^-1.sr^-1) | Deltaflux- | Deltaflux+"
	double *flux = (double *)malloc(type1->nb_ener*type1->format*sizeof(double *));
	
	for(int i = 0;i<type1->format;i++){ // header of file
		fscanf(flux_file,"%s",dummy);
	}
	for(int i = 0;i<type1->nb_ener;i++){
		for(int j = 0;j<type1->format;j++){
			fscanf(flux_file,"%lf",&(flux[i*type1->format+j]));
		}
	}
	fclose(flux_file);
	free(flux_name);
	
	printf("applying constraint...");
	fflush(stdout);
	
	double A; // normalization constant for the spectra
	
	switch(type1->type){
		case 0:{ // galactic flux
			// normalization of the PBH distribution for galactic constraint
			A = total_amplitude_DM(BH_masses,BH_params,BH_spectrum,parameters); // obtain the amplitude factor for extended mass/param distributions
			
			// compute the galactic contribution
			
			double *galactic = (double *)malloc(parameters->nb_fin_en*sizeof(double));
			
			galactic_contribution_type1(galactic,inst_secondary_spectra,A,type1,parameters,isatis);
			
			int nb_refined_en = 500;
			double *refined_ener = (double *)malloc(nb_refined_en*sizeof(double));
			for(int i = 0;i<nb_refined_en;i++){ // refined spectrum in the energy range [Emin/10,E_max*10] with the instrument energy range as reference
				refined_ener[i] = pow(10.,log10(flux[0*type1->format+0]/10.) + (log10(flux[(type1->nb_ener - 1)*type1->format+0]*10.) - log10(flux[0*type1->format+0]/10.))/(nb_refined_en - 1)*i);
			}
			double *galactic_refined = (double *)malloc(nb_refined_en*sizeof(double));
			
			refined_flux(galactic_refined,refined_ener,nb_refined_en,galactic,fin_ener,parameters->nb_fin_en);
			
			// computing the binned flux
			double *binned_galactic = (double *)malloc(type1->nb_ener*sizeof(double));
			
			binned_flux_type1(binned_galactic,galactic_refined,refined_ener,nb_refined_en,flux,type1);
			
			for(int i = 0;i<type1->nb_ener-1;i++){ // confidence_level*flux_plus denotes the x*superior_error_bar procedure
				if(isatis->confidence_level > 0){ // upper error bar
					if(binned_galactic[i] != 0. && (result == -1. || (result > 0. && 1./(binned_galactic[i]/((flux[i*type1->format+3] + isatis->confidence_level*flux[i*type1->format+5])*(flux[i*type1->format+2] + flux[i*type1->format+1]))) < result))){
						result = 1./(binned_galactic[i]/((flux[i*type1->format+3] + isatis->confidence_level*flux[i*type1->format+5])*(flux[i*type1->format+2] + flux[i*type1->format+1])));
					}
				}
				else{ // lower error bar
					if(binned_galactic[i] != 0. && (result == -1. || (result > 0. && 1./(binned_galactic[i]/((flux[i*type1->format+3] + isatis->confidence_level*flux[i*type1->format+4])*(flux[i*type1->format+2] + flux[i*type1->format+1]))) < result))){
						result = 1./(binned_galactic[i]/((flux[i*type1->format+3] + isatis->confidence_level*flux[i*type1->format+4])*(flux[i*type1->format+2] + flux[i*type1->format+1])));
					}
				}
			}
			
			free1D_double(galactic);
			free1D_double(galactic_refined);
			free1D_double(refined_ener);
			free1D_double(binned_galactic);
			break;
		}
		case 1:{ // extragalactic flux
			// compute the extra-galactic contribution
			int nb_stacked_en = 500;
			double Emin_stacked = 1e-6; // in GeV
			double Emax_stacked = 1e+19; // in GeV
			double *stacked_ener = (double *)malloc(nb_stacked_en*sizeof(double));
			double *extragalactic = (double *)malloc(nb_stacked_en*sizeof(double));
			for(int i = 0;i<nb_stacked_en;i++){
				stacked_ener[i] = pow(10.,log10(Emin_stacked) + (log10(Emax_stacked) - log10(Emin_stacked))/(nb_stacked_en - 1)*i);
				extragalactic[i] = 0.;
			}
			int today = 1; // we stop the redshifting at todays time
			int CMB = 1; // we start stacking after CMB
			
			// normalization of the PBH density function
			A = global_amplitude_DM(BH_masses,BH_params,BH_spectrum,parameters,isatis); // the amplitude factor for PBH=DM cosmologicaly
			
			extragalactic_contribution(extragalactic,stacked_ener,CMB,today,nb_stacked_en,fin_ener,dts,path,A,parameters,isatis);
			
			int nb_refined_en = 500;
			double *refined_ener = (double *)malloc(nb_refined_en*sizeof(double));
			for(int i = 0;i<nb_refined_en;i++){ // refined spectrum in the energy range [Emin/10,E_max*10] with the instrument energy range as reference
				refined_ener[i] = pow(10.,log10(flux[0*type1->format+0]/10.) + (log10(flux[(type1->nb_ener - 1)*type1->format+0]*10.) - log10(flux[0*type1->format+0]/10.))/(nb_refined_en - 1)*i);
			}
			double *extragalactic_refined = (double *)malloc(nb_refined_en*sizeof(double));
			
			refined_flux(extragalactic_refined,refined_ener,nb_refined_en,extragalactic,stacked_ener,nb_stacked_en);
			
			// computing the binned flux
			double *binned_extragalactic = (double *)malloc(type1->nb_ener*sizeof(double));
			
			binned_flux_type1(binned_extragalactic,extragalactic_refined,refined_ener,nb_refined_en,flux,type1);
			
			for(int i = 0;i<type1->nb_ener-1;i++){ // confidence_level*flux_plus denotes the x*superior_error_bar procedure
				if(isatis->confidence_level > 0){ // upper error bar
					if(binned_extragalactic[i] != 0. && (result == -1. || (result > 0. && 1./(binned_extragalactic[i]/((flux[i*type1->format+3] + isatis->confidence_level*flux[i*type1->format+5])*(flux[i*type1->format+2] + flux[i*type1->format+1]))) < result))){
						result = 1./(binned_extragalactic[i]/((flux[i*type1->format+3] + isatis->confidence_level*flux[i*type1->format+5])*(flux[i*type1->format+2] + flux[i*type1->format+1])));
					}
				}
				else{ // lower error bar
					if(binned_extragalactic[i] != 0. && (result == -1. || (result > 0. && 1./(binned_extragalactic[i]/((flux[i*type1->format+3] + isatis->confidence_level*flux[i*type1->format+4])*(flux[i*type1->format+2] + flux[i*type1->format+1]))) < result))){
						result = 1./(binned_extragalactic[i]/((flux[i*type1->format+3] + isatis->confidence_level*flux[i*type1->format+4])*(flux[i*type1->format+2] + flux[i*type1->format+1])));
					}
				}
			}
						
			free1D_double(extragalactic);
			free1D_double(stacked_ener);
			free1D_double(refined_ener);
			free1D_double(extragalactic_refined);
			free1D_double(binned_extragalactic);
			break;
		}
		default:{
			printf("\n\t\t [constraint_type1] : ERROR wrong constraint type %i !",type1->type);
			fflush(stdout);
			exit(0);
			break;
		}
	}
	
	free1D_double(flux);
	
	printf(" DONE");
	fflush(stdout);
	
	return result;
}

int main(int argc, char **argv){
	
	// START OF EXECUTION
	printf("[main] : start of execution\n");
	fflush(stdout);
	
	char runs[256];
	char params[256];
	
	if(argc < nb_params+1){
		printf("\n\t [main] : ERROR this program needs at least %i parameters!\n",nb_params);
		printf("[main] : end of execution");
		fflush(stdout);
		exit(0);
	}
	else{
		sscanf(argv[1],"%s",params);
		sscanf(argv[2],"%s",runs);
	}
	
	printf("[main] : reading Isatis params... ");
	fflush(stdout);
	
	struct Isatis isatis;
	
	read_Isatis(&isatis,params);
	
	printf("DONE\n");
	fflush(stdout);
	
	printf("[main] : reading runs paths... ");
	fflush(stdout);
	
	FILE *runs_file = fopen(runs,"r+");
	if(!runs_file){
		printf("\n\t [main] : ERROR could not open file %s !",runs);
		fflush(stdout);
		exit(0);
	}
	rewind(runs_file);
	
	char dummy[64];
	int nb_runs; // contains the number of runs to examine
	fscanf(runs_file,"%s",dummy);
	fscanf(runs_file,"%s",dummy);
	fscanf(runs_file,"%i",&(nb_runs));
	char **paths = (char **)malloc(nb_runs*sizeof(char *)); // contains the paths to the results of the runs
	double **constraints = (double **)malloc(nb_runs*sizeof(double *)); // contains the result constraints
	for(int i = 0;i<nb_runs;i++){
		paths[i] = (char *)malloc(256*sizeof(char));
		fscanf(runs_file,"%s",paths[i]);
		constraints[i] = (double *)malloc(tab_len*sizeof(double));
		for(int j = 0;j<tab_len;j++){
			constraints[i][j] = 0.;
		}
	}
	fclose(runs_file);
	
	printf("DONE\n");
	fflush(stdout);
	
	char *results_name = (char *)malloc(128*sizeof(char));
	sprintf(results_name,"results_photons_%s.txt",isatis.output);
	FILE *results_file = fopen(results_name,"w+");
	if(!(results_file)){
		printf("\n\t [main] : ERROR could not open file %s !");
		fflush(stdout);
		exit(0);
	}
	rewind(results_file);
	
	char **constraints_names = (char **)malloc(tab_len*sizeof(char *));
	for(int i = 0;i<tab_len;i++){
		constraints_names[i] = (char *)malloc(64*sizeof(char));
	}
	
	int counter = 0;
	
	int counter_existing = 0;
	// the existing X/gamma-ray constraints
	sprintf(constraints_names[counter],"COMPTEL_1502.06116");
	counter_existing++;
	counter++;
	sprintf(constraints_names[counter],"COMPTEL_1107.0200");
	counter_existing++;
	counter++;
	sprintf(constraints_names[counter],"EGRET_0405441");
	counter_existing++;
	counter++;
	sprintf(constraints_names[counter],"EGRET_9811211");
	counter_existing++;
	counter++;
	sprintf(constraints_names[counter],"Fermi-LAT_1410.3696");
	counter_existing++;
	counter++;
	sprintf(constraints_names[counter],"Fermi-LAT_1101.1381");
	counter_existing++;
	counter++;
	sprintf(constraints_names[counter],"INTEGRAL_1107.0200");
	counter_existing++;
	counter++;
	sprintf(constraints_names[counter],"HEAO_balloon_9903492");
	counter_existing++;
	counter++;
	
	int counter_prospective = 0;
	// the prospective X/gamma-ray constraints
	sprintf(constraints_names[counter],"AdEPT_2010.04797");
	counter_prospective++;
	counter++;
	sprintf(constraints_names[counter],"AMEGO_2010.04797");
	counter_prospective++;
	counter++;
	sprintf(constraints_names[counter],"AS-ASTROGRAM_2104.06168");
	counter_prospective++;
	counter++;
	sprintf(constraints_names[counter],"eASTROGRAM_2010.04797");
	counter_prospective++;
	counter++;
	sprintf(constraints_names[counter],"GECCO_2010.04797");
	counter_prospective++;
	counter++;
	sprintf(constraints_names[counter],"GRAMS_2010.04797");
	counter_prospective++;
	counter++;
	sprintf(constraints_names[counter],"MAST_1902.01491");
	counter_prospective++;
	counter++;
	sprintf(constraints_names[counter],"PANGU_2010.04797");
	counter_prospective++;
	counter++;
	sprintf(constraints_names[counter],"XGIS-THESEUS_2110.03333");
	counter_prospective++;
	counter++;
	
	fprintf(results_file,"%30s","run_path");
	for(int i = 0;i<counter;i++){
		fprintf(results_file,"%30s",constraints_names[i]);
	}
	fprintf(results_file,"\n");
	fflush(results_file);
	
	double *dts;
	double *BH_masses;
	double *BH_params;
	double *BH_spectrum;
	double *init_ener;
	double *fin_ener;
	double *inst_primary_spectra;
	double *inst_secondary_spectra;
	char *life_name;
	struct type_2 type2;
	struct type_1 type1;
	
	for(int i = 0;i<nb_runs;i++){
		
		// READING PARAMETERS
		printf("[main] : considering run... %i: %s\n",i+1,paths[i]); 
		printf("\t [main] : reading parameters... ");
		fflush(stdout);
		
		struct param parameters;
		char *parameters_name = (char *)malloc(256*sizeof(char));
		sprintf(parameters_name,"%s/results/%s/%s.txt",isatis.path,paths[i],paths[i]);
		
		read_params(&parameters,parameters_name);
		
		printf("DONE\n");
		fflush(stdout);
		
		// READING DATA
		printf("\t [main] : reading data...");
		fflush(stdout);
		
		BH_masses = (double *)malloc(parameters.BH_number*sizeof(double ));
		BH_params = (double *)malloc(parameters.param_number*sizeof(double));
		BH_spectrum = (double *)malloc(parameters.BH_number*parameters.param_number*sizeof(double));
		
		if(isatis.sessions == 0 || isatis.sessions == 2){
			init_ener = (double *)malloc(parameters.E_number*sizeof(double));
			fin_ener = (double *)malloc(parameters.nb_fin_en*sizeof(double));
			inst_primary_spectra = (double *)malloc(parameters.E_number*parameters.nb_init_part*sizeof(double));
			inst_secondary_spectra = (double *)malloc(parameters.nb_fin_en*parameters.nb_fin_part*sizeof(double));
		}
		
		if(isatis.sessions == 1 || isatis.sessions == 2){
			char *life_name = (char *)malloc(256*sizeof(char));
			sprintf(life_name,"%s/results/%s/life_evolutions.txt",isatis.path,paths[i]);
			FILE *life_file = fopen(life_name,"r+");
			if(!life_file){
				printf("\n\t\t [main] : ERROR could not open file %s !",life_name);
				fflush(stdout);
				exit(0);
			}
			rewind(life_file);
			for(int j = 0;j<16;j++){
				fscanf(life_file,"%s",dummy);
			}
			fscanf(life_file,"%i",&(parameters.nb_fin_times));
			fclose(life_file);
			dts = (double *)malloc(parameters.nb_fin_times*2*sizeof(double));
		}
		
		read_data(paths[i],BH_masses,BH_params,BH_spectrum,init_ener,fin_ener,inst_primary_spectra,inst_secondary_spectra,dts,&parameters,&isatis);
		
		printf("DONE\n");
		fflush(stdout);
		
		// APPLYING CONSTRAINTS
		printf("\t [main] : applying constraints...");
		fflush(stdout);
		
		printf("\n\t\t  Existing constraints :");
		fflush(stdout);
		counter = 0;
		for(int j = 0;j<counter_existing;j++){ // all the prospective X/gamma-ray constraints (type2)
			printf("\n\t\t  [%i] : %s : ",j+1,constraints_names[counter]);
			fflush(stdout);
			
			sprintf(type1.name,constraints_names[counter]);
			fill_type1(&type1);
			
			//constraints[i][counter] = constraint_type1(fin_ener,inst_secondary_spectra,BH_masses,BH_params,BH_spectrum,dts,paths[i],&type1,&parameters,&isatis);
			counter++;
		}
		
		printf("\n\t\t  Prospective constraints :");
		fflush(stdout);
		for(int j = 0;j<counter_prospective;j++){ // all the prospective X/gamma-ray constraints (type2)
			printf("\n\t\t  [%i] : %s : ",j+1,constraints_names[counter]);
			fflush(stdout);
			
			sprintf(type2.name,constraints_names[counter]);
			fill_type2(&type2);
			
			//constraints[i][counter] = constraint_type2(fin_ener,inst_secondary_spectra,BH_masses,BH_params,BH_spectrum,dts,paths[i],&type2,&parameters,&isatis);
			counter++;
		}
		
		fprintf(results_file,"%30s",paths[i]);
		for(int j = 0;j<counter_existing+counter_prospective;j++){
			fprintf(results_file,"%30.5e",constraints[i][j]);
		}
		fprintf(results_file,"\n");
		fflush(results_file);
		
		free1D_double(BH_masses);
		free1D_double(BH_params);
		free1D_double(init_ener);
		free1D_double(fin_ener);
		free1D_double(BH_spectrum);
		if(isatis.sessions == 0 || isatis.sessions == 2){
			free1D_double(inst_primary_spectra);
			free1D_double(inst_secondary_spectra);
		}
		if(isatis.sessions == 1 || isatis.sessions == 2){
			free1D_double(dts);
		}
		
		printf("\n\t DONE\n");
		fflush(stdout);
		
	}
	
	fclose(results_file);
	
	free2D_double(constraints,nb_runs);
	
	// warning remove "\n" for cygwin
	printf("[main] : end of execution\n");
	fflush(stdout);
	
	return 1;
}
