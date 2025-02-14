// This is the source file where general methods of
// BlackHawk are implemented.
// Last modification: 07 April 2024
// Authors: Jérémy Auffinger jeremy.auffinger.blackhawk@protonmail.com & Alexandre Arbey alexandre.arbey@ens-lyon.fr

#include "include.h"

// *******GENERAL FUNCTIONS*******
int read_params(struct param *parameters, char name[],int session){
	// This function reads the file 'parameters.txt' and fills
	// the parameters structure. It updates the value of success
	// to 1 (no error in the parameters) or 0 (errors). It also uses
	// the other reading functions to get the numerical tables info.
	
	FILE *param_file;
	param_file = fopen(name,"r");
	if(!param_file){
		printf("\n\t [read_params] : ERROR COULD NOT OPEN FILE '%s' !\n",name);
		fflush(stdout);
		exit(0);
	}
	rewind(param_file);
	
	char dummy[512];	// contains irrelevant characters of the parameters file (2 for each parameter, 1st is its name and 2nd is =)
	
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
			
			if(parameters->full_output != 0 && parameters->full_output != 1)
			{
				printf("\n\t [read_params] : ERROR WRONG OUTPUT CHOICE !\n");
				fflush(stdout);
				fclose(param_file);
				return 0;
			}
		}
		else if(!strcasecmp(dummy,"interpolation_method"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->interpolation_method));
			
			if(parameters->interpolation_method != 0 && parameters->interpolation_method != 1)
			{
				printf("\n\t [read_params] : ERROR WRONG INTERPOLATION METHOD !\n");
				fflush(stdout);
				fclose(param_file);
				return 0;
			}
		}
		else if(!strcasecmp(dummy,"metric"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->metric));
			if(parameters->metric < 0 || parameters->metric > 3)
			{
				printf("\n\t [read_params] : ERROR WRONG BLACK HOLE METRIC !\n");
				fflush(stdout);
				fclose(param_file);
				return 0;
			}
		}
		else if(!strcasecmp(dummy,"BH_number"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->BH_number));
			if(parameters->BH_number < 1)
			{
				printf("\n\t [read_params] : ERROR WRONG BLACK HOLE INITIAL MASSES NUMBER !\n");
				fflush(stdout);
				fclose(param_file);
				return 0;
			}
		}
		else if(!strcasecmp(dummy,"Mmin"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(parameters->Mmin));
			parameters->Mmin = parameters->Mmin*mass_conversion; // conversion from CGS to GeV
			
			if(parameters->Mmin < 0.)
			{
				printf("\n\t [read_params] : ERROR WRONG BLACK HOLE MASS BOUNDARIES !\n");
				fflush(stdout);
				fclose(param_file);
				return 0;
			}			
		}
		else if(!strcasecmp(dummy,"Mmax"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%lf",&(parameters->Mmax));
			parameters->Mmax = parameters->Mmax*mass_conversion; // conversion from CGS to GeV
		}
		else if(!strcasecmp(dummy,"param_number"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->param_number));
			if(parameters->param_number < 1)
			{
				printf("\n\t [read_params] : ERROR WRONG BLACK HOLE SECONDARY PARAMETER NUMBER !\n");
				fflush(stdout);
				fclose(param_file);
				return 0;
			}
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
			
			if(parameters->spectrum_choice < -1 || (parameters->spectrum_choice > 5 && parameters->spectrum_choice != 11))
			{
				printf("\n\t [read_params] : ERROR WRONG MASS SPECTRUM CHOICE %i!\n",parameters->spectrum_choice);
				fflush(stdout);
				fclose(param_file);
				return 0;
			}
		}
		else if(!strcasecmp(dummy,"spectrum_choice_param"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->spectrum_choice_param));
			
			if(parameters->spectrum_choice_param < 0 || parameters->spectrum_choice_param > 2)
			{
				printf("\n\t [read_params] : ERROR WRONG SPIN SPECTRUM CHOICE !\n");
				fflush(stdout);
				fclose(param_file);
				return 0;
			}
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
			
			if(parameters->tmin <= 0.)
			{
				printf("\n\t [read_params] : ERROR WRONG INITIAL TIME FOR BLACK HOLES EVOLUTION !\n");
				fflush(stdout);
				fclose(param_file);
				return 0;
			}	
			parameters->tmin = parameters->tmin*time_conversion; // conversion from CGS to GeV
			
			
			if(!(parameters->tmin_manual)){
				double gamma_BH = 0.2; // fraction of particle horizon mass that collapses into a BH during radiation-dominated phase (see Arbey et al. 2012.09867)
				parameters->tmin = 1./gamma_BH*parameters->Mmin/pow(Mp,2.);
			}
		}
		else if(!strcasecmp(dummy,"limit"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->limit));
			
			if(session == 0 && parameters->limit < 1)
			{
				printf("\n\t [read_params] : ERROR WRONG RECURSION LIMIT !\n");
				fflush(stdout);
				fclose(param_file);
				return 0;
			}
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
			
			if(parameters->E_number < 1)
			{
				printf("\n\t [read_params] : ERROR WRONG NUMBER OF PRIMARY EMISSION ENERGIES !\n");
				fflush(stdout);
				fclose(param_file);
				return 0;
			}
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
			
			if(parameters->grav != 0 && parameters->grav != 1)
			{
				printf("\n\t [read_params] : ERROR WRONG GRAVITON CHOICE !\n");
				fflush(stdout);
				fclose(param_file);
				return 0;
			}
		}
		else if(!strcasecmp(dummy,"add_DM"))
		{
			fscanf(param_file,"%s",dummy);
			fscanf(param_file,"%i",&(parameters->add_DM));
			
			if(parameters->add_DM != 0 && parameters->add_DM != 1)
			{
				printf("\n\t [read_params] : ERROR WRONG DM CHOICE !\n");
				fflush(stdout);
				fclose(param_file);
				return 0;
			}
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
			
			if(parameters->hadronization_choice < 0 || parameters->hadronization_choice > 4)
			{
				printf("\n\t [read_params] : ERROR WRONG HADRONIZATION CHOICE !\n");
				fflush(stdout);
				fclose(param_file);
				return 0;
			}
		}
	}
	fclose(param_file);
	
	if(session == 0){
		if(!(read_fM_infos(parameters) == 1)){
			return 0; // there was an error in the fM_tables infos file
		}
	}
	
	if(!(read_gamma_infos(parameters) == 1)){
		return 0; // there was an error in the gamma_tables infos file
	}
	
	if(parameters->primary_only == 0){
		if(!(read_hadronization_infos(parameters)) == 1){
			return 0; // there was an error in the hadronization tables infos file
		}
	}
	
	if(parameters->metric == 0 && (parameters->amin < parameters->param_min || parameters->amax > parameters->param_max))
	{
		printf("\n\t [read_params] : ERROR WRONG BH SPINS !\n");
		fflush(stdout);
		return 0;
	}
	
	if(parameters->metric == 1){
		if(parameters->epsilon_LQG < parameters->param_min || parameters->epsilon_LQG > parameters->param_max){
			printf("\n\t [read_params] : ERROR WRONG EPSILON LQG PARAMETER !\n");
			fflush(stdout);
			return 0;
		}
		if(parameters->a0_LQG != 0. && parameters->a0_LQG != 0.11){
			printf("\n\t [read_params] : ERROR WRONG a0 LQG PARAMETER !\n");
			fflush(stdout);
			return 0;
		}
		parameters->param_number = 1; // epsilon is the same for all BHs
		parameters->spectrum_choice_param = 0;
		if(parameters->add_DM && parameters->spin_DM == 1.5){
			printf("\n\t [read_params] : ERROR POLYMERIZED HR WAS NOT COMPUTED FOR SPIN 3/2 !\n");
			fflush(stdout);
			return 0;
		}
	}
	
	if(parameters->metric == 2){
		if(session == 0){
			printf("\n\t [read_params] : ERROR BH EVOLUTION IS NOT AVAILABLE WITH THIS METRIC !\n");
			fflush(stdout);
			return 0;
		}
		if(parameters->add_DM == 1 && parameters->spin_DM == 1.5){
			printf("\n\t [read_params] : ERROR CHARGED BH HR WAS NOT COMPUTED FOR SPIN 3/2 !\n");
			fflush(stdout);
			return 0;
		}
		if(parameters->Qmin < parameters->param_min || parameters->Qmax > parameters->param_max || (parameters->param_number > 1 && parameters->Qmin >= parameters->Qmax)){
			printf("\n\t [read_params] : ERROR WRONG BH CHARGE !\n");
			fflush(stdout);
			return 0;
		}
		printf("\n\t [read_params] : WARNING GREYBODY FACTORS ARE CORRECT FOR NEUTRAL PARTICLES ONLY !\n");
		fflush(stdout);
	}
	
	if(parameters->metric == 3){
		if(session == 0){
			printf("\n\t [read_params] : ERROR BH EVOLUTION IS NOT AVAILABLE WITH THIS METRIC !\n");
			fflush(stdout);
			return 0;
		}
		if(parameters->add_DM == 1 && parameters->spin_DM == 1.5){
			printf("\n\t [read_params] : ERROR HIGHER DIMENSIONAL BH HR WAS NOT COMPUTED FOR SPIN 3/2 !\n");
			fflush(stdout);
			return 0;
		}
		parameters->param_number = 1; // n is the same for all BHs
		parameters->spectrum_choice_param = 0;
		parameters->M_star = 1.;
		if(parameters->n < parameters->param_min || parameters->n > parameters->param_max){
			printf("\n\t [read_params] : ERROR WRONG NUMBER OF EXTRA DIMENSIONS!\n");
			fflush(stdout);
			return 0;
		}
	}
	
	if(parameters->spectrum_choice==1)
	{
		parameters->amplitude=amplitude_lognormal;
		parameters->stand_dev=stand_dev_lognormal;
		parameters->crit_mass=crit_mass_lognormal;
		
		parameters->amplitude = parameters->amplitude*mass_conversion/pow(leng_conversion,3.); // conversion from CGS to GeV
	}
	if(parameters->spectrum_choice==11)
	{
		parameters->amplitude=amplitude_lognormal2;
		parameters->stand_dev=stand_dev_lognormal;
		parameters->crit_mass=crit_mass_lognormal;
		
		parameters->amplitude = parameters->amplitude/pow(leng_conversion,3.); // conversion from CGS to GeV
	}
	else if(parameters->spectrum_choice==2)
	{
		parameters->amplitude=amplitude_powerlaw;
		parameters->eqstate=eqstate_powerlaw;
		
		double gamma = 2.*parameters->eqstate/(1.+parameters->eqstate);
		parameters->amplitude = parameters->amplitude*pow(mass_conversion,gamma-1.)/pow(leng_conversion,3.); // conversion from CGS to GeV
	
	}
	else if(parameters->spectrum_choice==3)
	{
		parameters->amplitude=amplitude_critical_collapse;
		parameters->crit_mass=crit_mass_critical_collapse;
		parameters->amplitude = parameters->amplitude*pow(mass_conversion,-2.85)/pow(leng_conversion,3.); // conversion from CGS to GeV
	}
	else if(parameters->spectrum_choice == 5)
	{
		parameters->amplitude = amplitude_uniform/pow(leng_conversion,3.); // conversion from CGS to GeV
	}
		
	if(parameters->spectrum_choice==1 && parameters->stand_dev < 0.)
	{
		printf("\n\t [read_params] : ERROR WRONG MASS SPECTRUM STANDARD DEVIATION !\n");
		fflush(stdout);
		return 0;
	}
	
	if((parameters->spectrum_choice==1 || parameters->spectrum_choice==3) && parameters->crit_mass < 0.)
	{
		printf("\n\t [read_params] : ERROR WRONG MASS SPECTRUM CRITICAL MASS !\n");
		fflush(stdout);
		return 0;
	}
	parameters->crit_mass = parameters->crit_mass*mass_conversion; // conversion from CGS to GeV
	
	if(parameters->amplitude < 0. && (parameters->spectrum_choice == 1 || parameters->spectrum_choice == 2 || parameters->spectrum_choice == 3 || parameters->spectrum_choice == 5)){
		printf("\n\t [read_params] : ERROR WRONG MASS SPECTRUM AMPLITUDE !\n");
		fflush(stdout);
		return 0;
	}
	
	if(parameters->spectrum_choice_param == 2){
		parameters->stand_dev_param = stand_dev_param_gaussian;
		if(parameters->stand_dev_param < 0.){
			printf("\n\t [read_params] : ERROR WRONG SPIN SPECTRUM STANDARD DEVIATION !\n");
			fflush(stdout);
			return 0;
		}
		parameters->mean_param = mean_param_gaussian;
	}
	
	if(parameters->BH_number > 1 && parameters->Mmin > parameters->Mmax)
	{
		printf("\n\t [read_params] : ERROR WRONG MASS SPECTRUM BOUNDARIES !\n");
		fflush(stdout);
		return 0;
	}
	
	if(parameters->add_DM)
	{
		if(parameters->m_DM < 0.){
			printf("\n\t [read_params] : ERROR WRONG DM MASS !\n");
			fflush(stdout);
			return 0;
		}
		if(parameters->spin_DM != 0. && parameters->spin_DM != 1. && parameters->spin_DM != 2. && parameters->spin_DM != 0.5 && parameters->spin_DM != 1.5){
			printf("\n\t [read_params] : ERROR WRONG DM SPIN !\n");
			fflush(stdout);
			return 0;
		}
		if(parameters->dof_DM < 0.){
			printf("\n\t [read_params] : ERROR WRONG DM DOFS !\n");
			fflush(stdout);
			return 0;
		}
	}
	if(session == 0)
	{
		if(parameters->Mmin < parameters->Mmin_fM || parameters->Mmax > parameters->Mmax_fM)
		{
			printf("\n\t [read_params] : ERROR MASS BOUNDARIES NOT COMPATIBLE WITH THE EVOLUTION TABLES !");
			printf("\n\t [read_params] : TABLES WERE COMPUTED FOR %.2e g < M < %.2e g\n",parameters->Mmin_fM/mass_conversion,parameters->Mmax_fM/mass_conversion);
			fflush(stdout);
			return 0;
		}
		
		if(parameters->add_DM && parameters->dof_DM > 1.){
			printf("\n\t [read_params] : WARNING THE EVOLUTION TABLES HAVE BEEN TABULATED FOR A SINGLE DM DOF !\n");
			fflush(stdout);
		}
	}
	
	if(session == 0){
		if(parameters->BH_remnant != 0 && parameters->BH_remnant != 1){
			printf("\n\t [read_params] : ERROR WRONG BH REMNANT CHOICE !\n");
			fflush(stdout);
			exit(0);
		}
		parameters->M_remnant = parameters->M_remnant*mass_conversion;
		if(parameters->M_remnant < Mp){
			printf("\n\t [read_params] : ERROR WRONG BH REMNANT MASS !\n");
			fflush(stdout);
			exit(0);
		}
	}
	
	if(parameters->Emin < 0. || parameters-> Emax < 0. || parameters->Emax < parameters->Emin)
	{
		printf("\n\t [read_params] : ERROR WRONG ENERGY BOUNDARIES !\n");
		fflush(stdout);
		exit(0);
	}
	
	if(parameters->primary_only == 0)
	{			
		if(parameters->Emin < parameters->Emin_hadro || parameters->Emax > parameters->Emax_hadro)
		{
			if(parameters->hadronization_choice==0)
			{
				printf("\n\n\t [read_params] : WARNING ENERGY BOUNDARIES ARE NOT COMPATIBLE WITH THE PYTHIA HADRONIZATION TABLES !");
				printf("\n\t [read_params] : PYTHIA HADRONIZATION TABLES WERE COMPUTED FOR %.2e GEV < E < %.2e GEV",parameters->Emin_hadro,parameters->Emax_hadro);
				printf("\n\t [read_params] : EXTRAPOLATION WILL BE USED\n\n");
				fflush(stdout);
			}
			else if(parameters->hadronization_choice==1)
			{
				printf("\n\n\t [read_params] : WARNING ENERGY BOUNDARIES ARE NOT COMPATIBLE WITH THE HERWIG HADRONIZATION TABLES !");
				printf("\n\t [read_params] : HERWIG HADRONIZATION TABLES WERE COMPUTED FOR %.2e GeV < E < %.2e GeV",parameters->Emin_hadro,parameters->Emax_hadro);
				printf("\n\t [read_params] : EXTRAPOLATION WILL BE USED\n\n");
				fflush(stdout);
			}
			else if(parameters->hadronization_choice==2)
			{
				printf("\n\n\t [read_params] : WARNING ENERGY BOUNDARIES ARE NOT COMPATIBLE WITH THE NEW PYTHIA HADRONIZATION TABLES !");
				printf("\n\t [read_params] : NEW PYTHIA HADRONIZATION TABLES WERE COMPUTED FOR %.2e GeV < E < %.2e GeV",parameters->Emin_hadro,parameters->Emax_hadro);
				printf("\n\t [read_params] : EXTRAPOLATION WILL BE USED\n\n");
				fflush(stdout);
			}
			else if(parameters->hadronization_choice==3){
				printf("\n\n\t [read_params] : ERROR EXTRAPOLATION IS NOT AVAILABLE WITH THE HAZMA TABLES !\n");
				printf("\n\t [read_params] : HAZMA HADRONIZATION TABLES WERE COMPUTED FOR %.2e GeV < E < %.2e GeV\n\n",parameters->Emin_hadro,parameters->Emax_hadro);
				fflush(stdout);
				exit(0);
			}
			else if(parameters->hadronization_choice==4){
				printf("\n\n\t [read_params] : ERROR EXTRAPOLATION IS NOT AVAILABLE WITH THE HDMSPECTRA TABLES !\n");
				printf("\n\t [read_params] : HDMSPECTRA HADRONIZATION TABLES WERE COMPUTED FOR %.2e GeV < E < %.2e GeV\n\n",parameters->Emin_hadro,parameters->Emax_hadro);
				fflush(stdout);
				exit(0);
			}
		}
		
	}
				
	if(parameters->spectrum_choice == 0 && parameters->BH_number > 1)
	{
		printf("\n\t [read_params] : WARNING THE DIRAC SPECTRUM WORKS ONLY WITH ONE BH MASS !\n");
		fflush(stdout);
		parameters->BH_number = 1;
	}
	if(!(parameters->spectrum_choice == -1) && (parameters->metric == 0 || parameters->metric == 2) && parameters->spectrum_choice_param == 0 && parameters->param_number > 1)
	{
		printf("\n\t [read_params] : WARNING THE DIRAC SPECTRUM WORKS ONLY WITH ONE BH ADDITIONAL PARAMETER !\n");
		fflush(stdout);
		parameters->param_number = 1;
	}
	
	if(parameters->hadronization_choice == 3){
		if(parameters->add_DM == 0){
			parameters->add_DM = 1;
			parameters->m_DM = m_pi0; // we have m_pi0 < m_pipm and we put back the constraint of E>m_pipm manually for the FSR 
			parameters->spin_DM = 0.;
			parameters->dof_DM = 1.;
		}
		else{
			printf("\n [read_params] : ERROR THE DM PARTICLE IS ALREADY DEFINED AS THE PION TO COMPUTE THE HADRONIZATION WITH HAZMA !\n");
			fflush(stdout);
			exit(0);
		}
	}
	
	return 1; // this is the only way of returning 1 -> all the parameters are okay
}

int read_fM_infos(struct param *parameters){
	// This function reads the information about the f(M,x_j)/g(M,a^*) numerical tables
	// and fills the parameters structure with the corresponding values.
	
	char file_name[500];
	sprintf(file_name,"./src/tables/fM_tables/infos.txt");
	FILE *infos_file = fopen(file_name,"r");
	if(!infos_file){
		printf("\n\t [read_fM_infos] : ERROR COULD NOT OPEN FILE '%s' !\n",file_name);
		fflush(stdout);
		exit(0);
	}
	rewind(infos_file);
	
	char dummy[64];
	while(EOF != fscanf(infos_file,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(infos_file,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		if(!strcasecmp(dummy,"Mmin_fM"))
		{
			fscanf(infos_file,"%s",dummy);
			fscanf(infos_file,"%lf",&(parameters->Mmin_fM));
		}
		else if(!strcasecmp(dummy,"Mmax_fM"))
		{
			fscanf(infos_file,"%s",dummy);
			fscanf(infos_file,"%lf",&(parameters->Mmax_fM));
		}
		else if(!strcasecmp(dummy,"nb_fM_masses"))
		{
			fscanf(infos_file,"%s",dummy);
			fscanf(infos_file,"%i",&(parameters->nb_fM_masses));
		}
		else if(!strcasecmp(dummy,"nb_fM_a"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->metric == 0){
				fscanf(infos_file,"%i",&(parameters->nb_fM_param));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_fM_epsilon"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->metric == 1){
				fscanf(infos_file,"%i",&(parameters->nb_fM_param));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
	}
	fclose(infos_file);
	
	parameters->Mmin_fM = parameters->Mmin_fM*mass_conversion;
	parameters->Mmax_fM = parameters->Mmax_fM*mass_conversion;
	
	if(parameters->Mmin_fM < 0. || parameters->Mmax_fM < 0. || parameters->Mmin_fM > parameters->Mmax_fM){
		printf("\n\t [read_fM_infos] : ERROR WRONG MASS BOUNDARIES !\n");
		fflush(stdout);
		return 0;
	}
	if(parameters->nb_fM_masses < 1 || parameters->nb_fM_param < 1){
		printf("\n\t [read_fM_infos] : ERROR WRONG TABLE PARAM DIMENSIONS !\n");
		fflush(stdout);
		return 0;
	}
	
	return 1; // this is the only way of returning 1: there was no error in the parameters
}

int read_gamma_infos(struct param *parameters){
	// This function reads the informations about the greybody factors and fits tables
	// and fills the parameters structure with the relevant informations.
	
	char file_name[500];
	sprintf(file_name,"./src/tables/gamma_tables/infos.txt");
	FILE *infos_file = fopen(file_name,"r");
	if(!infos_file){
		printf("\n\t [read_gamma_infos] : ERROR COULD NOT OPEN FILE '%s' !\n",file_name);
		fflush(stdout);
		exit(0);
	}
	rewind(infos_file);
	
	char dummy[64];
	while(EOF != fscanf(infos_file,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(infos_file,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		if(!strcasecmp(dummy,"particle_number"))
		{
			fscanf(infos_file,"%s",dummy);
			fscanf(infos_file,"%i",&(parameters->particle_number));
		}
		else if(!strcasecmp(dummy,"nb_gamma_a"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->metric == 0){
				fscanf(infos_file,"%i",&(parameters->nb_gamma_param));
				parameters->nb_gamma_spins = 5;
				parameters->nb_gamma_fits = 7;
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_gamma_epsilon"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->metric == 1){
				fscanf(infos_file,"%i",&(parameters->nb_gamma_param));
				parameters->nb_gamma_spins = 4;
				parameters->nb_gamma_fits = 3;
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_gamma_Q"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->metric == 2){
				fscanf(infos_file,"%i",&(parameters->nb_gamma_param));
				parameters->nb_gamma_spins = 4;
				parameters->nb_gamma_fits = 3;
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_gamma_n"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->metric == 3){
				fscanf(infos_file,"%i",&(parameters->nb_gamma_param));
				parameters->nb_gamma_spins = 4;
				parameters->nb_gamma_fits = 3;
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_gamma_x"))
		{
			fscanf(infos_file,"%s",dummy);
			fscanf(infos_file,"%i",&(parameters->nb_gamma_x));
		}
		else if(!strcasecmp(dummy,"a_min"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->metric == 0){
				fscanf(infos_file,"%lf",&(parameters->param_min));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"a_max"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->metric == 0){
				fscanf(infos_file,"%lf",&(parameters->param_max));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"epsilon_min"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->metric == 1){
				fscanf(infos_file,"%lf",&(parameters->param_min));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"epsilon_max"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->metric == 1){
				fscanf(infos_file,"%lf",&(parameters->param_max));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"Q_min"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->metric == 2){
				fscanf(infos_file,"%lf",&(parameters->param_min));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"Q_max"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->metric == 2){
				fscanf(infos_file,"%lf",&(parameters->param_max));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"n_min"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->metric == 3){
				fscanf(infos_file,"%lf",&(parameters->param_min));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"n_max"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->metric == 3){
				fscanf(infos_file,"%lf",&(parameters->param_max));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
	}
	fclose(infos_file);
	
	if(parameters->particle_number < 1){
		printf("\n\t [read_gamma_infos] : ERROR WRONG PARTICLE NUMBER !\n");
		fflush(stdout);
		return 0;
	}
	if(parameters->nb_gamma_param < 1 || parameters->nb_gamma_x < 1){
		printf("\n\t [read_gamma_infos] : ERROR WRONG TABLE X DIMENSIONS !\n");
		fflush(stdout);
		return 0;
	}
	
	if(parameters->metric == 0 && (parameters->param_min < 0. || parameters->param_min > 1. || parameters->param_max < 0. || parameters->param_max > 1. || parameters->param_min > parameters->param_max)){
		printf("\n\t [read_gamma_infos] : ERROR WRONG SPIN BOUNDARIES !\n");
		fflush(stdout);
		return 0;
	}
	if(parameters->metric == 1 && (parameters->param_min < 0. || parameters->param_min > 100. || parameters->param_max < 0. || parameters->param_max > 100. || parameters->param_min > parameters->param_max)){
		printf("\n\t [read_gamma_infos] : ERROR WRONG EPSILON BOUNDARIES !\n");
		fflush(stdout);
		return 0;
	}
	if(parameters->metric == 2 && (parameters->param_min < 0. || parameters->param_min > 1. || parameters->param_max < 0. || parameters->param_max > 1. || parameters->param_min > parameters->param_max)){
		printf("\n\t [read_gamma_infos] : ERROR WRONG CHARGE BOUNDARIES !\n");
		fflush(stdout);
		return 0;
	}
	if(parameters->metric == 3 && (parameters->param_min < 0 || parameters->param_max < 0 || parameters->param_min > parameters->param_max)){
		printf("\n\t [read_gamma_infos] : ERROR WRONG DIMENSION BOUNDARIES !\n");
		fflush(stdout);
		return 0;
	}
	
	return 1; // this is the only way of returning 1: there was no error in the parameters
}

int read_hadronization_infos(struct param *parameters){
	// This function reads the informations about the hadronization tables
	// and fills the parameters structure with the relevant informations.
	
	char file_name[500];
	sprintf(file_name,"./src/tables/hadronization_tables/infos.txt");
	FILE *infos_file = fopen(file_name,"r");
	if(!infos_file){
		printf("\n\t [read_hadronization_infos] : ERROR COULD NOT OPEN FILE '%s' !\n",file_name);
		fflush(stdout);
		exit(0);
	}
	rewind(infos_file);
	
	char dummy[64];
	while(EOF != fscanf(infos_file,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(infos_file,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		// PYHTIA TABLES PARAMETERS
		if(!strcasecmp(dummy,"Emin_hadro_pythia"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 0){
				fscanf(infos_file,"%lf",&(parameters->Emin_hadro));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"Emax_hadro_pythia"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 0){
				fscanf(infos_file,"%lf",&(parameters->Emax_hadro));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_init_en_pythia"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 0){
				fscanf(infos_file,"%i",&(parameters->nb_init_en));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_fin_en_pythia"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 0){
				fscanf(infos_file,"%i",&(parameters->nb_fin_en));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_init_part_pythia"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 0){
				fscanf(infos_file,"%i",&(parameters->nb_init_part));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_fin_part_pythia"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 0){
				fscanf(infos_file,"%i",&(parameters->nb_fin_part));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		// HERWIG TABLES PARAMETERS
		if(!strcasecmp(dummy,"Emin_hadro_herwig"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 1){
				fscanf(infos_file,"%lf",&(parameters->Emin_hadro));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"Emax_hadro_herwig"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 1){
				fscanf(infos_file,"%lf",&(parameters->Emax_hadro));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_init_en_herwig"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 1){
				fscanf(infos_file,"%i",&(parameters->nb_init_en));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_fin_en_herwig"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 1){
				fscanf(infos_file,"%i",&(parameters->nb_fin_en));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_init_part_herwig"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 1){
				fscanf(infos_file,"%i",&(parameters->nb_init_part));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_fin_part_herwig"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 1){
				fscanf(infos_file,"%i",&(parameters->nb_fin_part));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		// PYTHIA NEW TABLES PARAMETERS
		if(!strcasecmp(dummy,"Emin_hadro_pythia_new"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 2){
				fscanf(infos_file,"%lf",&(parameters->Emin_hadro));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"Emax_hadro_pythia_new"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 2){
				fscanf(infos_file,"%lf",&(parameters->Emax_hadro));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_init_en_pythia_new"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 2){
				fscanf(infos_file,"%i",&(parameters->nb_init_en));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_fin_en_pythia_new"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 2){
				fscanf(infos_file,"%i",&(parameters->nb_fin_en));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_init_part_pythia_new"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 2){
				fscanf(infos_file,"%i",&(parameters->nb_init_part));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_fin_part_pythia_new"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 2){
				fscanf(infos_file,"%i",&(parameters->nb_fin_part));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		// HAZMA TABLES PARAMETERS
		if(!strcasecmp(dummy,"Emin_hadro_hazma"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 3){
				fscanf(infos_file,"%lf",&(parameters->Emin_hadro));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"Emax_hadro_hazma"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 3){
				fscanf(infos_file,"%lf",&(parameters->Emax_hadro));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_init_en_hazma"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 3){
				fscanf(infos_file,"%i",&(parameters->nb_init_en));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_fin_en_hazma"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 3){
				fscanf(infos_file,"%i",&(parameters->nb_fin_en));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_init_part_hazma"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 3){
				fscanf(infos_file,"%i",&(parameters->nb_init_part));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_fin_part_hazma"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 3){
				fscanf(infos_file,"%i",&(parameters->nb_fin_part));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		// HDMSPECTRA TABLES PARAMETERS
		if(!strcasecmp(dummy,"Emin_hadro_hdmspectra"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 4){
				fscanf(infos_file,"%lf",&(parameters->Emin_hadro));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"Emax_hadro_hdmspectra"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 4){
				fscanf(infos_file,"%lf",&(parameters->Emax_hadro));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_init_en_hdmspectra"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 4){
				fscanf(infos_file,"%i",&(parameters->nb_init_en));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_fin_en_hdmspectra"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 4){
				fscanf(infos_file,"%i",&(parameters->nb_fin_en));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_init_part_hdmspectra"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 4){
				fscanf(infos_file,"%i",&(parameters->nb_init_part));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
		else if(!strcasecmp(dummy,"nb_fin_part_hdmspectra"))
		{
			fscanf(infos_file,"%s",dummy);
			if(parameters->hadronization_choice == 4){
				fscanf(infos_file,"%i",&(parameters->nb_fin_part));
			}
			else{
				fscanf(infos_file,"%s",dummy);
			}
		}
	}
	fclose(infos_file);
	
	if(parameters->Emin_hadro < 0. || parameters->Emax_hadro < 0. || parameters->Emin_hadro > parameters->Emax_hadro){
		printf("\n\t [read_hadronization_infos] : ERROR WRONG ENERGY BOUNDARIES !\n");
		fflush(stdout);
		return 0;
	}
	if(parameters->nb_init_en < 1 || parameters->nb_fin_en < 1){
		printf("\n\t [read_hadronization_infos] : ERROR WRONG TABLE DIMENSIONS !\n");
		fflush(stdout);
		return 0;
	}
	if(parameters->nb_init_part < 1 || parameters->nb_fin_part < 1){
		printf("\n\t [read_hadronization_infos] : ERROR WRONG NUMBER OF PARTICLES !\n");
		fflush(stdout);
		return 0;
	}
	
	return 1; // this is the only way of returning 1: there was no error in the parameters
}

int memory_estimation(struct param *parameters,int session){
	// This function computes and displays an estimation of the
	// used RAM and static memory throughout the desired run,
	// depending on the input parameters. The user has the choice
	// to cancel the run if the full_output parameter is set to 1.
	// This will return 0 and thus stop the BlackHawk execution.
	
	double mean_evolution = 1000.; // mean evolution_length of one BH
	double RAM_use; // contains the RAM estimation, in Bytes
	double static_use; // contains the static memory estimation, in Bytes
	switch(session){
		case 0:{ // session 0: BlackHawk_tot
			double RAM_use_steps[3]; // 3 distinct steps
			// definition of the RAM occupied by the arrays
			double init_masses = 8.*(parameters->BH_number);
			double init_params = 8.*(parameters->param_number);
			double spec_table = 8.*(parameters->BH_number)*(parameters->param_number);
			double fM_table = 8.*(parameters->nb_fM_param)*(parameters->nb_fM_masses);
			double gM_table = 8.*(parameters->nb_fM_param)*(parameters->nb_fM_masses);
			double fM_masses = 8.*(parameters->nb_fM_masses);
			double fM_a = 8.*(parameters->nb_fM_param);
			double evol_times = 8.*(parameters->BH_number)*(parameters->param_number);
			double rank = 4.*(parameters->BH_number)*(parameters->param_number);
			double sorted_times = 8.*(parameters->BH_number)*(parameters->param_number);
			double life_masses = 8.*(parameters->BH_number)*(parameters->param_number)*(parameters->limit);
			double life_params = 8.*(parameters->BH_number)*(parameters->param_number)*(parameters->limit);
			double life_times = 8.*(parameters->limit);
			double dts = 8.*(parameters->limit);
			double gammas = 8.*4.*(parameters->nb_gamma_param)*(parameters->nb_gamma_x);
			double gamma_param = 8.*(parameters->nb_gamma_param);
			double gamma_x = 8.*(parameters->nb_gamma_x);
			double fits = 8.*4.*(parameters->nb_gamma_param)*parameters->nb_gamma_fits;
			double dof = 8.*(parameters->particle_number + parameters->grav + parameters->add_DM);
			double masses_primary = 8.*(parameters->particle_number + parameters->grav + parameters->add_DM);
			double masses_secondary = 8.*(parameters->nb_fin_part);
			double spins = 8.*(parameters->particle_number + parameters->grav + parameters->add_DM);
			double times = 8.*mean_evolution*(parameters->BH_number)*(parameters->param_number);
			double energies = 8.*(parameters->E_number);
			double tables = 8.*(parameters->nb_fin_part)*(parameters->nb_init_en)*(parameters->nb_fin_en)*(parameters->nb_fin_part);
			double initial_energies = 8.*(parameters->nb_init_en);
			double final_energies = 8.*(parameters->nb_fin_en);
			double partial_hadronized_spectra = 8.*(parameters->nb_fin_part)*(parameters->E_number)*(parameters->nb_fin_en);
			double partial_primary_spectra = 8.*(parameters->particle_number + parameters->grav + parameters->add_DM)*(parameters->E_number);
			double partial_integrated_hadronized_spectra = 8.*(parameters->nb_fin_part)*(parameters->nb_fin_en);
			RAM_use_steps[0] = init_masses + init_params + spec_table; // step 1
			RAM_use_steps[1] = init_masses + init_params + spec_table + fM_table + fM_a + fM_masses + gM_table + evol_times + rank + sorted_times + life_masses + life_params + life_times + dts; // step 2
			RAM_use_steps[2] = spec_table + life_masses + life_params + life_times + dts + gammas + gamma_param + gamma_x + fits + dof + spins + masses_primary + times + energies + partial_primary_spectra; // step 3
			if(parameters->primary_only == 0){
				RAM_use_steps[2] += tables + initial_energies + final_energies + partial_hadronized_spectra + partial_integrated_hadronized_spectra + masses_secondary; // step 4
			}
			RAM_use = RAM_use_steps[ind_max(RAM_use_steps,3)]; // finding the maximum use
			// definition of the disk memory used by the output files
			double BH_spectrum = 12.*(parameters->BH_number + 1.)*(parameters->param_number + 1.);
			double life_evolutions = 12.*(parameters->BH_number)*(1. + 2.*parameters->param_number)*(parameters->BH_number)*(parameters->param_number)*mean_evolution;
			double primary_spectrum = 12.*(parameters->particle_number + parameters->grav + parameters->add_DM)*(parameters->E_number)*mean_evolution*(parameters->BH_number)*(parameters->param_number);
			double secondary_spectrum;
			if(parameters->primary_only == 0){
				secondary_spectrum = 12.*(parameters->nb_fin_part)*(parameters->nb_fin_en)*mean_evolution*(parameters->BH_number)*(parameters->param_number);
			}
			else{ // skipping hadronization
				secondary_spectrum = 0.;
			}
			static_use = BH_spectrum + life_evolutions + primary_spectrum + secondary_spectrum;
			break;
		}
		case 1:{ // session 1: BlackHawk_inst
			double RAM_use_steps[4]; // 4 distinctive steps
			double BH_masses = 8.*(parameters->BH_number);
			double life_params = 8.*(parameters->param_number);
			double spec_table = 8.*(parameters->BH_number)*(parameters->param_number);
			double gammas = 8.*4.*(parameters->nb_gamma_param)*(parameters->nb_gamma_x);
			double gamma_param = 8.*(parameters->nb_gamma_param);
			double gamma_x = 8.*(parameters->nb_gamma_x);
			double fits = 8.*4.*(parameters->nb_gamma_param)*(parameters->nb_gamma_fits);
			double dof = 8.*(parameters->particle_number + parameters->grav + parameters->add_DM);
			double masses1 = 8.*(parameters->particle_number + parameters->grav + parameters->add_DM);
			double masses2 = 8.*(parameters->nb_fin_part);
			double spins = 8.*(parameters->particle_number + parameters->grav + parameters->add_DM);
			double instantaneous_primary_spectra = 8.*(parameters->particle_number + parameters->grav + parameters->add_DM)*(parameters->E_number);
			double energies = 8.*(parameters->E_number);
			double tables = 8.*(parameters->nb_fin_part)*(parameters->nb_init_en)*(parameters->nb_fin_en)*(parameters->nb_fin_part);
			double initial_energies = 8.*(parameters->nb_init_en);
			double final_energies = 8.*(parameters->nb_fin_en);
			double instantaneous_hadronized_spectra = 8.*(parameters->nb_fin_part)*(parameters->E_number)*(parameters->nb_fin_en);
			double instantaneous_integrated_hadronized_spectra = 8.*(parameters->nb_fin_part)*(parameters->nb_fin_en);
			RAM_use_steps[0] = BH_masses + life_params + spec_table; // step 1
			RAM_use_steps[1] = BH_masses + spec_table + gammas + gamma_param + gamma_x + fits + dof + spins + masses1 + instantaneous_primary_spectra + energies; // step 2
			if(parameters->primary_only == 0){
				RAM_use_steps[2] = instantaneous_primary_spectra + energies + tables + initial_energies + final_energies + masses2 + instantaneous_hadronized_spectra; // step 3
				RAM_use_steps[3] = instantaneous_primary_spectra + energies + initial_energies + final_energies + instantaneous_hadronized_spectra + instantaneous_integrated_hadronized_spectra; // step 4
			}
			else{ // skipping hadronization
				RAM_use_steps[2] = 0.;
				RAM_use_steps[3] = 0.;
			}
			RAM_use = RAM_use_steps[ind_max(RAM_use_steps,4)]; // finding the maximum use
			// definition of the disk memory used by the output files
			double BH_spectrum = 12.*(parameters->BH_number + 1.)*(parameters->param_number +1.);
			double primary_spectrum = 12.*(parameters->E_number)*(parameters->particle_number + parameters->grav + parameters->add_DM);
			double secondary_spectrum;
			if(parameters->primary_only == 0){
				secondary_spectrum = 12.*(parameters->nb_fin_en)*(parameters->nb_fin_part);
			}
			else{ // skipping hadronization
				secondary_spectrum = 0.;
			}
			static_use = BH_spectrum + primary_spectrum + secondary_spectrum;
			break;
		}
		default:{
			printf("\n\t [memory_estimation] : ERROR WRONG SESSION !\n");
			fflush(stdout);
			exit(0);
			break;
		}
	}
	if(parameters->full_output){
		printf("\n\n\t Running this session will use around %.3f MB of RAM and %.3f MB of disc memory.\n",RAM_use/1.e+6,static_use/1.e+6); // conversion to MB
		printf("\n\t Do you want to continue? (type y or n) ");
		fflush(stdout);
		char answer;
		scanf("%c",&answer);
		if(answer == 'y'){ // in this case the run will continue
			printf("\n");
			return 1;
		}
		else if(answer == 'n'){ // in this case the run will stop
			return 0;
		}
		else{
			printf("\n\t [memory_estimation] : ERROR WRONG ANSWER !\n");
			fflush(stdout);
			exit(0);
		}
	}
	else{
		return 1;
	}
}
