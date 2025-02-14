// convert_tables is a side-program aimed at converting the computed hadronization tables into
// C arrays to make BlackHawk runs faster by including them at the compilation.
// Last modification: 07 April 2024
// Authors: Jérémy Auffinger jeremy.auffinger.blackhawk@protonmail.com & Alexandre Arbey alexandre.arbey@ens-lyon.fr

#include "src/include.h"

// &&&&&&&&&&&&&&&&&&&&& MAIN &&&&&&&&&&&&&&&&&&&&&&&&
// This is the convert_tables main routine. It reads
// the tables contained in pythia_tables, herwig_tables
// pythia_tables_new, hazma_tables and hdmspectra_tables
// and rewrites them into a header file 
// src/tables/hadronization_tables.h.
int main(int argc,char **argv)
{
	printf("\n\n");
	printf("\t ############################\n");
	printf("\t #      BLACKHAWK v2.3      #\n");
	printf("\t #     HAWKING SPECTRUM     #\n");
	printf("\t #    COMPUTATION DEVICE    #\n");
	printf("\t ############################\n\n");
	fflush(stdout);
	
	printf("[main] : STARTING EXECUTION...\n");
	
	struct param parameters;
	parameters.full_output = 1;
	
	// *******READING and CONVERTING THE HADRONIZATION TABLES IN FOLDERS "*_tables/"*******
	printf("[main] : CONVERTING HADRONIZATION TABLES...");
	fflush(stdout);
	
	FILE *table_file;
	table_file=fopen("./src/tables/hadronization_tables/hadronization_tables_pythia.h","w");
	fclose(table_file);
	table_file=fopen("./src/tables/hadronization_tables/hadronization_tables_herwig.h","w");
	fclose(table_file);
	table_file=fopen("./src/tables/hadronization_tables/hadronization_tables_pythianew.h","w");
	fclose(table_file);
	table_file=fopen("./src/tables/hadronization_tables/hadronization_tables_hazma.h","w");
	fclose(table_file);
	table_file=fopen("./src/tables/hadronization_tables/hadronization_tables_hdmspectra.h","w");
	fclose(table_file);

	for(int ie=0;ie<=4;ie++)
	{
		parameters.hadronization_choice=ie;
		
		switch(parameters.hadronization_choice){
			case 0:{ // PYTHIA tables (BBN epoch)
				parameters.nb_init_en = 250;
				parameters.nb_fin_en = 500;
				parameters.nb_init_part = 14;
				parameters.nb_fin_part = 11;
				break;
			}
			case 1:{ // HERWIG tables (BBN epoch)
				parameters.nb_init_en = 100;
				parameters.nb_fin_en = 100;
				parameters.nb_init_part = 14;
				parameters.nb_fin_part = 11;
				break;
			}
			case 2:{ // PYTHIA tables (present epoch)
				parameters.nb_init_en = 250;
				parameters.nb_fin_en = 500;
				parameters.nb_init_part = 14;
				parameters.nb_fin_part = 6;
				break;
			}
			case 3:{ // HAZMA tables (present epoch)
				parameters.nb_init_en = 250;
				parameters.nb_fin_en = 500;
				parameters.nb_init_part = 3;
				parameters.nb_fin_part = 2;
				break;
			}
			case 4:{ // HDMSpectra tables
				parameters.nb_init_en = 500;
				parameters.nb_fin_en = 1000;
				parameters.nb_init_part = 17;
				parameters.nb_fin_part = 6;
				break;
			}
			default:{
				printf("\n\t [main] : ERROR WRONG HADRONIZATION CHOICE !\n");
				fflush(stdout);
				return 0;
				break;
			}
		}
		
		
		double ****tables; // contains the branching ratios dN(E',E) of the hadronization tables
		tables = (double ****)malloc(parameters.nb_fin_part*sizeof(double ***));
		for(int i = 0;i<parameters.nb_fin_part;i++){
			tables[i] = (double ***)malloc(parameters.nb_init_en*sizeof(double **));
			for(int j = 0;j<parameters.nb_init_en;j++){
				tables[i][j] = (double **)malloc(parameters.nb_fin_en*sizeof(double *));
				for(int k = 0;k<parameters.nb_fin_en;k++){
					tables[i][j][k] = (double *)malloc(parameters.nb_init_part*sizeof(double));
				}
			}
		}
		double *initial_energies; // contains the tabulated initial energies
		double *final_energies; // contains the tabulated final energies
		initial_energies = (double *)malloc(parameters.nb_init_en*sizeof(double));
		final_energies = (double *)malloc(parameters.nb_fin_en*sizeof(double));
		
		read_hadronization_tables(tables,initial_energies,final_energies,&parameters);
		
		printf("\t writing...");
		fflush(stdout);
		
		convert_hadronization_tables(tables,initial_energies,final_energies,&parameters);
		
		free4D_double(tables,parameters.nb_fin_part,parameters.nb_init_en,parameters.nb_fin_en); // freeing arrays not used anymore
		free(initial_energies); // freeing arrays not used anymore
		free(final_energies); // freeing last arrays
		
		printf("\tDONE");
		fflush(stdout);
	}

	printf("\n\n\t DONE\n");
	fflush(stdout);
	
	printf("[main] : END OF EXECUTION\n");
	fflush(stdout);
	
	return 1; // this is the only possibility of returning 1: it means that all steps have succeeded
}
// &&&&&&&&&&&&&&&&&&&&&&&&& END OF MAIN &&&&&&&&&&&&&&&&&&&&&&&&&&&&&
