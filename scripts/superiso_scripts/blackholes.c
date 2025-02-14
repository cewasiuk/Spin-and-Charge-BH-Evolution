#include "src/include.h"
#include "src/include_dm.h"




/*******************************************************************/
/******************************************************************/
int main(int argc,char** argv)
{
	char name[500];

  	if(argc<3) 
  	{ 
    		printf(" This program needs 2 parameter:\n"
           	"name of the secondary spectrum file and mass of the blackhole in gram\n");
      		exit(1); 
  	} 
	else 
  	{
  		

	double mass=6.62e23*atof(argv[2]); //mass in GeV
	double fraction=1.;// fraction of PBH mass density over DM density, assuming that PBHs only accounts for a fraction od DM and follows DM density distribution in dwarf spheroidal galaxies
	
	printf("\n--------Fermi-lat gamma-rays-----------\n");
	int test=test_fermi(argv[1], mass, fraction);

	if(test==0)printf("Valid point\n");
	if(test==1)printf("Excluded point\n");
	
	
	printf("\n--------AMS-02 antiprotons-----------\n");
	int test_ams= indirect_ams02_calculator_decay(argv[1], mass, fraction);
	
	if(test_ams==0)printf("Valid point\n");
	if(test_ams==1)printf("Excluded point\n");
	
	}
	
	
	
	
	
	return 1;
}
