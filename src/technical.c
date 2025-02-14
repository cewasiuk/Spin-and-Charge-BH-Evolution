// This is the source file where some useful methods
// are implemented.
// Last modification: 15 September 2023
// Authors: Jérémy Auffinger jeremy.auffinger.blackhawk@protonmail.com & Alexandre Arbey alexandre.arbey@ens-lyon.fr

#include "include.h"

// *******FUNCTION THAT COMPUTES EXP(X)-1*******
double exp_adapt(double x){
	// This function computes exp(x) - 1 either with an expansion to 5th order
	
	if(x > 1e-5){
		return exp(x) - 1.;
	}
	else{
		return x + pow(x,2.)/2. + pow(x,3.)/6. + pow(x,4.)/24. + pow(x,5.)/120.;
	}
}

// *******FUNCTION THAT COMPUTES TRAPEZE INTEGRATION*******
double trapeze(double x1,double x2,double y1,double y2){
	// This function computes the trapeze integration of a function that takes
	// values y_i at point x_i.
	
	return 0.5*(x2-x1)*(y2+y1);
}

// *******FUNCTIONS THAT FREE MEMORY FROM MULTIDIMENTIONAL ARRAYS*******
void free1D_double(double *array){
	// This function frees the memory allocated to a 1D double array.
	
	if(array!=NULL) free(array);
	return;
}

void free1D_int(int *array){
	// This function frees the memory allocated to a 1D double array.
	
	if(array!=NULL) free(array);
	return;
}

void free2D_int(int **array,int l_1stD){
	// This function frees the memory allocated to a 2D int array
	// of 1st dimension l_1stD.
	
	if(array==NULL) return;

	for(int i = 0;i<l_1stD;i++){
		if(array[i]!=NULL) free(array[i]);
	}
	if(array!=NULL) free(array);
	return;
}

void free2D_double(double **array,int l_1stD){
	// This function frees the memory allocated to a 2D double array
	// of first dimension l_1stD.
	
	if(array==NULL) return;

	for(int i = 0;i<l_1stD;i++){
		if(array[i]!=NULL) free(array[i]);
	}
	if(array!=NULL) free(array);
	return;
}

void free2D_char(char **array,int l_1stD){
	// This function frees the memory allocated to a 2D char array
	// of first dimension l_1stD.
	
	if(array==NULL) return;

	for(int i = 0;i<l_1stD;i++){
		if(array[i]!=NULL) free(array[i]);
	}
	if(array!=NULL) free(array);
	return;
}

void free3D_int(int ***array,int l_1stD,int l_2ndD){
	// This function frees the memory allocated to a 3D double array
	// of 1st dimension l_1stD and 2nd dimension l_2ndD.
	
	if(array==NULL) return;

	for(int i = 0;i<l_1stD;i++){
		for(int j = 0;j<l_2ndD;j++){
			if(array[i][j]!=NULL) free(array[i][j]);
		}
		if(array[i]!=NULL) free(array[i]);
	}
	if(array!=NULL) free(array);
	return;
}

void free3D_double(double ***array,int l_1stD,int l_2ndD){
	// This function frees the memory allocated to a 3D double array
	// of 1st dimension l_1stD and 2nd dimension l_2ndD.
	
	if(array==NULL) return;

	for(int i = 0;i<l_1stD;i++){
		for(int j = 0;j<l_2ndD;j++){
			if(array[i][j]!=NULL) free(array[i][j]);
		}
		if(array[i]!=NULL) free(array[i]);
	}
	if(array!=NULL) free(array);
	return;
}

void free4D_double(double ****array,int l_1stD,int l_2ndD,int l_3rdD){
	// This function frees the memory allocated to a 4D double array
	// of 1st dimension l_1stD, 2nd dimension l_2ndD and 3rd dimension
	// l_3rdD.
	
	if(array==NULL) return;
	
	for(int i = 0;i<l_1stD;i++){
		for(int j = 0;j<l_2ndD;j++){
			for(int k = 0;k<l_3rdD;k++){
				if(array[i][j][k]!=NULL) free(array[i][j][k]);
			}
			if(array[i][j]!=NULL) free(array[i][j]);
		}
		if(array[i]!=NULL) free(array[i]);
	}
	if(array!=NULL) free(array);
	return;
}

// *******FUNCTION THAT FINDS A MAXIMUM*******
int ind_max(double *table,int llength){
	// This function returns the index of the maximum of
	// an int array.
	
	int index = 0;
	double temp = table[0];
	for(int i = 0;i<llength;i++){
		if(table[i] > temp){
			index = i;
			temp = table[i];
		}
	}
	return index;
}

// *******FUNCTIONS THAT PERFORM A FUSION SORTING*******
void fusion(double *table,int start1,int end1,int end2){
    double *table1;
	table1 = (double *)malloc((end1-start1+1)*sizeof(double));
    int start2 = end1+1;
    int count1 = start1;
    int count2 = start2;
    int i;
    for(i = start1;i<=end1;i++){
        table1[i-start1] = table[i];
    }
    for(i = start1;i<=end2;i++){
        if(count1 == start2){
			break;
        }
        else if(count2 == (end2+1)){
            table[i] = table1[count1-start1];
            count1++;
        }
        else if(table1[count1-start1]<table[count2]){
            table[i] = table1[count1-start1];
            count1++;
        }
        else{
            table[i] = table[count2];
            count2++;
        }
    }
    free1D_double(table1);
	return;
}

void sort_fusion_bis(double *table,int start,int end){
    if(start!=end){
        int middle=(end+start)/2;
        sort_fusion_bis(table,start,middle);
        sort_fusion_bis(table,middle+1,end);
        fusion(table,start,middle,end);
    }
	return;
}

void sort_fusion(double *table,int llength){
    if(llength>0){
        sort_fusion_bis(table,0,llength-1);
    }
	return;
}
