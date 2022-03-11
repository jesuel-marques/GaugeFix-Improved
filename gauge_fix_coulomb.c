//	gcc -o gauge_fix_coulomb gauge_fix_coulomb.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/math_ops.c source/gauge_fixing.c source/fourvector_field.c  -lm -O4 -march=skylake

#include <stdio.h>					//	Standard header files in C
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "SU3_parameters.h"			//	Simulation parameters
#include "source/SU3_global.h"		//	Global variables


#include "source/lattice.h"			//	Initialization functions and calculations of
									//	positionitions and links on the lattice.
#include "source/SU3_ops.h"

#include "source/gauge_fixing.h"	//	Specific functions about the gauge-fixing
#include "source/fourvector_field.h"


int main(){

	char filename[max_length_name];

	double complex * U = (double complex *) malloc(Volume * d * 3 * 3 * sizeof(double complex));
	double complex * A = (double complex *) malloc(Volume * d * 3 * 3 * sizeof(double complex));

	double complex * G = (double complex *) malloc(Volume * 3 * 3 * sizeof(double complex));
	//	Gauge transformation to bring U to Landau-gauge

	if ( U == NULL || A == NULL || G == NULL) {
        printf("Memory allocation failed for configuration or gauge-transformation");
        exit(1); 
    }


	SU2_identity = (double *) malloc(4 * sizeof(double));
	SU2_null = (double *) malloc(4 * sizeof(double));

	SU3_identity = (double complex *) malloc(3 * 3 * sizeof(double complex));
	SU3_null = (double complex *) malloc(3 * 3 * sizeof(double complex));

	if ( SU2_identity == NULL || SU2_null == NULL 
		|| SU3_identity == NULL || SU3_null == NULL) {
        printf("Memory allocation failed for SU2 or SU3 identity or null matrices");
        exit(1); 
    }

	SU3_initialize();
	
	for ( int config = 1; config <= 1; config++ ) {
		SU3_set_gauge_transf_to_identity(G);

	    //  load each configuration based on name template

		sprintf(filename, "configs/Config_1_beta_6.000_Nxyz_8_Nt_8.txt");
		SU3_load_U(filename, U);

		SU3_update_A(U, G, A);
        //  fix the gauge
        SU3_gauge_fix(U, G, A);

		// print the gauge fixed configuration based on template name
		// print_gaugefixedconfig();
	
	}

	free(U); 
	free(A);
	free(G);

	free(SU2_identity);	free(SU2_null);
	free(SU3_identity);	free(SU3_null);
	
	return 0;
}