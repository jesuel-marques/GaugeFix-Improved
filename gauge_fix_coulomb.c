 //	gcc -o gauge_fix_coulomb gauge_fix_coulomb.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/math_ops.c source/gauge_fixing.c source/fourvector_field.c  -lm -O4 -march=skylake
//	icc -o gauge_fix_coulomb gauge_fix_coulomb.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/math_ops.c source/gauge_fixing.c source/fourvector_field.c  -lm -fast-g

#include <stdio.h>					//	Standard header files in C
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include <omp.h>

#include "SU3_parameters.h"			//	Simulation parameters


#include "source/lattice.h"			//	Initialization functions and calculations of
									//	positions and links on the lattice.

#include "source/gauge_fixing.h"	//	Specific functions involved in the gauge-fixing


int main(){
	char config_filename[max_length_name];

	#pragma omp parallel for num_threads(8) private(config_filename)
			for (int config = 1; config <= 100; config++ ) {
				double complex * U = (double complex *) malloc(Volume * d * 3 * 3 * sizeof(double complex));
				double complex * G = (double complex *) malloc(Volume * 3 * 3 * sizeof(double complex));
				//	Gauge transformation to bring U to Landau-gauge

				if ( U == NULL || G == NULL) {
					//	Test if allocation was successful.
					printf("Memory allocation failed for configuration or gauge-transformation");
					exit(1); 

				}
				

				//  load each configuration based on name template

				sprintf(config_filename, "configs/NewFormConfig_%d_beta_5.700_Nxyz_8_Nt_8.txt",config);
				SU3_load_config(config_filename, U);
				
				//	Sets gauge transformation matrices to identity
				SU3_set_gauge_transf_to_identity(G);

				//  fix the gauge
				SU3_gauge_fix(U, G, config);

				// print the gauge fixed configuration based on template name
				// print_gaugefixedconfig();

				free(U);	//	Free allocated memory for the configuration
				free(G);	//	and gauge transformation.
			
			}
	return 0;
}