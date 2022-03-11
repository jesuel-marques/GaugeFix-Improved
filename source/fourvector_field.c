#include <complex.h>
#include <stdlib.h>

#include "../SU3_parameters.h"
#include "SU3_global.h"

#include "lattice.h"

#include "math_ops.h"
#include "SU3_ops.h"

					
void SU3_update_A(double complex * U, double complex * G, double complex * A){

//  Calculates the four-vector field from the links in Uaux,
//  which are the gauge-transformed links from the gauge-fixing
//  routines, and returns it in A

//  The formula is A = (U_aux - U_aux_dagger)/2i - trace_part

	pos_vec position;

	double complex * g_dagger_positionplusmu = (double complex *) malloc(3 * 3 * sizeof(double complex));
	
	double complex * U_aux = (double complex *) malloc(3 * 3 * sizeof(double complex));
    double complex * U_aux_dagger = (double complex *) malloc(3 * 3 * sizeof(double complex));
	double complex * U_minus_Udagger = (double complex *) malloc(3 * 3 * sizeof(double complex));

	double complex * trace_part = (double complex *) malloc(3 * 3 * sizeof(double complex));

	for ( position.t = 0; position.t < Nt; position.t++) {
		for ( position.i = 0; position.i < Nxyz; position.i++) {
			for ( position.j = 0; position.j < Nxyz; position.j++) {
				for ( position.k = 0; position.k < Nxyz; position.k++) {
					for (int mu = 0; mu < d; mu++) {
		
						//	U_aux_mu(x)=g(x)U_mu(x)g_dagger(x+mu)
						SU3_hermitean_conjugate(get_gauge_transf(G, hop_position_positive(position, mu)), g_dagger_positionplusmu);
						SU3_product_three(get_gauge_transf(G, position), get_link(U, position, mu), g_dagger_positionplusmu, U_aux);
						SU3_hermitean_conjugate(U_aux, U_aux_dagger);
						
						//	(U_aux - U_aux_dagger)/2i 
						SU3_subtraction(U_aux, U_aux_dagger, U_minus_Udagger);
						SU3_multiplication_by_scalar(U_minus_Udagger, -0.5 * I , get_link(A, position, mu));

						//  Subtract the trace part

						SU3_multiplication_by_scalar(SU3_identity, - (1.0 / 3.0) * SU3_trace(get_link(A,position,mu)) , trace_part);

						SU3_accumulate(trace_part, get_link(A, position, mu));

					}
				}
			}	
		}
	}

	free(g_dagger_positionplusmu);
	free(U_aux); free(U_aux_dagger);
	free(U_minus_Udagger);
	free(trace_part);
}

void SU3_divergence_A(double complex * A, pos_vec position, double complex * div_A){	
    
    //	Calculates the divergence of the field A on the lattice
    //  and returns it in div_A.

	double complex * term_divA = (double complex *) malloc(3 * 3 * sizeof(double complex));  //  Each term in the sum over all directions

	SU3_copy(SU3_null, div_A);
	for (int mu = 0; mu < d; mu++) {

		SU3_subtraction(get_link(A,position,mu), get_link(A, hop_position_negative(position, mu),mu), term_divA);

		SU3_accumulate(term_divA, div_A);	//	Sum of terms over all directions mu.

	}

	free(term_divA);
}