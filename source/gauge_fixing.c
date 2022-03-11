#include <stdio.h>						//	Standard header C files
#include <stdlib.h>
#include <math.h>

#include "../SU3_parameters.h"			//	Simulation parameters
#include "../SU3_gaugefixing_parameters.h"

#include "SU3_global.h"					//	Global variables definitions

#include "lattice.h"

#include "math_ops.h"
#include "SU2_ops.h"
#include "SU3_ops.h"
#include "fourvector_field.h"


void RST_initialize(double complex * R, double complex * S, double complex * T){

	for ( int e = 0; e < 3; e++) {
		for (int f = 0; f < 3; f++) {

			if (e == 2 && f == 2) {
				R[e * 3 + f] = 1.0;
			}
			else{
				R[e * 3 + f] = 0.0;
			}
			if (e == 1 && f == 1) {
				S[e * 3 + f] = 1.0;
			}
			else{
				S[e * 3 + f] = 0.0;
			}
			if (e == 0 && f == 0) {
				T[e * 3 + f] = 1.0;
			}
			else{
				T[e * 3 + f] = 0.0;
			}
		}
	}
}

void SU3_calculate_w(double complex * U, double complex * G, pos_vec position, double complex * w){	

	//	Calculates w(n) = g(n).h(n), following the notation in hep-lat/0301019v2
	//	returns result in w.
	//	h(n) = sum_mu U_mu(n).g_dagger(n+mu_hat)+U_dagger_mu(n-mu_hat).g_dagger(n-mu_hat)
	
	double complex * h = (double complex *) malloc(3 * 3 * sizeof(double complex));;	//	h(n)

	double complex * prod_front = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	For each mu there will be 
	double complex * prod_rear = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	two matrix products to be performed.

								//	Relevant positions to the 
	pos_vec position_minus_mu;	//	calculation of h(n).

	double complex * g_dagger_front = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	Matrices necessary to the calculation
	double complex * g_dagger_rear = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	of h(n).

	double complex * u_dagger_rear = (double complex *) malloc(3 * 3 * sizeof(double complex));	


	SU3_copy(SU3_null, h);		//	Initializing h(n)=0

	// h(n)	calculation

	for (int mu = 0; mu < d; mu++){

		position_minus_mu = hop_position_negative(position, mu);

		SU3_hermitean_conjugate(get_gauge_transf(G,hop_position_positive(position, mu)), g_dagger_front);
		SU3_hermitean_conjugate(get_gauge_transf(G, position_minus_mu), g_dagger_rear);

		SU3_hermitean_conjugate(get_link(U, position_minus_mu, mu), u_dagger_rear);

		SU3_product(get_link(U,position,mu), g_dagger_front, prod_front);
		SU3_accumulate(prod_front, h);

		SU3_product(u_dagger_rear, g_dagger_rear, prod_rear);
		SU3_accumulate(prod_rear, h);

	}
	//	w(x) = g(x).h(x)

	SU3_product(get_gauge_transf(G, position), h, w);

	free(h);

	free(prod_front); free(prod_rear);
	
	free(g_dagger_front); free(g_dagger_rear);
	
	free(u_dagger_rear);
}

// void SU3_local_update_U_aux(pos_vec position){

// 	//	Updates U_aux only at a given position

// 	double complex u_updated[3][3];	//	Matrices and positions needed to calculate U_aux.
	
// 	int versor_mu[d] = { 0, 0, 0, 0 };
// 	pos_vec position_plus_mu;
// 	double complex g_dagger_positionplusmu[3][3];

// 	int versor_minus_mu[d] = { 0, 0, 0, 0 };
// 	pos_vec position_minus_mu;
// 	double complex g_dagger[3][3];

// 	for (int mu = 0; mu < d ; mu++) {
 
// 		//	U'_mu(x)=g(x)U_mu(x)g_dagger(x+mu)

// 		versor_mu[mu] = 1;

// 		add_position_vector(position, versor_mu, position_plus_mu);

// 		SU3_hermitean_conjugate(G[position_plus_mu[0]][position_plus_mu[1]][position_plus_mu[2]][position_plus_mu[3]], g_dagger_positionplusmu);
		
// 		SU3_product_three(G[position.t][position.i][position.j][position.k], U[position.t][position.i][position.j][position.k][mu], g_dagger_positionplusmu, u_updated);

// 		//	Projetar no grupo SU(3) por segurança
// 		SU3_projection(u_updated, U_aux[position.t][position.i][position.j][position.k][mu]);

// 		versor_mu[mu] = 0; 

// 		//	U'_mu(x-mu)=g(x-mu)U_mu(x-mu)g_dagger(x)

// 		versor_minus_mu[mu] = -1;

// 		add_position_vector(position, versor_minus_mu, position_minus_mu);

// 		SU3_hermitean_conjugate(G[position.t][position.i][position.j][position.k], g_dagger);

// 		SU3_product_three(G[position_minus_mu[0]][position_minus_mu[1]][position_minus_mu[2]][position_minus_mu[3]],U[position_minus_mu[0]][position_minus_mu[1]][position_minus_mu[2]][position_minus_mu[3]][mu], g_dagger, u_updated);
		
// 		//	Projetar no grupo SU(3) por segurança
// 		SU3_projection(u_updated, U_aux[position_minus_mu[0]][position_minus_mu[1]][position_minus_mu[2]][position_minus_mu[3]][mu]);
		
// 		versor_minus_mu[mu] = 0; 

// 	}	
// }

double SU3_calculate_e2(double complex * A){

//	Função calcula e2 (definido na Dissertação Nelson, eq. 4.53, e hep-lat/0301019v2), 
//	utilizado para definir distancia até situação de gauge fixo.
	
	pos_vec position;

	double complex * div_A = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	Quadridivergência de A
	double * div_A_components = (double *) malloc(9 * sizeof(double));	//	Componentes da quadridivergência projetada nas matrizes de Gell-Mann

	double e2 = 0.0;

	for (position.t = 0; position.t < Nt; position.t++) {
		for (position.i = 0; position.i < Nxyz; position.i++) {
			for (position.j = 0; position.j < Nxyz; position.j++) {
				for (position.k = 0; position.k < Nxyz; position.k++) {

					//	Calcula quadridivergencia e decompõe na álgebra de SU(3)

					SU3_divergence_A(A, position, div_A);
					SU3_decompose_algebra(div_A, div_A_components);
					for(int a = 1; a <= 8; a++)	//	Soma normalizada dos quadrados das componentes de cor da divergência de A
						e2 += (double)pow2(div_A_components[a])/(Volume);

				}
			}
		}
	}

	free(div_A); free(div_A_components);
	return e2;					
}


void SU3_LosAlamos_common_block(double complex * w, double complex * A, double complex *R, double complex *S, double complex *T){

	//	Calculates the update matrix A from w(n)=g(n).h(n) as in the Los Alamos
	//	algorithm for SU(3), with a division of the update matrix in submatrices
	//	following the Cabbibo-Marinari trick. Actual A(n) is obtained after a number
	//	of "hits" to be performed one after another.



	double complex * Omega = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	First the conjugate to w, after that to Rw, then to SRw, ...
	double complex * Omega_dagger = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	Hermitean conjugate to Omega 
	
	double complex * w_inverse = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	Inverse of w(x)
	double complex * Rw = (double complex *) malloc(3 * 3 * sizeof(double complex));		//	First update to w(x)
	double complex * SRw = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	Second update to w(x)	
	double complex * TSRw = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	Third update to w(x)


	//	Omega matrices and projections to SU(2);

	double * w_R = (double *) malloc( 4 * sizeof(double));
	double * w_RSU2 = (double *) malloc( 4 * sizeof(double));

	double * w_S = (double *) malloc( 4 * sizeof(double));
	double * w_SSU2 = (double *) malloc( 4 * sizeof(double));

	double * w_T = (double *) malloc( 4 * sizeof(double)); 
	double * w_TSU2 = (double *) malloc( 4 * sizeof(double));


	SU3_hermitean_conjugate(w, Omega);


	for (int hits = 1; hits <= maxhits; hits++) {	
		//	Each hit contains the Cabbibo-Marinari subdivision

		SU3_hermitean_conjugate(Omega, Omega_dagger); 

		//	First submatrix

		w_R[0] = (creal(*(Omega + 0 * 3 + 0)) + creal(*(Omega + 1 * 3 + 1)));
		w_R[1] = (cimag(*(Omega + 0 * 3 + 1)) + cimag(*(Omega + 1 * 3 + 0)));
		w_R[2] = (creal(*(Omega + 0 * 3 + 1)) - creal(*(Omega + 1 * 3 + 0)));
		w_R[3] = (cimag(*(Omega + 0 * 3 + 0)) - cimag(*(Omega + 1 * 3 + 1)));


		SU2_projection(w_R, w_RSU2);
		
		R[0 * 3 + 0] = w_RSU2[0] + I * w_RSU2[3];
		R[0 * 3 + 1] = w_RSU2[2] + I * w_RSU2[1];
		R[1 * 3 + 0] = -w_RSU2[2] + I * w_RSU2[1];
		R[1 * 3 + 1] = w_RSU2[0] - I * w_RSU2[3];


		SU3_product(R, Omega_dagger, Rw);
		SU3_hermitean_conjugate(Rw, Omega);


		//	Second submatrix

		w_S[0] = (creal(*(Omega + 0 * 3 + 0)) + creal(*(Omega + 2 * 3 + 2)));
		w_S[1] = (cimag(*(Omega + 0 * 3 + 2)) + cimag(*(Omega + 2 * 3 + 0)));
		w_S[2] = (creal(*(Omega + 0 * 3 + 2)) - creal(*(Omega + 2 * 3 + 0)));
		w_S[3] = (cimag(*(Omega + 0 * 3 + 0)) - cimag(*(Omega + 2 * 3 + 2)));


		SU2_projection(w_S, w_SSU2);

		S[0 * 3 + 0] = w_SSU2[0] + I * w_SSU2[3];
		S[0 * 3 + 2] = w_SSU2[2] + I * w_SSU2[1];
		S[2 * 3 + 0] = -w_SSU2[2] + I * w_SSU2[1];
		S[2 * 3 + 2] = w_SSU2[0] - I * w_SSU2[3];


	    SU3_product(S, Rw, SRw);
	    SU3_hermitean_conjugate(SRw, Omega);


	    //	Third submatrix

		w_T[0] = (creal(*(Omega + 1 * 3 + 1)) + creal(*(Omega + 2 * 3 + 2)));
		w_T[1] = (cimag(*(Omega + 1 * 3 + 2)) + cimag(*(Omega + 2 * 3 + 1)));
		w_T[2] = (creal(*(Omega + 1 * 3 + 2)) - creal(*(Omega + 2 * 3 + 1)));
		w_T[3] = (cimag(*(Omega + 1 * 3 + 1)) - cimag(*(Omega + 2 * 3 + 2)));

		SU2_projection(w_T, w_TSU2);

		T[1 * 3 + 1] = w_TSU2[0] + I * w_TSU2[3];
		T[1 * 3 + 2] = w_TSU2[2] + I * w_TSU2[1];
		T[2 * 3 + 1] = -w_TSU2[2] + I * w_TSU2[1];
		T[2 * 3 + 2] = w_TSU2[0] - I * w_TSU2[3];

		SU3_product(T, SRw, TSRw);
		SU3_hermitean_conjugate(TSRw, Omega);

	}
	inverse_3_by_3(w, w_inverse);
	SU3_product(TSRw, w_inverse, A);	//	Update matrix to g. It is the
										//	accumulated updates from the hits.
										//	One has to multiply by w_inverse to 
										// 	eliminate w of the product.

										//	Think if there is a way to get directly
										//  the update to g and not have to do this.
}

void SU3_gaugefixing_overrelaxation(double complex * G, pos_vec position, double complex * w, double complex *R, double complex *S, double complex *T){	

//	Generalization of the algorithm described in hep-lat/0301019v2, using the 
//	Cabbibo-Marinari submatrices trick.
//	It updates the g at the given position.

	double complex * A = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	Final update to g coming from Los Alamos block.
						//	g_updated=A*g_old

	double complex factor1 = 1.0-omega_OR;	//	Factors used in the binomial
	double complex factor2 = omega_OR;		//	truncated expansion.

	double complex * omega_ORA = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	omega_OR times A
	double complex * A_OR = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	The over-relaxation update.
	double complex * g_updated = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	The updated value to be stored after
								//	projecting to SU(3)

	double complex * g_at_position;

	SU3_LosAlamos_common_block(w, A, R, S, T);

	//	The above function determines A which would be the naïve
	//	update to bring the local function to its mininum. However
	//	we found out that actually using overrelaxation, which 
	//	means using A^omega instead of A, where 1<omega<2 makes the
	//	converge faster. A^omega is calculated using the first two terms
	//	of the binomial expansion: A^omega=I+omega(A-I)+...=I(1-omega)+omega*A+...

	SU3_multiplication_by_scalar(SU3_identity, factor1, A_OR);	//	This is a constant, no need to calculate all the time
	SU3_multiplication_by_scalar(A, factor2, omega_ORA);

	g_at_position = get_gauge_transf(G, position);
	SU3_accumulate(omega_ORA, A_OR);
	SU3_product(A_OR, g_at_position, g_updated);
	SU3_projection(g_updated, g_at_position);

	free(A);

	free(omega_ORA); free(A_OR); free(g_updated);
}

int SU3_gauge_fix(double complex * U, double complex * G, double complex * A){	//	Fixa o calibre de Landau

	pos_vec position;	//	Position vector, to be looped over

	int	cont = 0;	//	Counter to the number of sweeps to fix config to Landau gauge
	double e2 = SU3_calculate_e2(A);;	//	Gauge-fixing index, which will be less than the tolerance,
				//	when the gauge-fixing is considered to be attained.
				//	Following the notation of hep-lat/0301019v2

	double complex * w  = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	w(n)=g(n)*h(n), following the notation of hep-lat/0301019v2


	double complex * R = (double complex *) malloc(3 * 3 * sizeof(double complex));
	double complex * S = (double complex *) malloc(3 * 3 * sizeof(double complex));
	double complex * T = (double complex *) malloc(3 * 3 * sizeof(double complex));

	if ( R == NULL || S == NULL || T == NULL) {
        printf("Memory allocation failed for Cabbibo-Marinari matrices");
        exit(1); 
    }

	RST_initialize(R, S, T);
	//	Measuring the starting distance to the gauge-fixing condition, e2=0

	do {	
		for (position.t = 0; position.t < Nt; position.t++) {
			for (position.i = 0; position.i < Nxyz; position.i++) {
				for (position.j = 0; position.j < Nxyz; position.j++) {
					for (position.k = 0; position.k < Nxyz; position.k++) {	
							// getchar();		
							SU3_calculate_w(U, G, position, w);	//	Calculating w(n)=g(n)*h(n)

							SU3_gaugefixing_overrelaxation(G, position, w, R, S, T);	//	The actual gauge-fixing algorithm

					}
                }
            }
        }

		cont++;


		if (cont % sweeps_to_measurement_e2 == 0) {	//	No need to calculate e2 all the time
													//	because it will take some hundreds of sweeps
													//	to fix the gauge.
			SU3_update_A(U, G, A);
			e2 = SU3_calculate_e2(A);	//	Indicates how far we are to the Landau-gauge 
										//	condition, e2=0
			
			printf("\n%d e2: %3.2e\n",cont,e2);

		}
		
	} while (e2 > tolerance);	//	As long as e2 is larger than the tolerance
								//	repeats the process iteratively.
	
							
	printf("Sweeps needed to gauge-fix: %d \n", cont);

	free(w);

	free(R); free(S); free(T);

	return cont;
}

//void SU3_print_gaugefixed_U(double complex * U, double complex * G, double complex *U_aux, char filename[max_length_name]){	

// //	Prints the current configuration to a file, after the process
// //	of gauge-fixing. The configuration printed is therefore the one
// //	stored in U_aux.

//	FILE *gauge_configuration_file;


//	SU3_update_U_aux(U, G, U_aux); 

//	gauge_configuration_file = fopen(filename, "w+");

//	fwrite(U_aux, sizeof(U_aux), 1, gauge_configuration_file);

//	fclose(gauge_configuration_file);
		
//}