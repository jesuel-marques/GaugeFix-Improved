#include <stdio.h>							//	Standard header C files
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "../SU3_parameters.h"				//	Simulation parameters
#include "../SU3_gaugefixing_parameters.h"	//	Gauge-fixing specific parameters

#include "lattice.h"						//	Initialization functions and 
											//	calculations of positions and
											//	links on the lattice.


#include "math_ops.h"						//	Math operations
#include "SU2_ops.h"						//	SU(2) operations
#include "SU3_ops.h"						//	SU(3) operations
#include "fourvector_field.h"				//	Calculation of A_mu(n) and related things

void SU3_calculate_w(double complex * U, double complex * G, pos_vec position, double complex * w){	

	//	Calculates w(n) = g(n).h(n), following the notation in hep-lat/0301019v2
	//	returns result in w.
	
	double complex * h = (double complex *) malloc(3 * 3 * sizeof(double complex));;	//	h(n)
	//	h(n) = sum_mu U_mu(n).g_dagger(n+mu_hat)+U_dagger_mu(n-mu_hat).g_dagger(n-mu_hat)


	double complex * prod_front = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	For each mu there will be 
	double complex * prod_rear = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	two matrix products to be performed.

	if ( h == NULL || prod_front == NULL || prod_rear == NULL ) {
		//	Test if allocation was successful.
        printf("Memory allocation failed in function SU3_calculate_w");
        exit(1); 

    }

								//	Relevant positions to the 
	pos_vec position_minus_mu;	//	calculation of h(n).

	double complex * g_dagger_front = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	Matrices necessary to the calculation
	double complex * g_dagger_rear = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	of h(n).

	double complex * u_dagger_rear = (double complex *) malloc(3 * 3 * sizeof(double complex));	

	if ( g_dagger_front == NULL || g_dagger_rear == NULL || u_dagger_rear == NULL) {
		//	Test if allocation was successful.
        printf("Memory allocation failed in function SU3_calculate_w");
        exit(1); 

    }

	SU3_set_to_null(h);		//	Initializing h(n)=0

	// h(n)	calculation

	for (int mu = 0; mu < d; mu++){
		//	h(n) = sum_mu U_mu(n).g_dagger(n+mu_hat)+U_dagger_mu(n-mu_hat).g_dagger(n-mu_hat)

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


void SU3_update_U(double complex * U, double complex * G){
	
	//	Updates U_aux, which must be done before calculating A

	double complex * u_updated = (double complex *) malloc(3 * 3 * sizeof(double complex));;	//	Matrices needed to calculate U_aux.

	double complex * g_dagger_positionplusmu = (double complex *) malloc(3 * 3 * sizeof(double complex));;

	if ( u_updated == NULL || g_dagger_positionplusmu == NULL) {
		//	Test if allocation was successful.
        printf("Memory allocation failed in function SU3_update_U");
        exit(1); 

    }

	pos_vec position;					

	for (position.t = 0; position.t < Nt; position.t++){
		for (position.i = 0; position.i < Nxyz; position.i++){
			for (position.j = 0; position.j < Nxyz; position.j++){
				for (position.k = 0; position.k < Nxyz; position.k++){
					for (int mu = 0; mu < d; mu++){

						//	U'_mu(x)=g(x)U_mu(x)g_dagger(x+mu)

						SU3_hermitean_conjugate(get_gauge_transf(G,hop_position_positive(position,mu)), g_dagger_positionplusmu);						
						
						SU3_product_three(get_gauge_transf(G,position),get_link(U,position, mu) , g_dagger_positionplusmu, u_updated);
						SU3_copy(u_updated, get_link(U,position, mu));
					}
                }
            }
        }
    }

	free(u_updated);
	free(g_dagger_positionplusmu);
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

double SU3_calculate_e2_local(double complex *U, pos_vec position){
	double e2_local=0.0;
	
	double complex * div_A = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	Quadridivergence of A
	double * div_A_components = (double *) malloc(9 * sizeof(double));	//	Components of the quadridergence projected onto the Gell-Mann matrices
	

	if ( div_A == NULL || div_A_components == NULL ) {
		//	Test if allocation was successful.
        printf("Memory allocation failed in function SU3_calculate_e2_local");
        exit(1); 

    }
	//	Calculates quadridivergence and projects to the SU(3) algebra

	SU3_divergence_A(U, position, div_A);
	SU3_decompose_algebra(div_A, div_A_components);

	for(int a = 1; a <= 8; a++){	
		//	Normalized sum of the squares of the color components of the divergence of A.
		
		e2_local += (double)pow2(div_A_components[a])/(Volume);
	}
	free(div_A); free(div_A_components);

	return e2_local;
}

double SU3_calculate_e2(double complex * U){

//	Calculates e2 (defined in hep-lat/0301019v2), 
//	used to find out distance to the gauge-fixed situation.
	
	pos_vec position;

	double e2 = 0.0;

	for (position.t = 0; position.t < Nt; position.t++) {
		for (position.i = 0; position.i < Nxyz; position.i++) {
			for (position.j = 0; position.j < Nxyz; position.j++) {
				for (position.k = 0; position.k < Nxyz; position.k++) {

					e2 += SU3_calculate_e2_local(U, position);

				}
			}
		}
	}

	return e2;					
}

void update_sub_LosAlamos(double complex * matrix_SU3, char submatrix, complex double * update_SU3){

	int a, b;
	double * matrix_SU2 = (double *) malloc( 4 * sizeof(double));


	if ( matrix_SU2 == NULL ) {
		//	Test if allocation was successful.
        printf("Memory allocation failed in function project_SU3_to_sub_SU2");
        exit(1); 

    }

	SU3_set_to_null(update_SU3);

	switch(submatrix){
		case 'R':
			a = 0;
			b = 1;
			update_SU3[ 2 * 3 + 2] = 1.0;
			break;

		case 'S':
			a = 0;
			b = 2;
			update_SU3[ 1 * 3 + 1] = 1.0;
			break;

		case 'T':
			a = 1;
			b = 2;
			update_SU3[ 0 * 3 + 0] = 1.0;
			break;

		default:
			printf("Invalid submatrix");
			exit(1);
	}

	matrix_SU2[0] = (creal(matrix_SU3[a * 3 + a]) + creal(matrix_SU3[b * 3 + b]));
	matrix_SU2[1] = (cimag(matrix_SU3[a * 3 + b]) + cimag(matrix_SU3[b * 3 + a]));
	matrix_SU2[2] = (creal(matrix_SU3[a * 3 + b]) - creal(matrix_SU3[b * 3 + a]));
	matrix_SU2[3] = (cimag(matrix_SU3[a * 3 + a]) - cimag(matrix_SU3[b * 3 + b]));

	SU2_projection(matrix_SU2);

	update_SU3[a * 3 + a] =  matrix_SU2[0] + I * matrix_SU2[3];
	update_SU3[a * 3 + b] =  matrix_SU2[2] + I * matrix_SU2[1];
	update_SU3[b * 3 + a] = -matrix_SU2[2] + I * matrix_SU2[1];
	update_SU3[b * 3 + b] =  matrix_SU2[0] - I * matrix_SU2[3];

	free(matrix_SU2);
}

void SU3_LosAlamos_common_block(double complex * w, double complex * A){

	//	Calculates the update matrix A from w(n)=g(n).h(n) as in the Los Alamos
	//	algorithm for SU(3), with a division of the update matrix in submatrices
	//	following the Cabbibo-Marinari trick. Actual A(n) is obtained after a number
	//	of "hits" to be performed one after another.

	double complex * update = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	First the conjugate to w, after that to Rw, then to SRw, ...

	double complex * Omega  = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	First the conjugate to w, after that to Rw, then to SRw, ...
	
	char submatrices[3] = {'R','S','T'};

	if ( update == NULL || Omega == NULL ) {
		//	Test if allocation was successful.
        printf("Memory allocation failed in function SU3_LosAlamos_common_block");
        exit(1); 

    }

	double complex * Aw   = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	Accumulated update for w(x)

	if ( Aw == NULL ) {
		//	Test if allocation was successful.
        printf("Memory allocation failed in function SU3_LosAlamos_common_block");
        exit(1); 

    }

	SU3_copy(w, Aw);

	for (int hits = 1; hits <= maxhits; hits++) {	
		//	Each hit contains the Cabbibo-Marinari subdivision
		for (int i = 0; i <= 2 ; i++){
			//	Submatrices are indicated by letter in subs

			SU3_hermitean_conjugate(Aw, Omega);
			update_sub_LosAlamos(Omega, submatrices[i], update);
			SU3_accumulate_left_product(update, Aw);

		}
	}
	free(Omega); free(update);

	double complex * w_inverse = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	Inverse of w(x)

	if ( w_inverse == NULL ) {
		//	Test if allocation was successful.
        printf("Memory allocation failed in function SU3_LosAlamos_common_block");
        exit(1); 

    }

	inverse_3_by_3(w, w_inverse);
	SU3_product(Aw, w_inverse, A);	//	Updates matrix to g. It is the
									//	accumulated updates from the hits.
									//	One has to multiply by w_inverse to 
									// 	eliminate w of the product.

	free(Aw);
	free(w_inverse); 	
}

void SU3_gaugefixing_overrelaxation(double complex * U, double complex * G, pos_vec position){	

//	Generalization of the algorithm described in hep-lat/0301019v2, using the 
//	Cabbibo-Marinari submatrices trick.
//	It updates the g at the given position.

	double complex * A = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	Final update to g coming from Los Alamos block.
						//	g_updated=A*g_old

	double complex * w  = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	w(n)=g(n)*h(n), following the notation of hep-lat/0301019v2

	if ( A == NULL || w == NULL ) {
		//	Test if allocation was successful.
        printf("Memory allocation failed in function SU3_gaugefixing_overrelaxation");
        exit(1); 

    }

	SU3_calculate_w(U, G, position, w);	//	Calculating w(n)=g(n)*h(n)
	SU3_LosAlamos_common_block(w, A);

	free(w);

	//	The above function determines A which would be the naïve
	//	update to bring the local function to its mininum. However
	//	we found out that actually using overrelaxation, which 
	//	means using A^omega instead of A, where 1<omega<2 makes the
	//	converge faster. A^omega is calculated using the first two terms
	//	of the binomial expansion: A^omega=I+omega(A-I)+...=I(1-omega)+omega*A+...


	double complex * omega_ORA = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	omega_OR times A
	double complex * A_OR = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	The over-relaxation update.

	if ( omega_ORA == NULL || A_OR == NULL ) {
		//	Test if allocation was successful.
        printf("Memory allocation failed in function SU3_gaugefixing_overrelaxation");
        exit(1); 

    }

	SU3_set_to_identity(A_OR);
	SU3_substitution_multiplication_by_scalar(1.0 - omega_OR, A_OR);
	SU3_multiplication_by_scalar(omega_OR, A, omega_ORA);
	SU3_accumulate(omega_ORA, A_OR);
	
	free(A); free(omega_ORA); 

	double complex * g_old;
	double complex * g_updated = (double complex *) malloc(3 * 3 * sizeof(double complex));	//	The updated value to be stored after
								//	projecting to SU(3)

	if ( g_updated == NULL ) {
		//	Test if allocation was successful.
        printf("Memory allocation failed in function SU3_gaugefixing_overrelaxation");
        exit(1); 

    }

	g_old = get_gauge_transf(G, position);
	SU3_product(A_OR, g_old, g_updated);
	SU3_projection(g_updated);
	SU3_copy(g_updated, g_old);

	free(A_OR);
	free(g_updated);
}

int SU3_gauge_fix(double complex * U, double complex * G, int config){
	//	Fix the gauge and follows the process by calculating e2;

	double e2 = SU3_calculate_e2(U);	
	//	Gauge-fixing index, which will be less than the tolerance,
	//	when the gauge-fixing is considered to be attained.
	//	Following the notation of hep-lat/0301019v2


	pos_vec position;	//	Position vector, to be looped over
	
	int	cont = 0;	
	//	Counter to the number of sweeps to fix config to Landau gauge

	do {	
		for (position.t = 0; position.t < Nt; position.t++) {
			for (position.i = 0; position.i < Nxyz; position.i++) {
				for (position.j = 0; position.j < Nxyz; position.j++) {
					for (position.k = 0; position.k < Nxyz; position.k++) {	
	
						SU3_gaugefixing_overrelaxation(U, G, position);	//	The actual gauge-fixing algorithm

					}
                }
            }
        }

		cont++;


		if (cont % sweeps_to_measurement_e2 == 0) {	//	No need to calculate e2 all the time
													//	because it will take some hundreds of sweeps
													//	to fix the gauge.
			
			SU3_update_U(U,G);						//	Updates the links with the accumulated 
													//	transformation and then set the transformation to
			SU3_set_gauge_transf_to_identity(G);	//	identity again.
							
			e2 = SU3_calculate_e2(U);	//	Indicates how far we are to the Landau-gauge 
										//	condition, e2=0. We have to update U before calculating this.
			
			printf("\nconfig: %d, sweep: %d, e2: %3.2e\n", config, cont, e2);

		}
		
	} while (e2 > tolerance);	//	As long as e2 is larger than the tolerance
								//	repeats the process iteratively.
	
							
	printf("Sweeps needed to gauge-fix: %d \n", cont);

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