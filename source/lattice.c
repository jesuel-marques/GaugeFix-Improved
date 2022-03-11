#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "../SU3_parameters.h"						//	Simulation parameters
#include "../SU3_gaugefixing_parameters.h"

#include "SU3_global.h"

#include "SU3_ops.h"

typedef struct {
	int t, i, j, k;
} pos_vec;

pos_vec add_position_vector(pos_vec u, pos_vec v){

//	Adds two position vectors v and u, taking into account the periodic boundary conditions,
//  and puts result in vplusu.

	pos_vec u_plus_v;

	if ((u.t + v.t) >= 0) {	//	time component
		u_plus_v.t = ((u.t + v.t) % Nt);
    }
	else {
		u_plus_v.t = (((u.t + v.t) % Nt + Nt) % Nt);
    }

	if ((u.i + v.i) >= 0) {	//	x component
		u_plus_v.i = ((u.i + v.i) % Nxyz);
    }
	else {
		u_plus_v.i = (((u.i + v.i) % Nxyz + Nxyz) % Nxyz);
    }

	if ((u.j + v.j) >= 0) {	//	y component
		u_plus_v.j = ((u.j + v.j) % Nxyz);
    }
	else {
		u_plus_v.j = (((u.j + v.j) % Nxyz + Nxyz) % Nxyz);
    }

	if ((u.k + v.k) >= 0) {	//	z component
		u_plus_v.k = ((u.k + v.k) % Nxyz);
    }
	else {
		u_plus_v.k = (((u.k + v.k) % Nxyz + Nxyz) % Nxyz);
    }

	return u_plus_v;
}

void print_pos_vec(pos_vec pos){

	printf("% d %d %d %d\n", pos.t, pos.i, pos.j, pos.k);

}

pos_vec hop_position_positive(pos_vec u, int mu){
	
	pos_vec u_plus_muhat;

	switch (mu){

		case 0 :
			u_plus_muhat.t = ((u.t + 1) % Nt);
			u_plus_muhat.i = u.i;
			u_plus_muhat.j = u.j;
			u_plus_muhat.k = u.k;
			break;
		
		case 1 :
			u_plus_muhat.t = u.t;
			u_plus_muhat.i = ((u.i + 1) % Nxyz);
			u_plus_muhat.j = u.j;
			u_plus_muhat.k = u.k;
			break;
		
		case 2 :
			u_plus_muhat.t = u.t;
			u_plus_muhat.i = u.i;
			u_plus_muhat.j = ((u.j + 1) % Nxyz);
			u_plus_muhat.k = u.k;
			break;

		case 3 :
			u_plus_muhat.t = u.t;
			u_plus_muhat.i = u.i;
			u_plus_muhat.j = u.j;
			u_plus_muhat.k = ((u.k + 1) % Nxyz);
			break;

			default :
			printf("mu outside of range");
			exit(1);
			
	}

	return u_plus_muhat;
}

pos_vec hop_position_negative(pos_vec u, int mu){
	
	pos_vec u_minus_muhat;

	switch (mu) {

		case 0:
			u_minus_muhat.t = (((u.t - 1) % Nt + Nt) % Nt);
			u_minus_muhat.i = u.i;
			u_minus_muhat.j = u.j;
			u_minus_muhat.k = u.k;
			break;
		
		case 1:
			u_minus_muhat.t = u.t;
			u_minus_muhat.i = (((u.i - 1) % Nxyz + Nxyz) % Nxyz);
			u_minus_muhat.j = u.j;
			u_minus_muhat.k = u.k;
			break;
		
		case 2:
			u_minus_muhat.t = u.t;
			u_minus_muhat.i = u.i;
			u_minus_muhat.j = (((u.j - 1) % Nxyz + Nxyz) % Nxyz);
			u_minus_muhat.k = u.k;
			break;

		case 3:
			u_minus_muhat.t = u.t;
			u_minus_muhat.i = u.i;
			u_minus_muhat.j = u.j;
			u_minus_muhat.k = (((u.k - 1) % Nxyz + Nxyz) % Nxyz);
			break;

		default :
			printf("mu outside of range");
			exit(1);
	}

	return u_minus_muhat;
}

double complex * get_link(double complex * U, pos_vec position, int mu){
	//	Does the pointer arithmetic to get a pointer to link at given position and mu
	return U + ((((( position.t * Nxyz + position.i) * Nxyz + position.j) * Nxyz + position.k) * d + mu) * 3 * 3);
}


void get_link_matrix(double complex * U, pos_vec position, int mu, int direction, double complex * u){
	
	if (direction == FRONT) {

		SU3_copy(get_link(U, position, mu), u);
		//	Link in the positive way is what is stored in U

	}
	else if (direction == REAR) {

		SU3_hermitean_conjugate(get_link(U, hop_position_negative(position, mu), mu), u);
		//	U_(-mu)(n)=(U_mu(n-mu))^\dagger

	}
}

double complex * get_gauge_transf(double complex * G, pos_vec position){
	//	Does the pointer arithmetic to get a pointer to link at given position and mu
	return G + (((( position.t * Nxyz + position.i) * Nxyz + position.j) * Nxyz + position.k) * 3 * 3);
}

void SU3_initialize(){
	
	SU2_identity[0] = 1.0;
	SU2_identity[1] = 0.0;
	SU2_identity[2] = 0.0;
	SU2_identity[3] = 0.0;

	SU2_null[0] = 0.0;
	SU2_null[1] = 0.0;
	SU2_null[2] = 0.0;
	SU2_null[3] = 0.0;

	for ( int e = 0; e < 3; e++) {
		for (int f = 0; f < 3; f++) {
				
			SU3_null[e * 3 + f] = 0.0;
			
			if ( e != f ) {
				SU3_identity[e * 3 + f] = 0.0;
			}
			else{
				SU3_identity[e * 3 + f] = 1.0;
			}
			
		}
	}
}

void SU3_set_gauge_transf_to_identity(double complex * G) {

	pos_vec position;

	for ( position.t = 0; position.t < Nt; position.t++) {
		for ( position.i = 0; position.i < Nxyz; position.i++) {
			for ( position.j = 0; position.j < Nxyz; position.j++) {
				for ( position.k = 0; position.k < Nxyz; position.k++) {
					SU3_copy(SU3_identity, get_gauge_transf(G, position));
				}
			}
		}
	}

}

double SU3_load_U(char filename[max_length_name], double complex *U){	

//	Carrega a partir de um arquivo uma configuração de elos no programa

	FILE *config_file;
	
	// SU3_initialize();	//	Antes de carregar coloca matrizes identidade em toda a rede por segurança

	printf("Loading: %s\n", filename);

	
	config_file = fopen(filename, "r");

	if (fread(U, Volume * d * 3 * 3 * sizeof(double complex), 1, config_file) == 1) {
		printf("U Loaded\n");
	}
	else {
		printf(" Configuration loading failed.\n");
	}

	fclose(config_file);

}

double SU3_copy_config(double complex *U, double complex *U_aux){

	pos_vec position;

	for ( position.t = 0; position.t < Nt; position.t++) {
		for ( position.i = 0; position.i < Nxyz; position.i++) {
			for ( position.j = 0; position.j < Nxyz; position.j++) {
				for ( position.k = 0; position.k < Nxyz; position.k++) {
					for ( int mu = 0; mu < d; mu++){
						SU3_copy(get_link(U, position, mu), get_link(U_aux,position,mu));
					}
				}
			}
		}
	}

}