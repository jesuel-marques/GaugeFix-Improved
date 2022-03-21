#ifndef LATTICE_H
#define LATTICE_H

typedef struct {
	int t;
	int i, j, k;
} pos_vec;	//	struct for position vectors


pos_vec add_position_vector(pos_vec u, pos_vec v);

pos_vec add_pos_vec(pos_vec u, pos_vec v);

void print_pos_vec(pos_vec u);

pos_vec hop_position_positive(pos_vec u, int mu);

pos_vec hop_position_negative(pos_vec u, int mu);


double complex * get_link(double complex * U, pos_vec position, int mu);

void get_link_matrix(double complex * U, pos_vec position, int mu, int direction, double complex * u);

double complex * get_gauge_transf(double complex * G, pos_vec position);


void SU3_load_config(char filename[max_length_name], double complex *U);

void SU3_copy_config(double complex *U, double complex *U_copy);

void SU3_set_gauge_transf_to_identity(double complex * G);

#endif