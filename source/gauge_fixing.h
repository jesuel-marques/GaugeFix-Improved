#ifndef GAUGEFIXING_H
#define GAUGEFIXING_H

void SU3_calculate_w(double complex * U, double complex * G, pos_vec position, double complex * w);
// void SU3_update_U_aux(double complex *U, double complex * G, double complex *U_aux);
// void SU3_local_update_U_aux(pos_vec position);
void SU3_update_U(double complex * U, double complex * G);

double SU3_calculate_e2_local(double complex *U, pos_vec position);
double SU3_calculate_e2(double complex * U);

void update_sub_LosAlamos(double complex * matrix_SU3, char submatrix, double * update_SU3);

void SU3_LosAlamos_common_block(double complex * w, double complex * A);
void SU3_gaugefixing_overrelaxation(double complex * U, double complex * G, pos_vec position);
int SU3_gauge_fix(double complex * U, double complex * G, int config);

// void SU3_print_gaugefixed_U(double complex * U, double complex * G, double complex *U_aux, char filename[max_length_name]);

#endif