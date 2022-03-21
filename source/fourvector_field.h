#ifndef FOURVECTORFIELD_H
#define FOURVECTORFIELD_H

void SU3_update_A_local(double complex * U, double complex * G, double complex * A, pos_vec position, int mu);

void SU3_divergence_A(double complex * A, pos_vec position, double complex * div_A);

#endif