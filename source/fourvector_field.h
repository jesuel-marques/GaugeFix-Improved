#ifndef FOURVECTORFIELD_H
#define FOURVECTORFIELD_H

void SU3_update_A(double complex * U, double complex * G, double complex * A);
void SU3_divergence_A(double complex * A, pos_vec position, double complex * div_A);

#endif