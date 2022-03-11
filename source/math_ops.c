#include "../SU3_parameters.h"

#include "SU3_global.h"			//	Definição de variáveis globais

#include "SU3_ops.h"

double pow2(double x){ 

//	Calculates the square of a number

	return x * x;
}

void inverse_3_by_3(double complex * x, double complex * x_inverse){

	//	Calculates the inverse of a 3 by 3 matrix x explicitely
    //  and returns result in x_inverse.

	double complex x_det;
	
	x_det = SU3_determinant(x);

	x_inverse[0 * 3 + 0] = (*(x + 1 * 3 + 1) * *(x + 2 * 3 + 2) - *(x + 1 * 3 + 2) * *(x + 2 * 3 + 1)) / x_det;
	x_inverse[0 * 3 + 1] = (*(x + 0 * 3 + 2) * *(x + 2 * 3 + 1) - *(x + 0 * 3 + 1) * *(x + 2 * 3 + 2)) / x_det;
	x_inverse[0 * 3 + 2] = (*(x + 0 * 3 + 1) * *(x + 1 * 3 + 2) - *(x + 0 * 3 + 2) * *(x + 1 * 3 + 1)) / x_det;

	x_inverse[1 * 3 + 0] = (*(x + 1 * 3 + 2) * *(x + 2 * 3 + 0) - *(x + 1 * 3 + 0) * *(x + 2 * 3 + 2)) / x_det;
	x_inverse[1 * 3 + 1] = (*(x + 0 * 3 + 0) * *(x + 2 * 3 + 2) - *(x + 0 * 3 + 2) * *(x + 2 * 3 + 0)) / x_det;
	x_inverse[1 * 3 + 2] = (*(x + 0 * 3 + 2) * *(x + 1 * 3 + 0) - *(x + 0 * 3 + 0) * *(x + 1 * 3 + 2)) / x_det;

	x_inverse[2 * 3 + 0] = (*(x + 1 * 3 + 0) * *(x + 2 * 3 + 1) - *(x + 1 * 3 + 1) * *(x + 2 * 3 + 0)) / x_det;
	x_inverse[2 * 3 + 1] = (*(x + 0 * 3 + 1) * *(x + 2 * 3 + 0) - *(x + 0 * 3 + 0) * *(x + 2 * 3 + 1)) / x_det;
	x_inverse[2 * 3 + 2] = (*(x + 0 * 3 + 0) * *(x + 1 * 3 + 1) - *(x + 0 * 3 + 1) * *(x + 1 * 3 + 0)) / x_det;
}