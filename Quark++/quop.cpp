/**********************************************
* Quantum operations  *
**********************************************/
#include "quop.h"
#include "qureg.h"

using namespace Quop;

Qureg Quop::kronecker(Q1, Q2, bool resultDense)
{
	int new_nqubit = q1.nqubit + q2.nqubit;
	size_t size1 = q1.size();
	size_t size2 = q2.size();

	Qureg qans = resultDense ?
		Qureg::create<true>(new_nqubit) :
		Qureg::create<false>(new_nqubit, size1 * size2);

	// If both of them are dense, then the resultant is also dense
	int ri = 0;
	int nqubit2 = q2.nqubit;
	for (int i1 = 0; i1 < size1; ++i1)
		for (int i2 = 0; i2 < size2; ++i2)
		{
			qubase newBase = (q1.get_base_d_s(i1) << nqubit2) | q2.get_base_d_s(i2);
			CX newAmp = q1.amp[i1] * q2.amp[i2];
			if (resultDense)
				qans.amp[newBase] = newAmp;
			else
				qans.add_base(newBase, newAmp);
		}

	return qans;
}

Qureg Quop::operator*(Q1, Q2)
{
	return kronecker(q1, q2, q1.dense && q2.dense);
}


/**********************************************/
/*********** Eigen  ***********/
/**********************************************/
MatrixXcf Quop::hadamard_mat(int nqubit)
{
	size_t size = 1 << nqubit;
	MatrixXcf mat(size, size);
	float coef = 1.0 / sqrt(size);
	for (size_t i = 0; i < size ; ++i)
		for (size_t j = 0; j < size; ++j)
			// binary dot product
			mat(i, j) = CX((bitwise_dot(i, j) ? -1 : 1) * coef);
	return mat;
}

MatrixXcf Quop::kronecker_mat(MatrixXcf& A, MatrixXcf& B)
{
	return MatrixXcf(1, 1);
}