#ifndef quop_h__
#define quop_h__

#include "qureg.h"
namespace Quop // quantum operations
{
	/*
	 *	Kronecker product
	 * resultDense: true to return a dense Qureg
	 */
	Qureg kronecker(Q1, Q2, bool resultDense);

	/*
	 *	Kronecker product
	 * Result Qureg is dense only when both q1 and q2 are dense.
	 */
	Qureg operator*(Q1, Q2);

	CX dot(Q1, Q2);

	void normalize(Q);

	/**********************************************/
	/*********** Eigen operations  ***********/
	/**********************************************/
	MatrixXcf hadamard_mat(int nqubit);

	MatrixXcf kronecker_mat(MatrixXcf& A, MatrixXcf& B);
}

#endif // quop_h__
