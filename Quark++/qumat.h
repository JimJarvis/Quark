#ifndef qumat_h__
#define qumat_h__

#include "qureg.h"
namespace Qumat // quantum matrices
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

	void normalize(Q);

	/**********************************************/
	/*********** Eigen operations  ***********/
	/**********************************************/
	MatrixXcf hadamard_mat(int nqubit);

	Matrix2cf hadamard_mat();

	MatrixXcf kronecker_mat(const MatrixXcf& A, const MatrixXcf& B);

	Matrix4cf cnot_mat();

	Matrix<CX, 8, 8> toffoli_mat();
}

#endif // qumat_h__
