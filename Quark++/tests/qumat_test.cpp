#include "tests.h"

// An inefficient recursive version of hadarmard
MatrixXcf hadamard_recursive(int nqubit)
{
	if (nqubit == 1)
		return hadamard_mat();
	else
		return kronecker_mat(
			hadamard_mat(), hadamard_recursive(nqubit - 1));
}

TEST(Qumat, EigenKroneckerHadamard)
{
	for (int nqubit : QubitRange(1))
		ASSERT_MAT(
			hadamard_recursive(nqubit),
			hadamard_mat(nqubit),
			_S + "Fail at qubit" + nqubit);
}

TEST(Qumat, QuregKronecker)
{
	for (int nqubit : QubitRange(1))
	{
		Qureg qd1 = rand_qureg_dense(nqubit, 1);
		Qureg qd2 = rand_qureg_dense(nqubit / 2 + 1, 2);

		size_t sparse = half_fill(nqubit);
		Qureg qs1 = rand_qureg_sparse(nqubit, sparse, 1, false);
		Qureg qs2 = rand_qureg_sparse(nqubit, sparse/2 + 1, 2);

		Qureg QQs[] = { move(qd1), move(qd2), move(qs1), move(qs2) };

		VectorXcf vec1, vec2, qvecProd;
		for (Qureg& q1 : QQs)
		for (Qureg& q2 : QQs)
		{
			if (&q1 == &q2) continue;

			vec1 = VectorXcf(q1);
			vec2 = VectorXcf(q2);
			qvecProd = kronecker(q1, q2, true);
			ASSERT_MAT(
				kronecker_mat(vec1, vec2), VectorXcf(qvecProd));
		}
	}
}

/*
 *	Shows that applying gates in any order shouldn't matter
 */
TEST(Qumat, KroneckerOrder)
{
	Matrix2cf m1 = Matrix2cf::Identity(2, 2);
	Matrix2cf m2 = rand_cxmat(2, 2);
	Matrix2cf m3 = rand_cxmat(2, 2);
	MatrixXcf res1 = (m1 & m1 & m3) * (m1 & m2 & m1);
	MatrixXcf res2 = (m1 & m2 & m1) * (m1 & m1 & m3);
	MatrixXcf res3 = m1 & (m2 & m3);
	ASSERT_MAT(res1, res2);
	ASSERT_MAT(res1, res3);
}