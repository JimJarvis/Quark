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

TEST(Qumat, Eigen_Kronecker_Hadamard)
{
	for (int nqubit : QubitRange(1))
		ASSERT_MAT(
		hadamard_recursive(nqubit),
		hadamard_mat(nqubit), _S + "Fail at qubit" + nqubit);
}

TEST(Qumat, Qureg_Kronecker)
{
	for (int nqubit : QubitRange(1))
	{
		Qureg qd1 = rand_qureg_dense(nqubit, 1);
		Qureg qd2 = rand_qureg_dense(nqubit / 2 + 1, 2);

		size_t sparse = half_fill(nqubit);
		Qureg qs1 = rand_qureg_sparse(nqubit, sparse, 1, false);
		Qureg qs2 = rand_qureg_sparse(nqubit, sparse/2 + 1, 2);

		Qureg QQs[] = { qd1, qd2, qs1, qs2 };

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