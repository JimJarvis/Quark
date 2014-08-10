#include "tests.h"

// Qubit from 1 to 7
auto qubitRange = Range<>(1, 8);

auto qubase_range = [](int nqubit) { return Range<qubase>(1 << nqubit); };

// An inefficient recursive version of hadarmard
MatrixXcf hadamard_recursive(int nqubit)
{
	if (nqubit == 1)
		return hadamard_mat();
	else
		return kronecker_mat(
			hadamard_mat(), hadamard_recursive(nqubit - 1));
}

TEST(Quop, Eigen_Kronecker_Hadamard)
{
	for (int nqubit : qubitRange)
		ASSERT_MAT(
		hadamard_recursive(nqubit),
		hadamard_mat(nqubit), _S + "Fail at qubit" + nqubit);
}

TEST(Quop, Qureg_Hadamard)
{
	Qureg qq;
	for (int nqubit : qubitRange)
	{
		MatrixXcf gold = hadamard_mat(nqubit);
		// Each column should agree with hadamard mat
		for (qubase base : qubase_range(nqubit))
		{
			Qureg QQs[2] = 
			{
				Qureg::create<true>(nqubit, base), 
				Qureg::create<false>(nqubit, 1 << nqubit, base)
			};
			for (Qureg qq : QQs)
			{
				hadamard(qq);
				ASSERT_MAT(gold.col(base), VectorXcf(qq),
						   _S + (qq.dense ? "Dense" : "Sparse")
						   + " disagree at base " + bits2str(base, nqubit));
			}
		}
	}
}

TEST(Quop, Qureg_Kronecker)
{
	for (int nqubit : qubitRange)
	{
		Qureg qd1 = Qureg::create<true>(nqubit);
		for (qubase base : qd1.base_iter_d())
			qd1.set_base_d(base, rand_cx(1));

		Qureg qd2 = Qureg::create<true>(nqubit / 2 + 1);
		for (qubase base : qd2.base_iter_d())
			qd2.set_base_d(base, rand_cx(2));

		size_t total = 1 << nqubit;
		size_t sparseCap = total / 2 + 1;
		Qureg qs1 = Qureg::create<false>(nqubit, sparseCap);
		for (qubase base : Range<qubase, false>(total, total - sparseCap))
			qs1.add_base(base, rand_cx(1));

		Qureg qs2 = Qureg::create<false>(nqubit, sparseCap);
		for (qubase base : Range<qubase>(sparseCap / 2 + 1))
			qs2.add_base(base, rand_cx(2));

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
				kronecker_mat(vec1, vec2), VectorXcf(qvecProd), "");
		}
	}
}