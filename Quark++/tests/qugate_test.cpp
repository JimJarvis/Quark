#include "tests.h"

TEST(Qugate, Hadamard)
{
	Qureg qq;
	for (int nqubit : QubitRange)
	{
		MatrixXcf gold = hadamard_mat(nqubit);
		// Each column should agree with hadamard mat
		for (qubase base : QubaseRange(nqubit))
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

TEST(Qugate, Cnot)
{
	for (int nqubit : QubitRange)
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