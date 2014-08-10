#include "tests.h"

TEST(Quop, Qureg_Hadamard)
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