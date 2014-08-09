#include "tests.h"

// An inefficient recursive version of hadarmard
MatrixXcf hadamard_recursive(int nqubit)
{
	static CX _sqrt2 = CX(1 / sqrt(2));
	static Matrix2cf HadamardMat;
	HadamardMat <<
		_sqrt2, _sqrt2,
		_sqrt2, -_sqrt2;
	if (nqubit == 1)
		return HadamardMat;
	else
		return kronecker_mat(
			HadamardMat, hadamard_recursive(nqubit - 1));
}

TEST(Quop, Kronecker_Hadamard)
{
	for (int i = 1; i < 8; ++i)
		ASSERT_MAT(
			hadamard_recursive(i), 
			hadamard_mat(i)) << "Fail at qubit " << i;
}

TEST(Quop, Hadamard_gate)
{
	Qureg qq;
	for (int nqubit : Range<>(1, 8))
	{
		MatrixXcf gold = hadamard_mat(nqubit);
		// Each column should agree with hadamard mat
		for (qubase base : Range<qubase>(1 << nqubit))
		{
			Qureg qqds[2] = 
			{
				Qureg::create<true>(nqubit, base), 
				Qureg::create<false>(nqubit, 1 << nqubit, base)
			};
			for (Qureg qq : qqds)
			{
				hadamard(qq);
				ASSERT_MAT(gold.col(base), VectorXcf(qq))
					<< (qq.dense ? "Dense" : "Sparse") 
					<< " disagree at base " << bits2str(base, nqubit);
			}
		}
	}
}