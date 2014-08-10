#include "tests.h"

// Qubit from 1 to 7
auto qubitRange = Range<>(1, 8);

auto qubase_range = [](int nqubit) { return Range<qubase>(1 << nqubit); };

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

TEST(Quop, Eigen_Kronecker_Hadamard)
{
	for (int nqubit : qubitRange)
		ASSERT_MAT(
		hadamard_recursive(nqubit),
		hadamard_mat(nqubit), ErrSS << "Fail at qubit " << nqubit);
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
			Qureg qqds[2] = 
			{
				Qureg::create<true>(nqubit, base), 
				Qureg::create<false>(nqubit, 1 << nqubit, base)
			};
			for (Qureg qq : qqds)
			{
				hadamard(qq);
				ASSERT_MAT(gold.col(base), VectorXcf(qq),
						   ErrSS << (qq.dense ? "Dense" : "Sparse")
						   << " disagree at base " << bits2str(base, nqubit));
			}
		}
	}
}

TEST(Quop, Qureg_Kronecker)
{
	for (int nqubit : qubitRange)
		for (int trial = 0; trial < 20 ; ++trial)
		{
			Qureg qd = Qureg::create<true>(nqubit);
			for (qubase base : qd.base_iter_d())
			{

			}
		}
}