#include "tests.h"

TEST(Qugate, Hadamard)
{
	Qureg qq;
	for (int nqubit : QubitRange(1))
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
				ASSERT_MAT(gold.col(base), VectorXcf(qq));
						   //_S + (qq.dense ? "Dense" : "Sparse")
						   //+ " disagree at base " + bits2str(base, nqubit));
			}
		}
	}
}

TEST(Qugate, GenericGate1)
{
	Vector2cf oldBitAmp, newBitAmp;
	for (int nqubit : QubitRange(2))
	{
		Qureg qd = rand_qureg_dense(nqubit, 1);
		Qureg qs1 = rand_qureg_sparse(nqubit, half_fill(nqubit), 2, false);
		Qureg qs2 = rand_qureg_sparse(nqubit, half_fill(nqubit) / 2 + 1, 1, true);
		Qureg QQs[] = { qd, qs1, qs2 };

		qubase t;
		for (Qureg& q : QQs)
		for (int tar : Range<>(nqubit))
		{
			Matrix2cf mat = rand_cxmat(2, 2);
			VectorXcf oldAmp = VectorXcf(q);
			generic_gate(q, mat, tar);
			t = q.to_bit(tar);
			VectorXcf newAmp = VectorXcf(q);
			for (qubase base : Range<>(1 << nqubit))
			{
				if (base & t) continue; // symmetry
				oldBitAmp << oldAmp(base), oldAmp(base ^ t);
				newBitAmp << newAmp(base), newAmp(base ^ t);
				ASSERT_MAT(mat * oldBitAmp, newBitAmp);
			}
		}
	}
}

TEST(Qugate, GenericGate2)
{
	Vector4cf oldBitAmp, newBitAmp;
	for (int nqubit : QubitRange(3))
	{
		Qureg qd = rand_qureg_dense(nqubit, 1);
		Qureg qs1 = rand_qureg_sparse(nqubit, half_fill(nqubit), 2, false);
		Qureg qs2 = rand_qureg_sparse(nqubit, half_fill(nqubit) / 2 + 1, 1, true);
		Qureg QQs[] = { qd, qs1, qs2 };

		qubase t1, t2;
		for (Qureg& q : QQs)
		for (int tar1 : Range<>(nqubit-1))
		for (int tar2 : Range<>(tar1+1, nqubit))
		{
			Matrix4cf mat = rand_cxmat(4, 4);
			VectorXcf oldAmp = VectorXcf(q);
			generic_gate(q, mat, tar1, tar2);
			t1 = q.to_bit(tar1); t2 = q.to_bit(tar2);
			VectorXcf newAmp = VectorXcf(q);
			for (qubase base : Range<>(1 << nqubit))
			{
				if ((base & t1) || (base & t2)) continue; // symmetry
				oldBitAmp << 
					oldAmp(base), oldAmp(base ^ t1), oldAmp(base ^ t2), oldAmp(base ^ (t1 | t2));
				newBitAmp << 
					newAmp(base), newAmp(base ^ t1), newAmp(base ^ t2), newAmp(base ^ (t1 | t2));
				ASSERT_MAT(mat * oldBitAmp, newBitAmp);
			}
		}
	}
}