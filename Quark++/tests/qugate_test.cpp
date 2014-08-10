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

TEST(Qugate, Cnot)
{
	for (int nqubit : QubitRange(2))
	{
		Qureg qd = rand_qureg_dense(nqubit, 1);
		Qureg qs1 = rand_qureg_sparse(nqubit, half_fill(nqubit), 2, false);
		Qureg qs2 = rand_qureg_sparse(nqubit, half_fill(nqubit)/2+1, 1, true);

		Qureg QQs[] = { qd, qs1, qs2 };

		VectorXcf vec, vecNew;
		qubase c, t; // ctrl and target
		for (Qureg& q : QQs)
		for (int trial : Range<>(20))
		{
			vec = VectorXcf(q);
			// generate two random bits

			c = rand_int(0, nqubit);
			t = rand_int(0, nqubit);
			if (t == c) // avoid collision
				// if last one
				t = c + 1 == nqubit ? c - 1 : c + 1;
			// Apply CNOT
			cnot(q, c, t);
			vecNew = VectorXcf(q);

			c = q.to_bit(c);
			t = q.to_bit(t);
			
			for (qubase base : Range<>(1 << nqubit))
			{
				ASSERT_CX_EQ(vec(base), 
							 vecNew((base & c) ? base ^ t : base), 
							 "base is " << bits2str(base));
			}
		}
	}
}

//TEST(Qugate, Toffoli)
//{
//	for (int nqubit : QubitRange(3))
//	{
//		Qureg qd = rand_qureg_dense(nqubit, 1);
//		Qureg qs1 = rand_qureg_sparse(nqubit, half_fill(nqubit), 2, false);
//		Qureg qs2 = rand_qureg_sparse(nqubit, half_fill(nqubit) / 2 + 1, 1, true);
//
//		Qureg QQs[] = { qd, qs1, qs2 };
//
//		VectorXcf vec, vecNew;
//		qubase c1, c2, t; // ctrl and target
//		for (Qureg& q : QQs)
//		for (int trial : Range<>(20))
//		{
//			vec = VectorXcf(q);
//			// generate two random bits
//			c1 = rand_int(0, nqubit);
//			c2 = rand_int(0, nqubit);
//			t = rand_int(0, nqubit);
//			if (t == c) // avoid collision
//				// if last one
//				t = c + 1 == nqubit ? c - 1 : c + 1;
//			// Apply CNOT
//			cnot(q, c, t);
//			vecNew = VectorXcf(q);
//
//			c = q.to_bit(c);
//			t = q.to_bit(t);
//
//			for (qubase base : Range<>(1 << nqubit))
//			{
//				ASSERT_CX_EQ(vec(base),
//							 vecNew((base & c) ? base ^ t : base),
//							 "base is " << bits2str(base));
//			}
//		}
//	}
//}