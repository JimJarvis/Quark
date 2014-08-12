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

TEST(Qugate, Generic_Gate)
{
	static Matrix2cf mat;
	mat << 2, CX(-1, .5), CX(.3, -.1), CX(2, -1);

	Vector2cf oldAmp, newAmp;
	for (int nqubit : QubitRange(2))
	{
		Qureg qd = rand_qureg_dense(nqubit, 1);
		Qureg qs1 = rand_qureg_sparse(nqubit, half_fill(nqubit), 2, false);
		Qureg qs2 = rand_qureg_sparse(nqubit, half_fill(nqubit) / 2 + 1, 1, true);

		Qureg QQs[] = { qd, qs1, qs2 };

		VectorXcf vec, vecNew;
		qubase t;
		for (Qureg& q : QQs)
		for (int tar : Range<>(nqubit))
		{
			vec = VectorXcf(q);
			generic_gate(q, mat, tar);
			t = q.to_bit(tar);
			vecNew = VectorXcf(q);
			for (qubase base : Range<>(1 << nqubit))
			{
				if (base & t) base ^= t;
				oldAmp << vec(base), vec(base ^ t);
				newAmp << vecNew(base), vecNew(base ^ t);
				ASSERT_MAT(mat * oldAmp, newAmp);
			}
		}
	}
}

TEST(Qugate, Cnot)
{
	// pre-alloc for 2 rand bits
	vector<int> randBitVec(2);
	Vector2cf oldAmp, newAmp;
	static Matrix2cf idrev2 = Matrix2cf::Identity(2, 2).colwise().reverse();
	static Matrix2cf mat;
	mat << 2, CX(-1, .5), CX(.3, -.1), CX(2, -1);

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
			rand_shuffle(rand_unique(randBitVec, 2, nqubit));
			c = randBitVec[0];
			t = randBitVec[1];

			// Apply CNOT
			//cnot(q, c, t);
			generic_control(q, mat, c, t);
			vecNew = VectorXcf(q);

			c = q.to_bit(c);
			t = q.to_bit(t);

			for (qubase base : Range<>(1 << nqubit))
			{
				if (base & t)
				{
					oldAmp << vec(base ^ t), vec(base);
					newAmp << vecNew(base ^ t), vecNew(base);
				}
				else
				{
					oldAmp << vec(base), vec(base ^ t);
					newAmp << vecNew(base), vecNew(base ^ t);
				}

				if (base & c)
					ASSERT_MAT(mat * oldAmp, newAmp);
				else
					ASSERT_MAT(oldAmp, newAmp, "Uncontrolled");
			}
		}
	}
}

TEST(Qugate, Toffoli)
{
	// pre-alloc for 3 rand bits
	vector<int> randBitVec(3);
	for (int nqubit : QubitRange(3))
	{
		Qureg qd = rand_qureg_dense(nqubit, 1);
		Qureg qs1 = rand_qureg_sparse(nqubit, half_fill(nqubit), 2, false);
		Qureg qs2 = rand_qureg_sparse(nqubit, half_fill(nqubit) / 2 + 1, 1, true);

		Qureg QQs[] = { qd, qs1, qs2 };

		VectorXcf vec, vecNew;
		qubase c1, c2, t; // ctrl and target
		for (Qureg& q : QQs)
		for (int trial : Range<>(20))
		{
			vec = VectorXcf(q);
			// generate two random bits
			rand_shuffle(rand_unique(randBitVec, 3, nqubit));
			c1 = randBitVec[0];
			c2 = randBitVec[1];
			t = randBitVec[2];

			// Apply Toffoli
			toffoli(q, c1, c2, t);

			vecNew = VectorXcf(q);

			c1 = q.to_bit(c1);
			c2 = q.to_bit(c2);
			t = q.to_bit(t);

			for (qubase base : Range<>(1 << nqubit))
			{
				ASSERT_CX_EQ(vec(base),
							 vecNew((base & c1) && (base & c2)
									 ? base ^ t : base),
							 "base is " << bits2str(base));
			}
		}
	}
}

TEST(Qugate, Ncnot)
{
	// pre-alloc for 3 rand bits
	const int NCNOT = 6;
	vector<int> randBitVec(NCNOT);
	for (int nqubit : Range<>(NCNOT, 9))
	{
		Qureg qd = rand_qureg_dense(nqubit, 1);
		Qureg qs1 = rand_qureg_sparse(nqubit, half_fill(nqubit), 2, false);
		Qureg qs2 = rand_qureg_sparse(nqubit, half_fill(nqubit) / 2 + 1, 1, true);

		Qureg QQs[] = { qd, qs1, qs2 };

		VectorXcf vec, vecNew;
		qubase t; // ctrl and target
		for (Qureg& q : QQs)
		for (int trial : Range<>(20))
		{
			vec = VectorXcf(q);
			// generate two random bits
			rand_shuffle(rand_unique(randBitVec, NCNOT, nqubit));
			t = randBitVec[randBitVec.size() - 1];

			// Apply generalized cnot
			ncnot(q, 
				  vector<int>(randBitVec.begin(), randBitVec.begin() + randBitVec.size()-1), 
				  t);

			vecNew = VectorXcf(q);

			vector<qubase> ctrlBasis;
			for (int i = 0; i < randBitVec.size() - 1; ++i)
				ctrlBasis.push_back(q.to_bit(randBitVec[i]));
			t = q.to_bit(t);

			bool process;
			for (qubase base : Range<>(1 << nqubit))
			{
				process = true;
				for (qubase& ctrl : ctrlBasis)
					if (!(base & ctrl))
					{
						process = false;
						break;
					}
				ASSERT_CX_EQ(vec(base),
							 vecNew(process ? base ^ t : base),
							 "base is " << bits2str(base));
			}
		}
	}
}