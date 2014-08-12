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
				if (base & t) base ^= t;
				oldBitAmp << oldAmp(base), oldAmp(base ^ t);
				newBitAmp << newAmp(base), newAmp(base ^ t);
				ASSERT_MAT(mat * oldBitAmp, newBitAmp);
			}
		}
	}
}

// Helper
// 'process': whether ctrl bit is 1 or not
void test_generic_ctrl(
	bool isCtrlOn, VectorXcf& oldAmp, VectorXcf& newAmp, qubase& base, qubase& t, Matrix2cf& mat)
{
	Vector2cf oldBitAmp, newBitAmp;
	if (base & t)
	{
		oldBitAmp << oldAmp(base ^ t), oldAmp(base);
		newBitAmp << newAmp(base ^ t), newAmp(base);
	}
	else
	{
		oldBitAmp << oldAmp(base), oldAmp(base ^ t);
		newBitAmp << newAmp(base), newAmp(base ^ t);
	}

	if (isCtrlOn)
		ASSERT_MAT(mat * oldBitAmp, newBitAmp, "Controlled");
	else
		ASSERT_MAT(oldBitAmp, newBitAmp, "Uncontrolled");
}

TEST(Qugate, SimpleCnot)
{
	// pre-alloc for 2 rand bits
	vector<int> randBitVec(2);
	static Matrix2cf idrev2 = Matrix2cf::Identity(2, 2).colwise().reverse();

	for (int nqubit : QubitRange(2))
	{
		Qureg qd = rand_qureg_dense(nqubit, 1);
		Qureg qs1 = rand_qureg_sparse(nqubit, half_fill(nqubit), 2, false);
		Qureg qs2 = rand_qureg_sparse(nqubit, half_fill(nqubit)/2+1, 1, true);

		Qureg QQs[] = { qd, qs1, qs2 };

		VectorXcf oldAmp, newAmp;
		qubase c, t; // ctrl and target
		for (Qureg& q : QQs)
		for (int trial : Range<>(20))
		{
			oldAmp = VectorXcf(q);
			// generate two random bits
			rand_shuffle(rand_unique(randBitVec, 2, nqubit));
			c = randBitVec[0]; t = randBitVec[1];

			cnot(q, c, t);
			newAmp = VectorXcf(q);

			c = q.to_bit(c); t = q.to_bit(t);
			for (qubase base : Range<>(1 << nqubit))
				test_generic_ctrl(base & c, oldAmp, newAmp, base, t, idrev2);
		}
	}
}

TEST(Qugate, GenericCnot)
{
	// pre-alloc for 2 rand bits
	vector<int> randBitVec(2);
	for (int nqubit : QubitRange(2))
	{
		Qureg qd = rand_qureg_dense(nqubit, 1);
		Qureg qs1 = rand_qureg_sparse(nqubit, half_fill(nqubit), 2, false);
		Qureg qs2 = rand_qureg_sparse(nqubit, half_fill(nqubit) / 2 + 1, 1, true);

		Qureg QQs[] = { qd, qs1, qs2 };

		VectorXcf oldAmp, newAmp;
		qubase c, t; // ctrl and target
		for (Qureg& q : QQs)
		for (int trial : Range<>(20))
		{
			Matrix2cf mat = rand_cxmat(2, 2);
			oldAmp = VectorXcf(q);
			// generate two random bits
			rand_shuffle(rand_unique(randBitVec, 2, nqubit));
			c = randBitVec[0]; t = randBitVec[1];

			generic_control(q, mat, c, t);
			newAmp = VectorXcf(q);

			c = q.to_bit(c); t = q.to_bit(t);
			for (qubase base : Range<>(1 << nqubit))
				test_generic_ctrl(base & c, oldAmp, newAmp, base, t, mat);
		}
	}
}

TEST(Qugate, SimpleToffoli)
{
	// pre-alloc for 3 rand bits
	vector<int> randBitVec(3);
	static Matrix2cf idrev2 = Matrix2cf::Identity(2, 2).colwise().reverse();

	for (int nqubit : QubitRange(3))
	{
		Qureg qd = rand_qureg_dense(nqubit, 1);
		Qureg qs1 = rand_qureg_sparse(nqubit, half_fill(nqubit), 2, false);
		Qureg qs2 = rand_qureg_sparse(nqubit, half_fill(nqubit) / 2 + 1, 1, true);

		Qureg QQs[] = { qd, qs1, qs2 };

		VectorXcf oldAmp, newAmp;
		qubase c1, c2, t; // ctrl and target
		for (Qureg& q : QQs)
		for (int trial : Range<>(20))
		{
			oldAmp = VectorXcf(q);
			// generate two random bits
			rand_shuffle(rand_unique(randBitVec, 3, nqubit));
			c1 = randBitVec[0]; c2 = randBitVec[1]; t = randBitVec[2];

			// Apply Toffoli
			toffoli(q, c1, c2, t);

			newAmp = VectorXcf(q);

			c1 = q.to_bit(c1); c2 = q.to_bit(c2); t = q.to_bit(t);
			for (qubase base : Range<>(1 << nqubit))
				test_generic_ctrl((base & c1) && (base & c2), 
						oldAmp, newAmp, base, t, idrev2);
		}
	}
}

TEST(Qugate, GenericToffoli)
{
	// pre-alloc for 3 rand bits
	vector<int> randBitVec(3);

	for (int nqubit : QubitRange(3))
	{
		Qureg qd = rand_qureg_dense(nqubit, 1);
		Qureg qs1 = rand_qureg_sparse(nqubit, half_fill(nqubit), 2, false);
		Qureg qs2 = rand_qureg_sparse(nqubit, half_fill(nqubit) / 2 + 1, 1, true);

		Qureg QQs[] = { qd, qs1, qs2 };

		VectorXcf oldAmp, newAmp;
		qubase c1, c2, t; // ctrl and target
		for (Qureg& q : QQs)
		for (int trial : Range<>(20))
		{
			Matrix2cf mat = rand_cxmat(2, 2);
			oldAmp = VectorXcf(q);
			// generate two random bits
			rand_shuffle(rand_unique(randBitVec, 3, nqubit));
			c1 = randBitVec[0]; c2 = randBitVec[1]; t = randBitVec[2];

			// Apply Toffoli
			generic_toffoli(q, mat, c1, c2, t);

			newAmp = VectorXcf(q);

			c1 = q.to_bit(c1); c2 = q.to_bit(c2); t = q.to_bit(t);
			for (qubase base : Range<>(1 << nqubit))
				test_generic_ctrl((base & c1) && (base & c2),
				oldAmp, newAmp, base, t, mat);
		}
	}
}

TEST(Qugate, SimpleNcnot)
{
	// pre-alloc for 3 rand bits
	const int NCNOT = 6;
	vector<int> randBitVec(NCNOT);
	static Matrix2cf idrev2 = Matrix2cf::Identity(2, 2).colwise().reverse();

	for (int nqubit : Range<>(NCNOT, 9))
	{
		Qureg qd = rand_qureg_dense(nqubit, 1);
		Qureg qs1 = rand_qureg_sparse(nqubit, half_fill(nqubit), 2, false);
		Qureg qs2 = rand_qureg_sparse(nqubit, half_fill(nqubit) / 2 + 1, 1, true);

		Qureg QQs[] = { qd, qs1, qs2 };

		VectorXcf oldAmp, newAmp;
		qubase t; // ctrl and target
		for (Qureg& q : QQs)
		for (int trial : Range<>(20))
		{
			oldAmp = VectorXcf(q);
			// generate two random bits
			rand_shuffle(rand_unique(randBitVec, NCNOT, nqubit));
			t = randBitVec[randBitVec.size() - 1];

			ncnot(q, 
				  vector<int>(randBitVec.begin(), randBitVec.begin() + randBitVec.size()-1), 
				  t);

			newAmp = VectorXcf(q);

			vector<qubase> ctrlBasis;
			for (int i = 0; i < randBitVec.size() - 1; ++i)
				ctrlBasis.push_back(q.to_bit(randBitVec[i]));
			t = q.to_bit(t);

			bool isCtrlOn;
			for (qubase base : Range<>(1 << nqubit))
			{
				isCtrlOn = true;
				for (qubase& ctrl : ctrlBasis)
					if (!(base & ctrl))
					{
						isCtrlOn = false;
						break;
					}
				test_generic_ctrl(isCtrlOn, oldAmp, newAmp, base, t, idrev2);
			}
		}
	}
}

TEST(Qugate, GenericNcnot)
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

		VectorXcf oldAmp, newAmp;
		qubase t; // ctrl and target
		for (Qureg& q : QQs)
		for (int trial : Range<>(20))
		{
			Matrix2cf mat = rand_cxmat(2, 2);
			oldAmp = VectorXcf(q);
			// generate two random bits
			rand_shuffle(rand_unique(randBitVec, NCNOT, nqubit));
			t = randBitVec[randBitVec.size() - 1];

			generic_ncontrol(q, mat,
				  vector<int>(randBitVec.begin(), randBitVec.begin() + randBitVec.size() - 1),
				  t);

			newAmp = VectorXcf(q);

			vector<qubase> ctrlBasis;
			for (int i = 0; i < randBitVec.size() - 1; ++i)
				ctrlBasis.push_back(q.to_bit(randBitVec[i]));
			t = q.to_bit(t);

			bool isCtrlOn;
			for (qubase base : Range<>(1 << nqubit))
			{
				isCtrlOn = true;
				for (qubase& ctrl : ctrlBasis)
					if (!(base & ctrl))
					{
						isCtrlOn = false;
						break;
					}
				test_generic_ctrl(isCtrlOn, oldAmp, newAmp, base, t, mat);
			}
		}
	}
}