#include "tests.h"

// Helper
// 'process': whether ctrl bit is 1 or not
void test_generic_ctrl(
	bool isCtrlOn, VectorXcf& oldAmp, VectorXcf& newAmp, qubase& base, qubase& t, Matrix2cf& mat)
{
	Vector2cf oldBitAmp, newBitAmp;
	if (base & t) return; // only process base with target 0, because of symmetry
	oldBitAmp << oldAmp(base), oldAmp(base ^ t);
	newBitAmp << newAmp(base), newAmp(base ^ t);

	if (isCtrlOn)
		ASSERT_MAT(mat * oldBitAmp, newBitAmp, "Controlled");
	else
		ASSERT_MAT(oldBitAmp, newBitAmp, "Uncontrolled");
}

TEST(Qugate, SimpleCnot)
{
	// pre-alloc for 2 rand bits
	vector<int> randBitVec(2);

	for (int nqubit : QubitRange(2))
	{
		Qureg qd = rand_qureg_dense(nqubit, 1);
		Qureg qs1 = rand_qureg_sparse(nqubit, half_fill(nqubit), 2, false);
		Qureg qs2 = rand_qureg_sparse(nqubit, half_fill(nqubit) / 2 + 1, 1, true);

		Qureg QQs[] = { move(qd), move(qs1), move(qs2) };

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

			c = q.to_qubase(c); t = q.to_qubase(t);
			for (qubase base : Range<>(1 << nqubit))
				test_generic_ctrl(base & c, oldAmp, newAmp, base, t, pauli_X_mat());
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

		Qureg QQs[] = { move(qd), move(qs1), move(qs2) };

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

			c = q.to_qubase(c); t = q.to_qubase(t);
			for (qubase base : Range<>(1 << nqubit))
				test_generic_ctrl(base & c, oldAmp, newAmp, base, t, mat);
		}
	}
}

TEST(Qugate, SimpleToffoli)
{
	// pre-alloc for 3 rand bits
	vector<int> randBitVec(3);

	for (int nqubit : QubitRange(3))
	{
		Qureg qd = rand_qureg_dense(nqubit, 1);
		Qureg qs1 = rand_qureg_sparse(nqubit, half_fill(nqubit), 2, false);
		Qureg qs2 = rand_qureg_sparse(nqubit, half_fill(nqubit) / 2 + 1, 1, true);

		Qureg QQs[] = { move(qd), move(qs1), move(qs2) };

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

			c1 = q.to_qubase(c1); c2 = q.to_qubase(c2); t = q.to_qubase(t);
			for (qubase base : Range<>(1 << nqubit))
				test_generic_ctrl((base & c1) && (base & c2),
				oldAmp, newAmp, base, t, pauli_X_mat());
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

		Qureg QQs[] = { move(qd), move(qs1), move(qs2) };

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

			c1 = q.to_qubase(c1); c2 = q.to_qubase(c2); t = q.to_qubase(t);
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

	for (int nqubit : Range<>(NCNOT, 9))
	{
		Qureg qd = rand_qureg_dense(nqubit, 1);
		Qureg qs1 = rand_qureg_sparse(nqubit, half_fill(nqubit), 2, false);
		Qureg qs2 = rand_qureg_sparse(nqubit, half_fill(nqubit) / 2 + 1, 1, true);

		Qureg QQs[] = { move(qd), move(qs1), move(qs2) };

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
				  vector<int>(randBitVec.begin(), randBitVec.begin() + randBitVec.size() - 1),
				  t);

			newAmp = VectorXcf(q);

			vector<qubase> ctrlBasis;
			for (int i = 0; i < randBitVec.size() - 1; ++i)
				ctrlBasis.push_back(q.to_qubase(randBitVec[i]));
			t = q.to_qubase(t);

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
				test_generic_ctrl(isCtrlOn, oldAmp, newAmp, base, t, pauli_X_mat());
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

		Qureg QQs[] = { move(qd), move(qs1), move(qs2) };

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
				ctrlBasis.push_back(q.to_qubase(randBitVec[i]));
			t = q.to_qubase(t);

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