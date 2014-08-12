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
			for (Qureg& qq : QQs)
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
		Qureg QQs[] = { move(qd), move(qs1), move(qs2) };

		qubase t;
		for (Qureg& q : QQs)
		for (int tar : Range<>(nqubit))
		{
			Matrix2cf mat = rand_cxmat(2, 2);
			VectorXcf oldAmp = VectorXcf(q);
			generic_gate(q, mat, tar);
			t = q.to_qubase(tar);
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
	for (int nqubit : QubitRange(3))
	{
		Qureg qd = rand_qureg_dense(nqubit, 1);
		Qureg qs1 = rand_qureg_sparse(nqubit, half_fill(nqubit), 1, false);
		Qureg qs2 = rand_qureg_sparse(nqubit, half_fill(nqubit) / 2 + 1, 1, true);
		Qureg QQs[] = { move(qd), move(qs1), move(qs2) };

		for (Qureg& q : QQs)
		for (int tar1 : Range<>(nqubit-1))
		for (int tar2 : Range<>(tar1+1, nqubit))
		{
			Matrix2cf mat1 = rand_cxmat(2, 2, .5);
			Matrix2cf mat2 = rand_cxmat(2, 2, .5);

			// Apply mat1 first, then mat2
			Qureg qc1 = q.clone();
			generic_gate(qc1, mat1, tar1);
			generic_gate(qc1, mat2, tar2);
			VectorXcf newAmp1 = VectorXcf(qc1);

			// Apply mat2 first, then mat1
			Qureg qc2 = q.clone();
			generic_gate(qc2, mat2, tar2);
			generic_gate(qc2, mat1, tar1);
			VectorXcf newAmp2 = VectorXcf(qc2);

			// Apply mat1 and mat2 at the same time by Matrix4cf
			Qureg qc3 = q.clone();
			generic_gate(qc3, (Matrix4cf) kronecker_mat(mat1, mat2), tar1, tar2);
			VectorXcf newAmp3 = VectorXcf(qc3);

			// Apply mat2 and mat1 at the same time by Matrix4cf
			generic_gate(q, (Matrix4cf) kronecker_mat(mat2, mat1), tar2, tar1);
			VectorXcf newAmp4 = VectorXcf(q);

			ASSERT_MAT(newAmp1, newAmp2, "mat1 then mat2 VS mat2 then mat1");
			ASSERT_MAT(newAmp1, newAmp3, "mat1 then mat2 VS mat1+mat2");
			ASSERT_MAT(newAmp3, newAmp4, "mat1+mat2 VS mat2+mat1");
		}
	}
}

TEST(Qugate, GenericGateN)
{
	const int NQUBIT = 5;
	vector<int> tars(NQUBIT);
	vector<Matrix2cf> mats(NQUBIT);
	for (int nqubit : QubitRange(NQUBIT))
	{
		Qureg qd = rand_qureg_dense(nqubit, 1);
		Qureg qs1 = rand_qureg_sparse(nqubit, half_fill(nqubit), 1, false);
		Qureg qs2 = rand_qureg_sparse(nqubit, half_fill(nqubit) / 2 + 1, 1, true);
		Qureg QQs[] = { move(qd), move(qs1), move(qs2) };

		for (Qureg& q : QQs)
		for (int trial : Range<>(20))
		{
			rand_shuffle(rand_unique(tars, NQUBIT, nqubit));
			Qureg qc1 = q.clone();
			MatrixXcf kroneckeredMat;
			bool first = true;
			for (int i = 0; i < NQUBIT; ++i)
			{
				Matrix2cf mat = rand_cxmat(2, 2, .5);
				if (first)
				{
					kroneckeredMat = mat;
					first = false;
				}
				else
					kroneckeredMat = kroneckeredMat & mat;
				// Apply single-qubit gate
				generic_gate(qc1, mat, tars[i]);
			}
			VectorXcf newAmp1 = VectorXcf(qc1);

			generic_gate(q, kroneckeredMat, tars);
			VectorXcf newAmp2 = VectorXcf(q);

			ASSERT_MAT(newAmp1, newAmp2);
		}
	}
}

TEST(Qugate, PauliXYZ)
{
	std::function<void(Qureg&, int)> 
		pauliFuncs[] = { pauli_X, pauli_Y, pauli_Z };

	std::function<Matrix2cf()> 
		verifyMats[] = { pauli_X_mat, pauli_Y_mat, pauli_Z_mat };

	Vector2cf oldBitAmp, newBitAmp;

	for (int fi = 0; fi < 3; ++fi)
	{
		auto pauliFunc = pauliFuncs[fi];
		auto verifyMat = verifyMats[fi];

		for (int nqubit : QubitRange(2))
		{
			Qureg qd = rand_qureg_dense(nqubit, 1);
			Qureg qs1 = rand_qureg_sparse(nqubit, half_fill(nqubit), 2, false);
			Qureg qs2 = rand_qureg_sparse(nqubit, half_fill(nqubit) / 2 + 1, 1, true);
			Qureg QQs[] = { move(qd), move(qs1), move(qs2) };

			qubase t;
			for (Qureg& q : QQs)
			for (int tar : Range<>(nqubit))
			{
				VectorXcf oldAmp = VectorXcf(q);
				pauliFunc(q, tar);
				t = q.to_qubase(tar);
				VectorXcf newAmp = VectorXcf(q);
				for (qubase base : Range<>(1 << nqubit))
				{
					if (base & t) continue; // symmetry
					oldBitAmp << oldAmp(base), oldAmp(base ^ t);
					newBitAmp << newAmp(base), newAmp(base ^ t);
					ASSERT_MAT(verifyMat() * oldBitAmp, newBitAmp, 
							   fi == 0 ? "pauliX" : fi == 1 ? "pauliY" : "pauliZ");
				}
			}
		}
	}
}