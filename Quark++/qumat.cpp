/**********************************************
* Quantum operations  *
**********************************************/
#include "qumat.h"
#include "qureg.h"

using namespace Qumat;

Qureg Qumat::kronecker(Q1, Q2, bool resultDense)
{
	int new_nqubit = q1.nqubit + q2.nqubit;
	size_t size1 = q1.size();
	size_t size2 = q2.size();

	Qureg qans = resultDense ?
		Qureg::create<true>(new_nqubit) :
		Qureg::create<false>(new_nqubit, size1 * size2);
	if (resultDense)  qans.set_base_d(qubase(0), CX(0));

	// If both of them are dense, then the resultant is also dense
	int ri = 0;
	int nqubit2 = q2.nqubit;
	for (int i1 = 0; i1 < size1; ++i1)
		for (int i2 = 0; i2 < size2; ++i2)
		{
			qubase newBase = (q1.get_base_d_s(i1) << nqubit2) | q2.get_base_d_s(i2);
			CX newAmp = q1.amp[i1] * q2.amp[i2];
			if (resultDense)
				qans.amp[newBase] = newAmp;
			else
				qans.add_base(newBase, newAmp);
		}

	return qans;
}

/**********************************************/
/*********** Eigen  ***********/
/**********************************************/
MatrixXcf Qumat::hadamard_mat(int nqubit)
{
	size_t size = 1 << nqubit;
	MatrixXcf mat(size, size);
	float coef = 1.0 / sqrt(size);
	for (size_t i = 0; i < size ; ++i)
		for (size_t j = 0; j < size; ++j)
			// binary dot product
			mat(i, j) = CX((bitwise_dot(i, j) ? -1 : 1) * coef);
	return mat;
}

Matrix2cf Qumat::hadamard_mat()
{
	static CX _sqrt2 = CX(1 / sqrt(2));
	static Matrix2cf HadamardMat;
	INIT_ONCE(
		HadamardMat <<
			_sqrt2, _sqrt2,
			_sqrt2, -_sqrt2;
	)
	return HadamardMat;
}

MatrixXcf Qumat::kronecker_mat(const MatrixXcf& A, const MatrixXcf& B)
{
	size_t ar = A.rows();
	size_t ac = A.cols();
	size_t br = B.rows();
	size_t bc = B.cols();
	MatrixXcf mat(ar * br, ac * bc);
	for (int i = 0; i < ar ; ++i)
		for (int j = 0; j < ac; ++j)
			mat.block(i * br, j * bc, br, bc) = B * A(i, j);
	return mat;
}

Matrix4cf Qumat::cnot_mat()
{
	static Matrix4cf CnotMat;
	INIT_ONCE(
		CnotMat <<
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 0, 1,
			0, 0, 1, 0; 
	)
	return CnotMat;
}

Matrix<CX, 8, 8> Qumat::toffoli_mat()
{
	static Matrix<CX, 8, 8> ToffoliMat;
	static MatrixXcf zero6_2 = MatrixXcf::Zero(6, 2);
	INIT_ONCE(
	ToffoliMat << 
		MatrixXcf::Identity(6, 6), zero6_2, 
		zero6_2.transpose(), MatrixXcf::Identity(2, 2).colwise().reverse();
	)
	return ToffoliMat;
}

MatrixXcf Qumat::toffoli_mat(int nctrl)
{
	static MatrixXcf ToffoliMat;
	INIT_ONCE(
	size_t size = 1 << (nctrl + 1);
	ToffoliMat = MatrixXcf::Identity(size, size);
	ToffoliMat.block(size - 2, size - 2, 2, 2) 
		= MatrixXcf::Identity(2, 2).colwise().reverse();
	)
	return ToffoliMat;
}

MatrixXcf Qumat::generic_control_mat(int nctrl, const Matrix2cf& mat)
{
	size_t size = 1 << (nctrl + 1);
	MatrixXcf GenericCtrlMat = MatrixXcf::Identity(size, size);
	GenericCtrlMat.block(size - 2, size - 2, 2, 2) = mat;
	return GenericCtrlMat;
}

Matrix2cf Qumat::pauli_X_mat()
{
	static Matrix2cf pauliXMat;
	INIT_ONCE(
		pauliXMat = Matrix2cf::Identity(2, 2).colwise().reverse();
	)
	return pauliXMat;
}

Matrix2cf Qumat::pauli_Y_mat()
{
	static Matrix2cf pauliYMat;
	INIT_ONCE(
		pauliYMat <<
			0, CX(0, -1),
			CX(0, 1), 0
	)
	return pauliYMat;
}

Matrix2cf Qumat::pauli_Z_mat()
{
	static Matrix2cf pauliZMat;
	INIT_ONCE(
		pauliZMat = Matrix2cf::Identity(2, 2);
		pauliZMat(1, 1) = -1;
	)
	return pauliZMat;
}

Matrix2cf Qumat::rot_X_mat(float theta)
{
	theta *= 0.5;
	float c = cos(theta);
	float s = sin(theta);
	static Matrix2cf RotXMat;
	RotXMat <<
		c, CX(0, -s),
		CX(0, -s), c;
	return RotXMat;
}

Matrix2cf Qumat::rot_Y_mat(float theta)
{
	theta *= 0.5;
	float c = cos(theta);
	float s = sin(theta);
	static Matrix2cf RotYMat;
	RotYMat <<
		c, -s,
		s, c;
	return RotYMat;
}

Matrix2cf Qumat::rot_Z_mat(float theta)
{
	CX x = expi(theta / 2);
	static Matrix2cf RotZMat;
	RotZMat <<
		conj(x), 0,
		0, x;
	return RotZMat;
}

Matrix2cf Qumat::phase_scale_mat(float theta)
{
	CX x = expi(theta);
	static Matrix2cf PhaseScaleMat;
	PhaseScaleMat <<
		x, 0,
		0, x;
	return PhaseScaleMat;
}

Matrix2cf Qumat::phase_shift_mat(float theta)
{
	CX x = expi(theta);
	static Matrix2cf PhaseScaleMat;
	PhaseScaleMat <<
		1, 0,
		0, x;
	return PhaseScaleMat;
}

Matrix4cf Qumat::control_phase_shift_mat(float theta)
{
	static Matrix4cf CtrlPhaseShiftMat;
	CtrlPhaseShiftMat = Matrix4cf::Identity(4, 4);
	CtrlPhaseShiftMat(3, 3) = expi(theta);
	return CtrlPhaseShiftMat;
}

Matrix4cf Qumat::swap_mat()
{
	static Matrix4cf SwapMat;
	INIT_ONCE(
		SwapMat <<
			1, 0, 0, 0,
			0, 0, 1, 0,
			0, 1, 0, 0,
			0, 0, 0, 1;
	)
	return SwapMat;
}

Matrix<CX, 8, 8> Qumat::cswap_mat()
{
	static Matrix<CX, 8, 8> CswapMat;
	INIT_ONCE(
		CswapMat = Matrix<CX, 8, 8>::Identity(8, 8);
		CswapMat.block(5, 5, 2, 2) = Matrix2cf::Identity(2, 2).colwise().reverse();
	)
	return CswapMat;
}
