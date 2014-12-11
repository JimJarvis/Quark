/**********************************************/
/*********** Specific for quarklang compiler ***********/
/**********************************************/
#ifndef quarklang_h__
#define quarklang_h__
#include "utils.h"
#include "qureg.h"
#include "qugate.h"

// [1, 2] & [3,4,5]
template <typename T>
vector<T> concat_vector(vector<T> vec1, vector<T> vec2)
{
	vector<T> ans = vec1; ans.insert(ans.end(), vec2.begin(), vec2.end()); return ans;
}

template<typename T>
Matrix<T, Dynamic, Dynamic> matrix_literal(int col, vector<T> vec)
{
	int row = vec.size() / col;
	Matrix<T, Dynamic, Dynamic> mlit = Matrix<T, Dynamic, Dynamic>::Zero(row, col);
	for (int i = 0; i < row; ++i)
	for (int j = 0; j < col; ++j)
		mlit(i, j) = vec[col * i + j];
	return mlit;
}

// keyword 'in'
template<typename T>
bool membership_in(T elem, vector<T> vec)
{
	return std::find(vec.begin(), vec.end(), elem) != vec.end();
}

template<typename T>
int len(vector<T>& vec) { return vec.size(); }
template<typename T>
int len(vector<T>&& vec) { return vec.size(); }

int qsize(Qureg& q) { return q.nqubit; }
// int qsize(Qureg&& q) { return q.nqubit; }

template<typename T>
int rowdim(Matrix<T, Dynamic, Dynamic>& mat) { return mat.rows(); }
template<typename T>
int rowdim(Matrix<T, Dynamic, Dynamic>&& mat) { return mat.rows(); }

template<typename T>
int coldim(Matrix<T, Dynamic, Dynamic>& mat) { return mat.cols(); }
template<typename T>
int coldim(Matrix<T, Dynamic, Dynamic>&& mat) { return mat.cols(); }

//**** Fraction getter
int num(Frac f) { return f.num; }
int denom(Frac f) { return f.denom; }

///////************** Qureg specific **************///////
float prefix_prob(Q, int nbit, int64_t prefix)
{
	return q.prefix_prob(nbit, prefix);
}

Qureg qclone(Q)
{
	return q.clone();
}

///////************** Qugate adaptor **************///////
void generic_gate_1(Q, const Matrix2cf mat, int tar)
{
	Qugate::generic_gate(q, mat, tar);
}

void generic_gate_2(Q, const Matrix4cf mat, int tar1, int tar2)
{
	Qugate::generic_gate(q, mat, tar1, tar2);
}
/*
*	Works with arbitrary number of target qubits
*/
void generic_gate_n(Q, const MatrixXcf mat, vector<int> tars)
{
	Qugate::generic_gate(q, mat, tars);
}

///////************** Single-qubit gates **************///////
//void hadamard(Q_, int tar)
//{
//	Qugate::hadamard(q, tar);
//}
//
//void hadamard_top(Q, int topSize);

//void pauli_X(Q, int tar);
//void pauli_Y(Q, int tar);
//void pauli_Z(Q, int tar);
//
//void rot_X(Q, float theta, int tar);
//void rot_Y(Q, float theta, int tar);
//void rot_Z(Q, float theta, int tar);
//
//void phase_scale(Q, float theta, int tar);
//void phase_shift(Q, float theta, int tar);
//
/////////************** Controlled gates **************///////
//void generic_control(Q, const Matrix2cf&, int ctrl, int tar);
//void cnot(Q, int ctrl, int tar);
//
//void generic_toffoli(Q, const Matrix2cf&, int ctrl1, int ctrl2, int tar);
//void toffoli(Q, int ctrl1, int ctrl2, int tar);
//
//void generic_ncontrol(Q, const Matrix2cf&, vector<int>& ctrls, int tar);
//void ncnot(Q, vector<int>& ctrls, int tar);
//
//void control_phase_shift(Q, float theta, int ctrl, int tar);
//
/////////************** Swap gates **************///////
//void swap(Q, int tar1, int tar2);
//
//void cswap(Q, int ctrl, int tar1, int tar2);
//
/////////************** Special gates **************///////
//// tarSize: number of qubits to be operated on
//void qft(Q, int tarStart, int tarSize);
//inline void qft(Q) { qft(q, 0, q.nqubit); }
//
//// diag([2; 0; 0; ...; 0]) - I, invert amplitude unless the state is 0^n 
//void grover_diffuse(Q, int tarStart, int tarSize);
//inline void grover_diffuse(Q) { grover_diffuse(q, 0, q.nqubit); }



#endif // quarklang_h__

