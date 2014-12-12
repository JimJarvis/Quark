/**********************************************/
/*********** Specific for quarklang compiler ***********/
/**********************************************/
#ifndef quarklang_h__
#define quarklang_h__
#include "utils.h"
#include "qureg.h"
#include "qugate.h"

// [1, 2] & [3,4,5]
#define GEN_CONCAT_VECTOR(type1, type2) \
template <typename T> \
vector<T> concat_vector(type1 vec1, type2 vec2) \
{ \
	vector<T> ans = vec1; \
	ans.insert(ans.end(), vec2.begin(), vec2.end()); \
	return ans; \
}

GEN_CONCAT_VECTOR(vector<T>&, vector<T>&)
GEN_CONCAT_VECTOR(vector<T>&&, vector<T>&)
GEN_CONCAT_VECTOR(vector<T>&, vector<T>&&)
GEN_CONCAT_VECTOR(vector<T>&&, vector<T>&&)

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
bool membership_in(T elem, vector<T>& vec)
{
	return std::find(vec.begin(), vec.end(), elem) != vec.end();
}
template<typename T>
bool membership_in(T elem, vector<T>&& vec)
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
void generic_1gate(Q, const Matrix2cf mat, int tar)
{
	Qugate::generic_gate(q, mat, tar);
}

void generic_2gate(Q, const Matrix4cf mat, int tar1, int tar2)
{
	Qugate::generic_gate(q, mat, tar1, tar2);
}
/*
*	Works with arbitrary number of target qubits
*/
void generic_ngate(Q, const MatrixXcf mat, vector<int> tars)
{
	Qugate::generic_gate(q, mat, tars);
}

#endif // quarklang_h__