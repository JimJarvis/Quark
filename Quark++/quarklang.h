/**********************************************/
/*********** Specific for quarklang compiler ***********/
/**********************************************/
#ifndef quarklang_h__
#define quarklang_h__
#include "utils.h"
#include "qureg.h"

template <typename T>
vector<T> concat_vector(vector<T> vec1, vector<T> vec2)
{
	vector<T> ans = vec1;
	ans.insert(ans.end(), vec2.begin(), vec2.end());
	return ans;
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

template<typename T>
int len(vector<T>& vec) { return vec.size(); }
template<typename T>
int len(vector<T>&& vec) { return vec.size(); }

int qsize(Qureg& q) { return q.nqubit; }
int qsize(Qureg&& q) { return q.nqubit; }

template<typename T>
int rowdim(Matrix<T, Dynamic, Dynamic>& mat) { return mat.rows(); }
template<typename T>
int rowdim(Matrix<T, Dynamic, Dynamic>&& mat) { return mat.rows(); }

template<typename T>
int coldim(Matrix<T, Dynamic, Dynamic>& mat) { return mat.cols(); }
template<typename T>
int coldim(Matrix<T, Dynamic, Dynamic>&& mat) { return mat.cols(); }

#endif // quarklang_h__

