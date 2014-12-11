#include "qureg.h"
#include "qumat.h"
#include "qugate.h"
#include "quarklang.h"
using namespace Qumat;
using namespace Qugate;

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

int main(int argc, char **argv)
{
	pr(matrix_literal<float>(3, vector<float>{6, 4, 5, 10, 11, 12}));
	pr(pow(complex<float>(4), 3));

	pr(concat_vector(vector<float>{3, 2, 1}, vector<float>{5, 4, 10}));
	pr( (vector<vector<float>>{{1, 5, 6, -3}, {10, 3, 2}, { 6 }}) );
}

