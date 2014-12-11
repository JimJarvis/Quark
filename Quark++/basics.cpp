#include "qureg.h"
#include "qumat.h"
#include "qugate.h"
#include "quarklang.h"

using namespace Qumat;
using namespace Qugate;

int main(int argc, char **argv)
{
	pr(matrix_literal<float>(3, vector<float>{6, 4, 5, 10, 11, 12}));
	pr(pow(complex<float>(4), 3));

	pr(concat_vector(vector<float>{3, 2, 1}, vector<float>{5, 4, 10}));
	pr( (vector<vector<float>>{{1, 5, 6, -3}, {10, 3, 2}, { 6 }}) );
	auto v = (vector<vector<float>>{{1, 5, 6, -3}, {10, 3, 2}, { 6 }});
	pr(len(v));
	pr(len(vector<vector<float>>{{1, 5, 6, -3}, {10, 3, 2}, { 6 }}));
	pr( (vector<vector<float>>{vector<float>{1, 5, 6, -3}, vector<float>{10, 3, 2}, vector<float>{ 6 }}) );
}

