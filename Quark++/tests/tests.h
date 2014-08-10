#ifndef tests_h__
#define tests_h__

#include "gtest/gtest.h"
#include "../utils.h"
#include "../qureg.h"
#include "../qugate.h"
#include "../qumat.h"
using namespace Qumat;
using namespace Qugate;
using namespace Eigen;

///////************** Test conventions **************///////
#define QubitRange Range<>(1, 8)
#define QubaseRange [](int nqubit) { return Range<qubase>(1 << nqubit); }

/*
 *	Error stream must be terminated by "Eend"
 * ErrSS << 23 << "dudulu" << Eend
 */
inline void ASSERT_MAT(const MatrixXcf& m1, const MatrixXcf& m2, const string& errstr, float tol = TOL)
{
	size_t r, c;
	ASSERT_EQ(r = m1.rows(), m2.rows()) << "Row dims should agree";
	ASSERT_EQ(c = m1.cols(), m2.cols()) << "Col dims should agree";

	auto genError = [&](size_t i, size_t j)
	{
		ostringstream oss;
		oss << "Disagree at [" << i << ", " << j << "]: "
			<< m1(i, j) << " vs " << m2(i, j) << endl;
		return oss.str() + errstr;
	};

	for (size_t i = 0; i < r ; ++i)
		for (size_t j = 0; j < c; ++j)
		{
			CX a1 = m1(i, j);
			CX a2 = m2(i, j);
			ASSERT_NEAR(a1.real(), a2.real(), tol) << genError(i, j);
			ASSERT_NEAR(a1.imag(), a2.imag(), tol) << genError(i, j);
		}
}

#endif // tests_h__
