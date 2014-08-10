#ifndef tests_h__
#define tests_h__

#include "gtest/gtest.h"
#include "../utils.h"
#include "../qureg.h"
#include "../qugate.h"
#include "../quop.h"
using namespace Quop;
using namespace Qugate;
using namespace Eigen;

// Dummy struct for ErrorStream
struct Endl {};
struct ErrorStream
{
	ostringstream erross;
	bool clear = false;
	template<typename T>
	ErrorStream& operator<<(const T& obj)
	{
		if (clear)
		{
			erross.clear();
			clear = false;
		}
		erross << obj;
		return *this;
	}
	
	template<>
	ErrorStream& operator<< <Endl>(const Endl& obj)
	{
		clear = true;
		return *this;
	}

	string str() { return erross.str(); }
};

// Global
static ErrorStream ErrSS;
static Endl Eend;

/*
 *	Debugging matrix comparison
 */
void ASSERT_MAT(const MatrixXcf& m1, const MatrixXcf& m2, ErrorStream& erross, float tol = TOL)
{
	size_t r, c;
	ASSERT_EQ(r = m1.rows(), m2.rows()) << "Row dims should agree";
	ASSERT_EQ(c = m1.cols(), m2.cols()) << "Col dims should agree";

	auto genError = [&](size_t i, size_t j)
	{
		ostringstream oss;
		oss << "Disagree at [" << i << ", " << j << "]: "
			<< m1(i, j) << " vs " << m2(i, j) << endl;
		return oss.str() + erross.str();
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
