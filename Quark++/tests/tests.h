#ifndef tests_h__
#define tests_h__

#include "gtest/gtest.h"
#include "../utils.h"
#include "../qureg.h"
#include "../qugate.h"
#include "../qumat.h"
#include "../algor.h"
#include "../quarklang.h"
using namespace Qumat;
using namespace Qugate;
using namespace Eigen;

///////************** Test conventions **************///////
#define QubitRange(start) Range<>(start, 8)
#define QubaseRange(nqubit) Range<qubase>(1 << nqubit)

// The MSVC has a bug when parsing '__VA_ARGS__'. Workaround:
#define VA_EXPAND(x) x
// always return the fifth argument in place
#define VARARG_INDEX(_0, _1, _2, _3, _4, _5, ...) _5
// how many variadic parameters?
#define VARARG_COUNT(...) VA_EXPAND(VARARG_INDEX(__VA_ARGS__, 5, 4, 3, 2, 1))
#define VARARG_HELPER2(base, count, ...) base##_##count(__VA_ARGS__)
#define VARARG_HELPER(base, count, ...) VARARG_HELPER2(base, count, __VA_ARGS__)
#define VARARG(base, ...) VARARG_HELPER(base, VARARG_COUNT(__VA_ARGS__), __VA_ARGS__)
// Define DEBUG_MSG_1 or _2 or _n to define a debug message printout macro with n args
// Warning: intelliSense might underline this as syntax error. Ignore it and compile. 
#define ASSERT_CX_EQ(...) VARARG(ASSERT_CX_EQ,	 __VA_ARGS__)

#define ASSERT_CX_EQ_2(cx1, cx2) \
	ASSERT_NEAR(cx1.real(), cx2.real(), TOL); \
	ASSERT_NEAR(cx1.imag(), cx2.imag(), TOL);

#define ASSERT_CX_EQ_3(cx1, cx2, errpipe) \
	ASSERT_NEAR(cx1.real(), cx2.real(), TOL) << errpipe; \
	ASSERT_NEAR(cx1.imag(), cx2.imag(), TOL) << errpipe;

#define ASSERT_CX_EQ_4(cx1, cx2, errpipe, tol) \
	ASSERT_NEAR(cx1.real(), cx2.real(), tol) << errpipe; \
	ASSERT_NEAR(cx1.imag(), cx2.imag(), tol) << errpipe;

/*
*	Error stream must be terminated by "Eend"
* ErrSS << 23 << "dudulu" << Eend
*/
inline void ASSERT_MAT(const MatrixXcf& m1, const MatrixXcf& m2, const string& errstr = "", float tol = TOL)
{
	size_t r, c;
	ASSERT_EQ(r = m1.rows(), m2.rows()) << "Row dims should agree";
	ASSERT_EQ(c = m1.cols(), m2.cols()) << "Col dims should agree";

	CX a1, a2;
	for (size_t i = 0; i < r; ++i)
	for (size_t j = 0; j < c; ++j)
	{
		a1 = m1(i, j); a2 = m2(i, j);
		// ignore intellisense error here
		ASSERT_CX_EQ(a1, a2, 
					 "Disagree at [" << i << ", " << j << "]: " 
					 << m1(i, j) << " vs " << m2(i, j) << endl << errstr, tol);
	}
}

///////************** Generate random qubits **************///////
inline Qureg rand_qureg_dense(int nqubit, float symm)
{
	Qureg q = Qureg::create<true>(nqubit);
	DENSE_ITER(base)
		q.set_base_d(base, rand_cx(symm));
	return q;
}

/*
 *	'sparse': how many amplitudes to fill out
 */
inline Qureg rand_qureg_sparse(int nqubit, size_t sparse, float symm, bool forwardFill = true)
{
	size_t total = 1 << nqubit;
	if (sparse == 0)
		sparse = total / 2 + 1;
	Qureg qs = Qureg::create<false>(nqubit, sparse);
	if (forwardFill)
		for (qubase base : Range<qubase>(sparse))
			qs.add_base(base, rand_cx(symm));
	else
		for (qubase base : Range<qubase, false>(total, total - sparse))
			qs.add_base(base, rand_cx(symm));
	return qs;
}

/*
 *	Fill half of the sparse amplitudes
 */
inline size_t half_fill(int nqubit)
{
	return (1 << nqubit) / 2 + 1;
}

#endif // tests_h__
