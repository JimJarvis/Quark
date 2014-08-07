/**********************************************
* Quantum operations  *
**********************************************/
#include "quop.h"
#include "qureg.h"

using namespace Quop;

Qureg Quop::kronecker(Q1, Q2)
{
	int new_nqubit = q1.nqubit + q2.nqubit;
	size_t size1 = q1.size();
	size_t size2 = q2.size();
	bool sparse = !q1.dense || !q2.dense;

	Qureg qans = sparse ? 
				Qureg::create<false>(new_nqubit, size1 * size2) :
				Qureg::create<true>(new_nqubit);

	// If both of them are dense, then the resultant is also dense
	int ri = 0;
	int nqubit2 = q2.nqubit;
	for (int i1 = 0; i1 < size1; ++i1)
		for (int i2 = 0; i2 < size2; ++i2)
		{
			qans.amp[ri] = q1.amp[i1] * q2.amp[i2];
			if (sparse)
				qans.get_base(ri) = 
					(q1.get_base<true>(i1) << nqubit2) | q2.get_base<true>(i2);
			++ri;
		}

	return qans;
}