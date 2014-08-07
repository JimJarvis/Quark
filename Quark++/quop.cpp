/**********************************************
* Quantum operations  *
**********************************************/
#include "quop.h"
#include "qureg.h"

using namespace Quop;

Qureg Quop::kronecker(Q1, Q2, bool resultDense)
{
	int new_nqubit = q1.nqubit + q2.nqubit;
	size_t size1 = q1.size();
	size_t size2 = q2.size();

	Qureg qans = resultDense ?
		Qureg::create<true>(new_nqubit) :
		Qureg::create<false>(new_nqubit, size1 * size2);

	// If both of them are dense, then the resultant is also dense
	int ri = 0;
	int nqubit2 = q2.nqubit;
	for (int i1 = 0; i1 < size1; ++i1)
		for (int i2 = 0; i2 < size2; ++i2)
		{
			qubase newBase = (q1.get_base<true>(i1) << nqubit2) | q2.get_base<true>(i2);
			CX newAmp = q1.amp[i1] * q2.amp[i2];
			if (resultDense)
				qans.amp[newBase] = newAmp;
			else
				qans.add_base(newBase, newAmp);
		}

	return qans;
}

Qureg Quop::operator*(Q1, Q2)
{
	return kronecker(q1, q2, q1.dense && q2.dense);
}
