/**********************************************
* Quantum operations  *
**********************************************/
#include "quop.h"
#include "qureg.h"

using namespace Quop;

Qureg Quop::kronecker(Q1, Q2)
{
	int new_nqubit = q1.nqubit + q2.nqubit;
	int size1 = q1.size;
	int size2 = q2.size;
	bool sparse = q1.isSparse() || q2.isSparse();

	Qureg qans = sparse ? 
				Qureg(new_nqubit, size1 * size2) :
				Qureg(new_nqubit);

	// If both of them are dense, then the resultant is also dense
	int ri = 0;
	int nqubit2 = q2.nqubit;
	for (int i1 = 0; i1 < size1; ++i1)
		for (int i2 = 0; i2 < size2; ++i2)
		{
			qans.amp[ri] = q1.amp[i1] * q2.amp[i2];
			if (sparse)
				qans.basis[ri] = 
					(q1.getBase(i1) << nqubit2) | q2.getBase(i2);
			++ri;
		}

	return qans;
}