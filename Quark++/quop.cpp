/**********************************************
* Quantum operations  *
**********************************************/
#include "quop.h"
#include "qureg.h"

using namespace Quop;

Qureg Quop::kronecker(Qureg& reg1, Qureg& reg2)
{
	int new_nqubit = reg1.nqubit + reg2.nqubit;
	int size1 = reg1.size;
	int size2 = reg2.size;
	bool sparse = reg1.isSparse() || reg2.isSparse();

	Qureg reg = sparse ? 
				Qureg(new_nqubit, size1 * size2) :
				Qureg(new_nqubit);

	// If both of them are dense, then the resultant is also dense
	int ri = 0;
	int nqubit2 = reg2.nqubit;
	for (int i1 = 0; i1 < size1; ++i1)
		for (int i2 = 0; i2 < size2; ++i2)
		{
			reg.amp[ri] = reg1.amp[i1] * reg2.amp[i2];
			if (sparse)
				reg.basis[ri] = 
					(reg1.getBase(i1) << nqubit2) | reg2.getBase(i2);
			++ri;
		}

	return reg;
}