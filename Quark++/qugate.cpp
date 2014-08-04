/**********************************************
* Quantum gates  *
**********************************************/

#include "qureg.h"
#include "qugate.h"

using namespace Qugate;

void Qugate::generic_gate(Q, Matrix2cf mat, int tar)
{
	qubase t = 1 << tar;
	auto& amp = q.amp;
	CX a0, a1;
	qubase base1;
	if (q.isDense())
	{
		for (qubase base = 0; base < q.size() ; ++base)
			// only process base with 0 at the given target
			if (base & t)
			{
				base1 = base ^ t;
				a0 = amp[base];
				a1 = amp[base1];
				amp[base] = a0 * mat(0, 0) + a1 * mat(0, 1);
				amp[base1] = a0 * mat(1, 1) + a1 * mat(1, 0);
			}
	}
}

void Qugate::hadamard(Q, int tar)
{
	static CX _sqrt2 = CX(1 / sqrt(2));
	static Matrix2cf HadamardMat;
	HadamardMat << 
		_sqrt2, _sqrt2, 
		_sqrt2, -_sqrt2;
	generic_gate(q, HadamardMat, tar);
}

void Qugate::hadamard(Q)
{
	for (int qbit = 0; qbit < q.nqubit; ++qbit)
		hadamard(q, qbit);
}