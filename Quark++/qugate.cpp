/**********************************************
* Quantum gates  *
**********************************************/

#include "qureg.h"
#include "qugate.h"

using namespace Qugate;

template<> void
Qugate::generic_gate(Qdense, Matrix2cf& mat, int tar)
{
	qubase t = 1 << tar;
	auto& amp = q.amp;
	CX a0, a1;
	qubase base0, base1;
	for (base0 = 0; base0 < q.size() ; ++base0)
		// only process base with 0 at the given target
		if (!(base0 & t))
		{
			base1 = base0 ^ t;
			a0 = amp[base0];
			a1 = amp[base1];
			amp[base0] = a0 * mat(0, 0) + a1 * mat(0, 1);
			amp[base1] = a0 * mat(1, 0) + a1 * mat(1, 1);
		}
}

template<> void
Qugate::generic_gate(Qsparse, Matrix2cf& mat, int tar)
{
	qubase t = 1 << tar;
	auto& amp = q.amp;
	CX a0, a1;
	qubase base0, base1;
	size_t oldSize = q.size();
	// Add new states to the end, if any
	for (int i = 0; i < oldSize; ++i)
	{
		base0 = q.basis[i];
		base1 = base0 ^ t;
		bool counterpart = q.contains_base(base1);
		// We always process this basis if it's 0 at target bit
		if (!(base0 & t))
		{
			a0 = q.amp_sparse(base0);
			// we get the amplitude and don't add new base
			if (counterpart)
			{
				a1 = q.amp_sparse(base1);
				q.amp_sparse(base0) = a0 * mat(0, 0) + a1 * mat(0, 1);
				q.amp_sparse(base1) = a0 * mat(1, 0) + a1 * mat(1, 1);
			}
			else // amp of base1 is 0 and we have to add new base
			{
				q.amp_sparse(base0) = a0 * mat(0, 0);
				q.add_base(base1, a0 * mat(1, 0));
			}
		}
		// Otherwise base0 is target 1 and base1 is target 0
		// we process target 1 only if its 0 counterpart has not been processed
		else if (!counterpart)
		{
			a1 = q.amp_sparse(base0);
			// a0 == 0
			q.add_base(base1, a1 * mat(0, 1));
			q.amp_sparse(base0) = a1 * mat(1, 1);
		}
	}
}

_DENSE_ void 
Qugate::hadamard(Q, int tar)
{
	static CX _sqrt2 = CX(1 / sqrt(2));
	static Matrix2cf HadamardMat;
	HadamardMat << 
		_sqrt2, _sqrt2, 
		_sqrt2, -_sqrt2;
	generic_gate(q, HadamardMat, tar);
}
template void Qugate::hadamard(Qdense, int tar);
template void Qugate::hadamard(Qsparse, int tar);

_DENSE_ void 
Qugate::hadamard(Q)
{
	for (int qi = 0; qi < q.nqubit; ++qi)
		hadamard(q, qi);
}
template void Qugate::hadamard(Qdense);
template void Qugate::hadamard(Qsparse);