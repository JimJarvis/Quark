/**********************************************
* Quantum gates  *
**********************************************/

#include "qureg.h"
#include "qugate.h"

using namespace Qugate;

void Qugate::generic_gate(Q, Matrix2cf& mat, int tar)
{
	qubase t = 1 << tar;
	CX a0, a1;
	qubase base1;
	if (q.dense)
	{
		auto& amp = q.amp;
		for (qubase base0 : q.base_iter_d())
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
	else // sparse
	{
		// Add new states to the end, if any
		for (qubase base0 : q.base_iter())
		{
			base1 = base0 ^ t;
			bool counterpart = q.contains_base(base1);
			// We always process this basis if it's 0 at target bit
			if (!(base0 & t))
			{
				a0 = q[base0];
				// we get the amplitude and don't add new base
				if (counterpart)
				{
					a1 = q[base1];
					q[base0] = a0 * mat(0, 0) + a1 * mat(0, 1);
					q[base1] = a0 * mat(1, 0) + a1 * mat(1, 1);
				}
				else // amp of base1 is 0 and we have to add new base
				{
					q[base0] = a0 * mat(0, 0);
					q.add_base(base1, a0 * mat(1, 0));
				}
			}
			// Otherwise base0 is target 1 and base1 is target 0
			// we process target 1 only if its 0 counterpart has not been processed
			else if (!counterpart)
			{
				a1 = q[base0];
				// a0 == 0
				q.add_base(base1, a1 * mat(0, 1));
				q[base0] = a1 * mat(1, 1);
			}
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
	for (int qi = 0; qi < q.nqubit; ++qi)
		hadamard(q, qi);
}

void Qugate::cnot(Q, int ctrl, int tar)
{
	qubase c = 1 << ctrl;
	qubase t = 1 << tar;
	qubase base1; // base1 is flipped base
	if (q.dense)
	{
		auto& amp = q.amp;
		for (qubase base : q.base_iter_d())
			if ((base & c) && (base & t)) // base & t: don't flip (swap)  twice
				std::swap(amp[base ^ t], amp[base]);
	}
	else // sparse
	{
		// Add new states to the end, if any
		for (qubase base : q.base_iter())
		{
			if (base & c)
			{
				base1 = base ^ t;
				if (q.contains_base(base1))
				{
					// don't flip (swap) twice
					if (base & t)
						std::swap(q[base], q[base1]);
				}
				else
				{
					q.add_base(base1, q[base]);
					q[base] = 0;
				}
			}
		}
	}
}

void Qugate::generic_control(Q, Matrix2cf& mat, int ctrl, int tar)
{
	qubase c = 1 << ctrl;
	qubase t = 1 << tar;
	qubase base1; // base1 is flipped base
	CX a, a1;
	if (q.dense)
	{
		auto& amp = q.amp;
		for (qubase base : q.base_iter_d())
			if ((base & c) && (base & t)) // base & t: don't flip (swap)  twice
			{
				base1 = base ^ t;
				a = amp[base];
				a1 = amp[base1];
				amp[base] = a * mat(0, 0) + a1 * mat(0, 1);
				amp[base1] = a * mat(1, 0) + a1 * mat(1, 1);
			}
	}
	else // sparse
	{
		// Add new states to the end, if any
		for (qubase base : q.base_iter())
		{
			if (base & c)
			{
				base1 = base ^ t;
				if (q.contains_base(base1))
				{
					// don't flip (swap) twice
					if (base & t)
					{
						a = q[base];
						a1 = q[base1];
						q[base] = a * mat(0, 0) + a1 * mat(0, 1);
						q[base1] = a * mat(1, 0) + a1 * mat(1, 1);
					}
				}
				else
				{
					a = q[base];
					q.add_base(base1, a * mat(1, 0));
					q[base] = a * mat(0, 0);
				}
			}
		}
	}
}