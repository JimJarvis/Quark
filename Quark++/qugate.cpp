/**********************************************
* Quantum gates  *
**********************************************/

#include "qureg.h"
#include "qugate.h"
#include "qumat.h"

using namespace Qugate;
using namespace Qumat;

/**********************************************/
/*********** Single-qubit gates  ***********/
/**********************************************/

void Qugate::generic_gate(Q, Matrix2cf& mat, int tar)
{
	qubase t = q.to_bit(tar);
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
	generic_gate(q, hadamard_mat(), tar);
}

void Qugate::hadamard(Q)
{
	for (int qi = 0; qi < q.nqubit; ++qi)
		hadamard(q, qi);
}


/**********************************************/
/*********** Multi-controlled gates  ***********/
/**********************************************/
/*
	*	Macro to avoid code copying
	*/
#define ControlDenseNot \
	std::swap(amp[base ^ t], amp[base]);
#define ControlSparseNot \
{ \
	base1 = base ^ t; \
	if (q.contains_base(base1)) \
	{ \
		/* don't flip (swap) twice */ \
		if (base & t) \
			std::swap(q[base], q[base1]); \
	} \
	else \
	{ \
		q.add_base(base1, q[base]); \
		q[base] = 0; \
	} \
}
#define ControlDenseGeneric \
{ \
	base1 = base ^ t; \
	a = amp[base]; \
	a1 = amp[base1]; \
	amp[base] = a * mat(0, 0) + a1 * mat(0, 1); \
	amp[base1] = a * mat(1, 0) + a1 * mat(1, 1); \
}
#define ControlSparseGeneric \
{ \
	base1 = base ^ t; \
	if (q.contains_base(base1)) \
	{ \
		/* don't flip (swap) twice */ \
		if (base & t) \
		{ \
			a = q[base]; \
			a1 = q[base1]; \
			q[base] = a * mat(0, 0) + a1 * mat(0, 1); \
			q[base1] = a * mat(1, 0) + a1 * mat(1, 1); \
		} \
	} \
	else \
	{ \
		a = q[base]; \
		q.add_base(base1, a * mat(1, 0)); \
		q[base] = a * mat(0, 0); \
	} \
}

void Qugate::cnot(Q, int ctrl, int tar)
{
	qubase c = q.to_bit(ctrl);
	qubase t = q.to_bit(tar);
	qubase base1; // base1 is flipped base
	if (q.dense)
	{
		auto& amp = q.amp;
		for (qubase base : q.base_iter_d())
			if ((base & c) && (base & t)) // base & t: don't flip (swap)  twice
				ControlDenseNot
	}
	else // sparse
	{
		// Add new states to the end, if any
		for (qubase base : q.base_iter())
		{
			if (base & c)
				ControlSparseNot
		}
	}
}

void Qugate::generic_control(Q, Matrix2cf& mat, int ctrl, int tar)
{
	qubase c = q.to_bit(ctrl);
	qubase t = q.to_bit(tar);
	qubase base1; // base1 is flipped base
	CX a, a1;
	if (q.dense)
	{
		auto& amp = q.amp;
		for (qubase base : q.base_iter_d())
			if ((base & c) && (base & t)) // base & t: don't flip (swap)  twice
				ControlDenseGeneric
	}
	else // sparse
	{
		// Add new states to the end, if any
		for (qubase base : q.base_iter())
		{
			if (base & c)
				ControlSparseGeneric
		}
	}
}

void Qugate::toffoli(Q, int ctrl1, int ctrl2, int tar)
{
	qubase c1 = q.to_bit(ctrl1);
	qubase c2 = q.to_bit(ctrl2);
	qubase t = q.to_bit(tar);
	qubase base1; // base1 is flipped base
	if (q.dense)
	{
		auto& amp = q.amp;
		for (qubase base : q.base_iter_d())
		if ((base & c1) && (base & c2) && (base & t)) // base & t: don't flip (swap)  twice
			ControlDenseNot
	}
	else // sparse
	{
		// Add new states to the end, if any
		for (qubase base : q.base_iter())
		{
			if ((base & c1) && (base & c2))
				ControlSparseNot
		}
	}
}

void Qugate::generic_toffoli(Q, Matrix2cf& mat, int ctrl1, int ctrl2, int tar)
{
	qubase c1 = q.to_bit(ctrl1);
	qubase c2 = q.to_bit(ctrl2);
	qubase t = q.to_bit(tar);
	qubase base1; // base1 is flipped base
	CX a, a1;
	if (q.dense)
	{
		auto& amp = q.amp;
		for (qubase base : q.base_iter_d())
		if ((base & c1) && (base & c2) && (base & t)) // base & t: don't flip (swap)  twice
			ControlDenseGeneric
	}
	else // sparse
	{
		// Add new states to the end, if any
		for (qubase base : q.base_iter())
		{
			if ((base & c1) && (base & c2))
				ControlSparseGeneric
		}
	}
}