#ifndef quop_h__
#define quop_h__

#include "qureg.h"
namespace Quop // quantum operations
{
	/*
	 *	Kronecker product
	 */
	Qureg kronecker(Q1, Q2);

	CX dot(Q1, Q2);

	void normalize(Q);
}

#endif // quop_h__
