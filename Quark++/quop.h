#ifndef quop_h__
#define quop_h__

#include "qureg.h"
namespace Quop // quantum operations
{
	/*
	 *	Kronecker product
	 */
	Qureg kronecker(Qureg& q1, Qureg& q2);

	CX dot(Qureg& q1, Qureg& q2);

	void normalize(Qureg& q);
}

#endif // quop_h__
