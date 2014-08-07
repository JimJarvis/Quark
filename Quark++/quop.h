#ifndef quop_h__
#define quop_h__

#include "qureg.h"
namespace Quop // quantum operations
{
	/*
	 *	Kronecker product
	 */
	template<bool denseReturn, bool dense1, bool dense2>
	Qureg<denseReturn> kronecker(Q1, Q2);

	_DENSE_2_
	CX dot(Q1, Q2);

	_DENSE_
	void normalize(Q);
}

#endif // quop_h__
