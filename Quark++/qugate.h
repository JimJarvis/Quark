#ifndef qugate_h__
#define qugate_h__

#include "qureg.h"

namespace Qugate
{
	_DENSE_ void generic_gate(Q, Matrix2cf& mat, int tar);
	_DENSE_ void generic_gate(Q, Matrix4cf& mat, int tar1, int tar2);

	_DENSE_ void hadamard(Q, int tar);
	_DENSE_ void hadamard(Q);

	_DENSE_ void cnot(Q, int ctrl, int tar);
	_DENSE_ void toffoli(Q, int ctrl1, int ctrl2, int tar);

}

#endif // qugate_h__
