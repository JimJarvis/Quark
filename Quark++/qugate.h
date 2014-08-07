#ifndef qugate_h__
#define qugate_h__

#include "qureg.h"

namespace Qugate
{
	void generic_gate(Q, Matrix2cf& mat, int tar);
	void generic_gate(Q, Matrix4cf& mat, int tar1, int tar2);

	void generic_control(Q, Matrix2cf& mat, int ctrl, int tar);
	void cnot(Q, int ctrl, int tar);
	void toffoli(Q, int ctrl1, int ctrl2, int tar);

	void hadamard(Q, int tar);
	void hadamard(Q);
}

#endif // qugate_h__
