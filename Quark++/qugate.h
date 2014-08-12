#ifndef qugate_h__
#define qugate_h__

#include "qureg.h"

namespace Qugate
{
	void generic_gate(Q, Matrix2cf&, int tar);

	void generic_gate(Q, Matrix4cf&, int tar1, int tar2);

	void generic_gate(Q, MatrixXcf&, vector<int>& tars);

	void generic_control(Q, Matrix2cf&, int ctrl, int tar);
	void cnot(Q, int ctrl, int tar);

	void generic_toffoli(Q, Matrix2cf&, int ctrl1, int ctrl2, int tar);
	void toffoli(Q, int ctrl1, int ctrl2, int tar);

	void generic_ncontrol(Q, Matrix2cf&, vector<int>& ctrls, int tar);
	void ncnot(Q, vector<int>& ctrls, int tar);


	void hadamard(Q, int tar);
	void hadamard(Q);
}

#endif // qugate_h__
