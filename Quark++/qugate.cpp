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
// Helper for generic gate
INLINE void generic_dense_update(
	vector<CX>& amp, qubase& base0, qubase& t, Matrix2cf& mat)
{
	if (!(base0 & t))
	{
		qubase base1 = base0 ^ t;
		CX a0 = amp[base0];
		CX a1 = amp[base1];
		amp[base0] = a0 * mat(0, 0) + a1 * mat(0, 1);
		amp[base1] = a0 * mat(1, 0) + a1 * mat(1, 1);
	}
}

INLINE void generic_sparse_update(
	Q, qubase& base0, qubase& t, Matrix2cf& mat)
{
	qubase base1 = base0 ^ t;
	CX a0, a1;
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

// Helper for cnot family and pauli_X
INLINE void cnot_sparse_update(Q, qubase& base, qubase& t)
{
	qubase base1 = base ^ t;
	if (q.contains_base(base1))
	{
		/* don't flip (swap) twice */
		if (base & t)
			std::swap(q[base], q[base1]);
	}
	else
	{
		q.add_base(base1, q[base]);
		q[base] = 0;
	}
}

void Qugate::generic_gate(Q, Matrix2cf& mat, int tar)
{
	qubase t = q.to_qubase(tar);
	if (q.dense)
	{
		auto& amp = q.amp;
		for (qubase base0 : q.base_iter_d())
			// only process base with 0 at the given target
			generic_dense_update(amp, base0, t, mat);
	}
	else // sparse
		// Add new states to the end, if any
		for (qubase base0 : q.base_iter())
			generic_sparse_update(q, base0, t, mat);
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

/*
 *	Explicitly expand the code for optimization
 */
void Qugate::pauli_X(Q, int tar)
{
	qubase t = q.to_qubase(tar);
	if (q.dense)
	{
		auto& amp = q.amp;
		for (qubase base0 : q.base_iter_d())
			// only process base with 0 at the given target
			if (!(base0 & t))
				swap(amp[base0], amp[base0 ^ t]);
	}
	else // sparse
		// Add new states to the end, if any
	for (qubase base0 : q.base_iter())
		cnot_sparse_update(q, base0, t);
}

void Qugate::pauli_Y(Q, int tar)
{
	generic_gate(q, pauli_Y_mat(), tar);
}

/*
*	Explicitly expand the code for optimization
*/
void Qugate::pauli_Z(Q, int tar)
{
	qubase t = q.to_qubase(tar);
	if (q.dense)
	{
		auto& amp = q.amp;
		for (qubase base0 : q.base_iter_d())
			// only process base with 0 at the given target
			if (!(base0 & t))
				amp[base0 ^ t] *= -1;
	}
	else // sparse
		// Add new states to the end, if any
	for (qubase base0 : q.base_iter())
	{
		qubase base1 = base0 ^ t;
		if (q.contains_base(base1))
		{
			if (!(base0 & t))
				q[base1] *= -1;
		}
		else
		{
			if (base0 & t)
				q[base0] *= -1;
		}
	}

}

void Qugate::rot_X(Q, float theta, int tar)
{
	generic_gate(q, rot_X_mat(theta), tar);
}

void Qugate::rot_Y(Q, float theta, int tar)
{
	generic_gate(q, rot_Y_mat(theta), tar);
}

void Qugate::rot_Z(Q, float theta, int tar)
{
	generic_gate(q, rot_Z_mat(theta), tar);
}

void Qugate::phase_scale(Q, float theta, int tar)
{
	const CX phase = expi(theta);
	qubase t = q.to_qubase(tar);
	if (q.dense)
	{
		auto& amp = q.amp;
		for (qubase base0 : q.base_iter_d())
			amp[base0] *= phase;
	}
	else // sparse
		for (qubase base0 : q.base_iter())
			q[base0] *= phase;
}

/**********************************************/
/*********** Multi-qubit gates  ***********/
/**********************************************/
void Qugate::generic_gate(Q, Matrix4cf& mat, int tar1, int tar2)
{
	qubase t1 = q.to_qubase(tar1);
	qubase t2 = q.to_qubase(tar2);
	Vector4cf a, newa;
	vector<qubase> basis(4);
	if (q.dense)
	{
		auto& amp = q.amp;
		for (qubase base0 : q.base_iter_d())
			// only process base with 00 at the given targets
			if (!(base0 & t1) && !(base0 & t2))
			{
				basis[0] = base0;
				basis[1] = base0 ^ t2;
				basis[2] = base0 ^ t1;
				basis[3] = base0 ^ (t1 | t2);
				for (int i = 0; i < 4; ++i)
					a(i) = amp[basis[i]];
				newa = mat * a;
				for (int i = 0; i < 4; ++i)
					amp[basis[i]] = newa(i);
			}
	}
	else // sparse
		// Not-so-efficient implementation: pretend to be dense
	for (qubase base0 : q.base_iter_d())
		if (!(base0 & t1) && !(base0 & t2))
		{
			basis[0] = base0;
			basis[1] = base0 ^ t2;
			basis[2] = base0 ^ t1;
			basis[3] = base0 ^ (t1 | t2);
			for (int i = 0; i < 4 ; ++i)
			{
				qubase base = basis[i];
				// if any of them doesn't exist, add
				if (!q.contains_base(base))
				{
					q.add_base(base, CX(0));
					a(i) = CX(0);
				}
				else
					a(i) = q[base];
			}
			newa = mat * a;
			for (int i = 0; i < 4; ++i)
				q[basis[i]] = newa(i);
		}
}

// Helper: get a vector of shifted qubase
INLINE vector<qubase> to_qubasis(Q, vector<int>& tars)
{
	vector<qubase> basis;
	basis.reserve(tars.size());
	for (int tar : tars)
		basis.push_back(q.to_qubase(tar));
	return basis;
}

// Helper for generic_gate: is the bits of 'base' at 'tar' positions all zero?
INLINE bool is_tar_all_zero(qubase& base, vector<qubase>& tarBasis)
{
	for (qubase& tar : tarBasis)
		if (base & tar) 
			return false;
	return true;
}

// Helper for generic_gate: bit positions of an int decides which target bits to flip
// fill out 'basis' parameter
INLINE void flipped_basis(vector<qubase>& basis, qubase& base0, vector<qubase>& tarBasis)
{
	// basis.size() is 2^n, where n = tarBasis.size()
	const int n = tarBasis.size();
	for (size_t i = 0; i < basis.size(); ++i)
	{
		qubase base = base0;
		for (int b = 0; b < n ; ++b)
			if (i & (1 << b))
				base ^= tarBasis[n-b-1]; // Most significant bit 
		basis[i] = base;
	}
}

void Qugate::generic_gate(Q, MatrixXcf& mat, vector<int>& tars)
{
	if (mat.rows() != 1 << tars.size())
		throw QuantumException(
			"Unitary matrix must have row/col width 2^(number of target bits)");

	int N = mat.rows();
	vector<qubase> tarBasis = to_qubasis(q, tars);
	VectorXcf a(N), newa(N);
	vector<qubase> basis(N);
	if (q.dense)
	{
		auto& amp = q.amp;
		for (qubase base0 : q.base_iter_d())
			// only process base with 00 at the given targets
		if (is_tar_all_zero(base0, tarBasis))
		{
			flipped_basis(basis, base0, tarBasis);
			for (int i = 0; i < N; ++i)
				a(i) = amp[basis[i]];
			newa = mat * a;
			for (int i = 0; i < N; ++i)
				amp[basis[i]] = newa(i);
		}
	}
	else // sparse
		// Not-so-efficient implementation: pretend to be dense
	for (qubase base0 : q.base_iter_d())
	if (is_tar_all_zero(base0, tarBasis))
	{
		flipped_basis(basis, base0, tarBasis);
		for (int i = 0; i < N; ++i)
		{
			qubase base = basis[i];
			// if any of them doesn't exist, add
			if (!q.contains_base(base))
			{
				q.add_base(base, CX(0));
				a(i) = CX(0);
			}
			else
				a(i) = q[base];
		}
		newa = mat * a;
		for (int i = 0; i < N; ++i)
			q[basis[i]] = newa(i);
	}
}

/**********************************************/
/*********** Multi-controlled gates  ***********/
/**********************************************/
void Qugate::cnot(Q, int ctrl, int tar)
{
	qubase c = q.to_qubase(ctrl);
	qubase t = q.to_qubase(tar);
	if (q.dense)
	{
		auto& amp = q.amp;
		for (qubase base : q.base_iter_d())
			if ((base & c) && (base & t)) // base & t: don't flip (swap)  twice
				std::swap(amp[base ^ t], amp[base]);
	}
	else // sparse
		// Add new states to the end, if any
		for (qubase base : q.base_iter())
			if (base & c)
				cnot_sparse_update(q, base, t);
}

void Qugate::generic_control(Q, Matrix2cf& mat, int ctrl, int tar)
{
	qubase c = q.to_qubase(ctrl);
	qubase t = q.to_qubase(tar);
	if (q.dense)
	{
		auto& amp = q.amp;
		for (qubase base : q.base_iter_d())
			if (base & c) // base & t: don't flip (swap)  twice
				generic_dense_update(amp, base, t, mat);
	}
	else // sparse
		// Add new states to the end, if any
		for (qubase base : q.base_iter())
			if (base & c)
				generic_sparse_update(q, base, t, mat);
}

void Qugate::toffoli(Q, int ctrl1, int ctrl2, int tar)
{
	qubase c1 = q.to_qubase(ctrl1);
	qubase c2 = q.to_qubase(ctrl2);
	qubase t = q.to_qubase(tar);
	if (q.dense)
	{
		auto& amp = q.amp;
		for (qubase base : q.base_iter_d())
		if ((base & c1) && (base & c2) && (base & t)) // base & t: don't flip (swap)  twice
			std::swap(amp[base ^ t], amp[base]);
	}
	else // sparse
		// Add new states to the end, if any
		for (qubase base : q.base_iter())
			if ((base & c1) && (base & c2))
				cnot_sparse_update(q, base, t);
}

void Qugate::generic_toffoli(Q, Matrix2cf& mat, int ctrl1, int ctrl2, int tar)
{
	qubase c1 = q.to_qubase(ctrl1);
	qubase c2 = q.to_qubase(ctrl2);
	qubase t = q.to_qubase(tar);
	if (q.dense)
	{
		auto& amp = q.amp;
		for (qubase base : q.base_iter_d())
		if ((base & c1) && (base & c2)) // base & t: don't flip (swap)  twice
			generic_dense_update(amp, base, t, mat);
	}
	else // sparse
		// Add new states to the end, if any
		for (qubase base : q.base_iter())
			if ((base & c1) && (base & c2))
				generic_sparse_update(q, base, t, mat);
}

// Helper: are the control bits on?
INLINE bool is_ctrl_on(qubase base, vector<qubase>& ctrlBasis)
{
	for (qubase& ctrl : ctrlBasis)
		if (!(base & ctrl))
			return false;
	return true;
}

void Qugate::ncnot(Q, vector<int>& ctrls, int tar)
{
	vector<qubase> ctrlBasis;
	ctrlBasis.reserve(ctrls.size());
	for (int ctrl : ctrls)
		ctrlBasis.push_back(q.to_qubase(ctrl));
	qubase t = q.to_qubase(tar);
	if (q.dense)
	{
		auto& amp = q.amp;
		for (qubase base : q.base_iter_d())
			if (is_ctrl_on(base, ctrlBasis) && (base & t))
				std::swap(amp[base ^ t], amp[base]);
	}
	else // sparse
		// Add new states to the end, if any
		for (qubase base : q.base_iter())
			if (is_ctrl_on(base, ctrlBasis))
				cnot_sparse_update(q, base, t);
}

void Qugate::generic_ncontrol(Q, Matrix2cf& mat, vector<int>& ctrls, int tar)
{
	vector<qubase> ctrlBasis = to_qubasis(q, ctrls);
	qubase t = q.to_qubase(tar);
	if (q.dense)
	{
		auto& amp = q.amp;
		for (qubase base : q.base_iter_d())
			if (is_ctrl_on(base, ctrlBasis))
				generic_dense_update(amp, base, t, mat);
	}
	else // sparse
		// Add new states to the end, if any
		for (qubase base : q.base_iter())
			if (is_ctrl_on(base, ctrlBasis))
				generic_sparse_update(q, base, t, mat);
}