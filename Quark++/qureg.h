#ifndef qureg_h__
#define qureg_h__

#include "utils.h"

class Qureg
{
public:
	int nqubit; // number of qubits
	int size; // number of non-zeros in CX_VEC
	vector<CX> amp; // amplitudes
	// non-zero basis, e.g. |00110> and |10100>
	// if basis is empty, we iterate over all 2^nqubit basis
	vector<qubase> basis; 

	/*
	* Init to a specific basis state. All other basis have amp 0
	*/
	Qureg(int _nqubit, qubase initBase) : 
		nqubit(_nqubit), 
		size(1),
		amp(vector<CX>(1)),
		basis(vector<qubase>(1))
	{
		amp[0] = 1;
		basis[0] = initBase;
	}

	/*
	 *	Init to a dense register of n qubits with all amp = 0.
	 */
	Qureg(int _nqubit) :
		nqubit(_nqubit),
		size(1 << nqubit),
		amp(vector<CX>(size)),
		basis(vector<qubase>(0)) { }

	/*
	 *	Init to a sparse register of size N with all amp = 0
	 */
	Qureg(int _nqubit, int _size) :
		nqubit(_nqubit),
		size(_size),
		amp(vector<CX>(size)),
		basis(vector<qubase>(size)) { }

	/*
	 *	Convert to string
	 */
	operator string();

	friend ostream& operator<<(ostream& os, Qureg& qureg)
	{
		return os << string(qureg);
	}

	/*
	 *	True if we explicitly store the basis
	 */
	bool isSparse() { return !basis.empty(); }
};


#endif // qureg_h__
