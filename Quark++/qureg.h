#ifndef qureg_h__
#define qureg_h__

#include "utils.h"

/* convenient for function args */
#define Q Qureg& q
#define Q1 Qureg& q1
#define Q2 Qureg& q2

class Qureg
{
public:
	int nqubit; // number of qubits
	vector<CX> amp; // amplitudes
	// non-zero basis, e.g. |00110> and |10100>
	// if basis is empty, we iterate over all 2^nqubit basis
	vector<qubase> basis; 

	/*
	 *	Default dummy ctor
	 */
	Qureg() {}

	// Reserved memory for single-basis ctor
	static const int RSV_SIZE = 64;
	/*
	* Init to a single specific basis state that has amp 1. All other amps = 0.
	* Sparse if we don't store the basis explicitly, default true
	*/
	Qureg(int _nqubit, qubase initBase, bool sparse = true) : 
		nqubit(_nqubit)
	{
		amp = vector<CX>(sparse ? 1 : 1 << nqubit);
		basis = vector<qubase>();

		if (sparse)
		{
			amp.reserve(RSV_SIZE);
			basis.reserve(RSV_SIZE);
			amp.push_back(1);
			basis.push_back(initBase);
		}
		else
			amp[initBase] = 1;
	}

	/*
	 *	Init to a dense register of n qubits with |00..0> amp 1 and all other amp = 0.
	 */
	Qureg(int _nqubit)
	{
		*this = Qureg(_nqubit, 0, false);
	}


	/*
	 *	Init to a sparse register of size N with all amp = 0
	 */
	Qureg(int _nqubit, int _size) :
		nqubit(_nqubit),
		amp(vector<CX>(_size)),
		basis(vector<qubase>(_size)) { }

	/*
	 *	Size of complex amplitude vector
	 */
	int size() { return amp.size(); }

	/*
	 *	Convert to string
	 */
	operator string();

	friend ostream& operator<<(ostream& os, Q)
	{
		return os << string(q);
	}

	/*
	 *	Add scratch bits to 'this'. (Add to most significant bit)
	 */
	Qureg& operator+=(int scratch_nqubit);

	/*
	 *	True if we explicitly store the basis
	 */
	bool isSparse() { return !basis.empty(); }
	/*
	 *	True if we don't explicitly store the basis
	 * == !isSparse()
	 */
	bool isDense() { return basis.empty(); }

	/*
	 *	If dense, getBase(i) == i
	 */
	qubase getBase(int i) { return isDense() ? i : basis[i]; }

	///////***** Quop *****///////
	friend Qureg operator*(Q1, Q2);

	///////***** Qugate *****///////
};


#endif // qureg_h__
