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
	bool dense; // if we don't store basis[] explicitly
	vector<CX> amp; // amplitudes
	// non-zero basis, e.g. |00110> and |10100>
	// if basis is empty, we iterate over all 2^nqubit basis
	vector<qubase> basis; 
	// Maps a qubase to an index in amp[] array
	unordered_map<qubase, unsigned long> basemap;

	/*
	 *	Dense: we don't store the basis explicitly
	 * Init and set a specific base to amplitude 1. All other amps = 0.
	 */
	Qureg(int _nqubit, qubase initBase = 0) :
		nqubit(_nqubit),
		amp(vector<CX>(1 << nqubit)),
		dense(true)
	{
		amp[initBase] = 1;
	}

	/*
	* Sparse: we store the basis explicitly
	* If init true, we add initBase to amp[] with value 1
	* If init false, amp/basis[] will be empty and 'initBase' ignored
	* reservedSize for internal allocation
	*/
	Qureg(int _nqubit, size_t reservedSize, bool init, qubase initBase = 0) :
		nqubit(_nqubit),
		amp(vector<CX>()),
		basis(vector<qubase>()),
		dense(false)
	{
		amp.reserve(reservedSize);
		basis.reserve(reservedSize);
		basemap = unordered_map<qubase, unsigned long>(reservedSize);
		if (init)
		{
			basis.push_back(initBase);
			amp.push_back(CX(1));
			basemap[initBase] = 0; // indexed at 0
		}
	}

	/*
	 *	Size of complex amplitude vector
	 */
	int size() { return amp.size(); }

	/*
	 *	Test if a base already exists in basis[]
	 */
	bool contains_base(qubase base)
	{
		return basemap.find(base) == basemap.end();
	}

	/*
	 *	Add a base. Processes hashmap
	 * Sparse ONLY. If base already exists, update the amplitude
	 */
	void add_base(qubase base, CX a)
	{
		if (contains_base(base))
			amp[basemap[base]] = a;
		else
		{
			basis.push_back(base);
			amp.push_back(a);
		}
	}


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
	 *	If dense, get_base(i) == i
	 */
	qubase get_base(int i) { return dense ? i : basis[i]; }

	///////***** Quop *****///////
	friend Qureg operator*(Q1, Q2);

	///////***** Qugate *****///////
};


#endif // qureg_h__
