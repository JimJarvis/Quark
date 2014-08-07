#ifndef qureg_h__
#define qureg_h__

#include "utils.h"

/* convenient for function args */
#define Q Qureg<dense>& q
// dense
#define Qdense Qureg<true>& q
// sparse
#define Qsparse Qureg<false>& q
#define Q1 Qureg<dense1>& q1
#define Q2 Qureg<dense2>& q2
#define _DENSE_ template<bool dense>
#define _DENSE_2_ template<bool dense1, bool dense2>


///////***** Dummy type *****///////
template<bool dense, bool fullyDefined = true> class Qureg;

///////***** Common stuff between dense/sparse *****///////
template<bool dense>
class Qureg<dense, false>
{
public:
	int nqubit; // number of qubits
	vector<CX> amp; // amplitudes

	Qureg(int _nqubit, size_t ampSize) :
		nqubit(_nqubit),
		amp(vector<CX>(ampSize))
	{ }

	/*
	 *	Size of complex amplitude vector
	 */
	size_t size() { return amp.size(); }

	qubase get_base(size_t i);

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
	Qureg<dense>& operator+=(int scratchNqubit);

	///////***** Quop *****///////
	//friend Qureg operator*(Q1, Q2);
};

// more dummy
template<bool dense>
class Qureg<dense, true> : public Qureg<dense, false> {};

/**********************************************
* DENSE: we don't store basis vectors explicitly  *
**********************************************/
template<>
class Qureg<true> : public Qureg<true, false>
{
public:
	/*
	 * Init and set a specific base to amplitude 1. All other amps = 0.
	 */
	Qureg(int nqubit, qubase initBase = 0);

	qubase get_base(size_t i) { return i; }

	operator string();
	Qureg<true>& operator+=(int scratchNqubit);
};

/**********************************************
* SPARSE: we store basis vectors explicitly  *
**********************************************/
template<>
class Qureg<false> : public Qureg<false, false>
{
public:
	// non-zero basis, e.g. |00110> and |10100>
	// if basis is empty, we iterate over all 2^nqubit basis
	vector<qubase> basis; 
	// Maps a qubase to an index in amp[] array
	unordered_map<qubase, size_t> basemap;

	/*
	* If init true, we add initBase to amp[] with value 1
	* If init false, amp/basis[] will be empty and 'initBase' ignored
	* reservedSize for internal allocation
	* amp/basis will be empty
	*/
	Qureg(int nqubit, size_t reservedSize);

	/*
	* If init true, we add initBase to amp[] with value 1
	* If init false, amp/basis[] will be empty and 'initBase' ignored
	* reservedSize for internal allocation
	*/
	Qureg(int nqubit, size_t reservedSize, qubase initBase);

	qubase get_base(size_t i) { return basis[i]; }

	/*
	 *	Test if a base already exists in basis[]
	 */
	bool contains_base(qubase base)
	{
		return basemap.find(base) != basemap.end();
	}

	/*
	 *	Add a base. Processes hashmap
	 * Sparse ONLY. 
	 * Check true: if base already exists, update the amplitude
	 */
	template<bool check>
	void add_base(qubase base, CX a);
	// default: unchecked
	void add_base(qubase base, CX a) { add_base<false>(base, a); }

	/*
	 * Sparse ONLY: read index from basemap and get amplitude
	 */
	CX& amp_sparse(qubase base) { return amp[basemap[base]]; }

	operator string();
	Qureg<false>& operator+=(int scratchNqubit);

	///////***** Qugate *****///////
};


#endif // qureg_h__
