#ifndef qureg_h__
#define qureg_h__

#include "utils.h"

/* convenient for function args */
#define Q Qureg& q
#define Q1 Qureg& q1
#define Q2 Qureg& q2

class Qureg
{
private:
	// Maps a qubase to an index in amp[] array
	unordered_map<qubase, size_t> basemap;
	// non-zero basis, e.g. |00110> and |10100>
	// if basis is empty, we iterate over all 2^nqubit basis
	vector<qubase> basis; 

	/*
	 * Private ctor
	 *	Dense: we don't store the basis explicitly
	 *    Set initBase to amplitude 1. All other amps = 0.
	 * Sparse: we store the basis explicitly
	 *    If init true, we add initBase to amp[] with value 1
	 *    If init false, amp/basis[] will be empty and 'initBase' ignored
	 *    reservedSize for internal allocation
	 */
	Qureg(bool dense, int nqubit, qubase initBase, size_t reservedSize, bool init);

public:
	int nqubit; // number of qubits
	bool dense; // if we don't store basis[] explicitly
	vector<CX> amp; // amplitudes

	/**********************************************
	* Creation ctors  *
	**********************************************/
	template<bool dense>
	static Qureg create(int nqubit, unsigned long long = 0);
	/*
	 *	Create a dense Qureg, amp[initBase] = 1 while all others 0
	 */
	template<>
	static Qureg create<true>(int nqubit, qubase initBase)
	{
		return Qureg(true, nqubit, initBase, 0, false);
	}
	/*
	 *	Create a sparse Qureg, all amps = 0. 
	 */
	template<>
	static Qureg create<false>(int nqubit, size_t reservedSize)
	{
		return Qureg(false, nqubit, 0, reservedSize, false);
	}
	template<bool dense>
	static Qureg create(int nqubit, size_t reservedSize, qubase initBase);
	/*
	 *	Create a dense Qureg, amp[initBase] = 1 while all others 0
	 */
	template<>
	static Qureg create<false>(int nqubit, size_t reservedSize, qubase initBase)
	{
		return Qureg(false, nqubit, initBase, reservedSize, true);
	}

	/**********************************************/
	/*********** Sparse ONLY  ***********/
	/**********************************************/
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
	 * 'checkExists' true: if base already exists, update the amplitude. Default false
	 */
	template<bool checkExists>
	void add_base(qubase base, CX a);
	// default: unchecked
	void add_base(qubase base, CX a) { add_base<false>(base, a); }

	/*
	 * Read index from basemap and get amplitude
	 */
	CX& operator[](qubase base) { return amp[basemap[base]]; }

	/**********************************************/
	/*********** Common part  ***********/
	/**********************************************/
	/*
	*	Size of complex amplitude vector
	*/
	size_t size() { return amp.size(); }

	/*
	 *	Get base stored at an internal index
	 * 'checkDense' true: checks whether we are dense or not. Default false
	 */
	template<bool checkDense>
	qubase& get_base(size_t i) { return checkDense && dense ? i : basis[i]; }
	qubase& get_base(size_t i) { return get_base<false>(i); }

	operator string();

	friend ostream& operator<<(ostream& os, Q)
	{
		return os << string(q);
	}

	/*
	*	Add scratch bits to 'this'. (Add to most significant bit)
	*/
	Qureg& operator+=(int scratchNqubit);

	///////************** Qugate **************///////
};


#endif // qureg_h__
