#ifndef qureg_h__
#define qureg_h__

#include "utils.h"

/* convenient for function args */
#define Q Qureg& q
#define Q1 Qureg& q1
#define Q2 Qureg& q2

/*
 * Qubits start from least significant bit. 
 */
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

	/*
	 *	Dummy ctor for reference declaration
	 */
	Qureg() {}

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
	/*********** Dense ONLY  ***********/
	/**********************************************/
	INLINE void set_base_d(qubase base, CX a)
	{
		amp[base] = a;
	}

	/*
	 * Iterate over all basis from 0 to 1<<nqubit
	 */
	INLINE Range<qubase> base_iter_d() { return Range<qubase>(1 << nqubit); }

	/**********************************************/
	/*********** Sparse ONLY  ***********/
	/**********************************************/
	/*
	 *	Test if a base already exists in basis[]
	 */
	INLINE bool contains_base(qubase base)
	{
		return basemap.find(base) != basemap.end();
	}

	/*
	 *	Add a base. Processes hashmap. Sparse ONLY. 
	 */
	INLINE void add_base(qubase base, CX a)
	{
		basemap[base] = amp.size();
		basis.push_back(base);
		amp.push_back(a);
	}

	/*
	 *	Get base stored internally at index i
	 */
	INLINE qubase& get_base(size_t i) { return basis[i]; }

	/*
	 * Read index from basemap and get amplitude
	 */
	INLINE CX& operator[](qubase base) { return amp[basemap[base]]; }

	/*
	 *	For-each loop over basis[]. You can append to basis[] as you iterate
	 */
	INLINE VecRange<qubase> base_iter() { return VecRange<qubase>(basis); }

	/*
	 *	Sort the basis vectors. 
	 * The amp array is NOT sorted, cause basemap takes care of indices
	 */
	Qureg& sort()
	{
		if (!dense)
			std::sort(basis.begin(), basis.end());
		return *this;
	}

	/*
	 *	Remove near-zero amplitude, tolerance defined by TOL
	 */
	Qureg& purge();

	/**********************************************/
	/*********** Common part  ***********/
	/**********************************************/
	/*
	*	Size of complex amplitude vector
	*/
	INLINE size_t size() { return amp.size(); }

	/*
	 *	Get base stored at an internal index
	 * dense and sparse
	 */
	INLINE qubase& get_base_d_s(size_t i) { return dense ? i : basis[i]; }

	/*
	 *	Get a bit string representing a target qubit
	 * Most significant bit
	 */
	INLINE qubase to_qubase(int tar)
	{
		return qubase(1) << (nqubit - 1 - tar);
	}

	/*
	 *	If nonZeroOnly true, prints only states with non-zero amp
	 * default true
	 */
	string to_string(bool nonZeroOnly = true);

	operator string()
	{
		return to_string(true);
	}

	friend ostream& operator<<(ostream& os, Q)
	{
		return os << string(q);
	}

	/*
	*	Add scratch bits to 'this'. (Add to most significant bit)
	*/
	Qureg& operator+=(int scratchNqubit);

	/*
	 *	Convert to a column vector of 2^nqubit size
	 */
	operator VectorXcf();

	///////************** Qugate **************///////
};


#endif // qureg_h__
