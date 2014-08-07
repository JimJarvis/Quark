/**********************************************
* Quantum Register  *
**********************************************/

#include "qureg.h"
#include "quop.h"
#include "qugate.h"
using namespace Quop;
using namespace Qugate;

/*
* Dense: we don't store the basis explicitly
* Init and set a specific base to amplitude 1. All other amps = 0.
*/
Qureg<true>::Qureg(int _nqubit, qubase initBase) :
	Qureg<true, false>(_nqubit, 1 << nqubit)
{
	amp[initBase] = 1;
}

/*
* Sparse: we store the basis explicitly
* If init true, we add initBase to amp[] with value 1
* If init false, amp/basis[] will be empty and 'initBase' ignored
* reservedSize for internal allocation
*/
Qureg<false>::Qureg(int _nqubit, size_t reservedSize) :
	Qureg<false, false>(_nqubit, 0),
	basis(vector<qubase>()),
	basemap(unordered_map<qubase, size_t>(reservedSize))
{
	amp.reserve(reservedSize);
	basis.reserve(reservedSize);
}

Qureg<false>::Qureg(int _nqubit, size_t reservedSize, qubase initBase) :
	Qureg<false, false>(_nqubit, 1),
	basis(vector<qubase>(1)),
	basemap(unordered_map<qubase, size_t>(reservedSize))
{
	amp.reserve(reservedSize);
	basis.reserve(reservedSize);
	basis[0] = initBase;
	amp[0] = CX(1);
	basemap[initBase] = 0; // indexed at 0
}

template<bool check>
void Qureg<false>::add_base(qubase base, CX a)
{
	if (check && contains_base(base))
		amp[basemap[base]] = a;
	else
	{
		basemap[base] = amp.size();
		basis.push_back(base);
		amp.push_back(a);
	}
}
template void Qureg<false>::add_base<true>(qubase, CX);
template void Qureg<false>::add_base<false>(qubase, CX);


#define BIT_PRINT
#ifdef BIT_PRINT
#define PRINT_KET(ket) bits2str<4>(ket)
#else
#define PRINT_KET(ket) (ket)
#endif // BIT_PRINT

_DENSE_
Qureg<dense, false>::operator string()
{
	ostringstream oss;
	oss << setprecision(3) << "Qureg[";
	for (int i = 0; i < size() ; ++i)
	{
		oss << "|" << PRINT_KET(get_base(i)) << "> ";
		CX a = amp[i];
		oss << a.real() << "+"
			<< a.imag() << "i"
			<< " (" << abs(a) << ")";
		if (i != size() - 1)
		{
			oss << ", ";
			if (i % 4 == 3)
				oss << "\n";
		}
	}
	oss << "]";
	return oss.str();
}

template Qureg<true, false>::operator string();
template Qureg<false, false>::operator string();

//Qureg<false>::operator string()
//{
//	ostringstream oss;
//	oss << setprecision(3) << "Qureg[";
//	for (int i = 0; i < size() ; ++i)
//	{
//		oss << "|" << PRINT_KET(get_base(i)) << "> ";
//		CX a = amp[i];
//		oss << a.real() << "+"
//			<< a.imag() << "i"
//			<< " (" << abs(a) << ")";
//		if (i != size() - 1)
//		{
//			oss << ", ";
//			if (i % 4 == 3)
//				oss << "\n";
//		}
//	}
//	oss << "]";
//	return oss.str();
//}

Qureg<true>& Qureg<true>::operator+=(int scratch_nqubit)
{
	nqubit += scratch_nqubit;
	if (true)
		amp.resize(1 << nqubit, CX(0));
	return *this;
}
Qureg<false>& Qureg<false>::operator+=(int scratch_nqubit)
{
	nqubit += scratch_nqubit;
	if (false)
		amp.resize(1 << nqubit, CX(0));
	return *this;
}
//template Qureg<true>& Qureg<true>::operator+=(int);
//template Qureg<false>& Qureg<false>::operator+=(int);

///////***** Quop *****///////
//_DENSE_
//Qureg<dense> operator*(Q1, Q2)
//{
//	return kronecker(q1, q2);
//}