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
Qureg::Qureg(int _nqubit, qubase initBase) :
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
Qureg::Qureg(int _nqubit, size_t reservedSize, bool init, qubase initBase) :
	nqubit(_nqubit),
	amp(vector<CX>()),
	basis(vector<qubase>()),
	dense(false)
{
	amp.reserve(reservedSize);
	basis.reserve(reservedSize);
	basemap = unordered_map<qubase, size_t>(reservedSize);
	if (init)
	{
		basis.push_back(initBase);
		amp.push_back(CX(1));
		basemap[initBase] = 0; // indexed at 0
	}
}

template<bool Check>
void Qureg::add_base(qubase base, CX a)
{
	if (Check && contains_base(base))
		amp[basemap[base]] = a;
	else
	{
		basemap[base] = amp.size();
		basis.push_back(base);
		amp.push_back(a);
	}
}
template void Qureg::add_base<true>(qubase, CX);
template void Qureg::add_base<false>(qubase, CX);


#define BIT_PRINT
#ifdef BIT_PRINT
#define PRINT_KET(ket) bits2str<4>(ket)
#else
#define PRINT_KET(ket) (ket)
#endif // BIT_PRINT

Qureg::operator string()
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

Qureg& Qureg::operator+=(int scratch_nqubit)
{
	nqubit += scratch_nqubit;
	if (dense)
		amp.resize(1 << nqubit, CX(0));
	return *this;
}

///////***** Quop *****///////
Qureg operator*(Q1, Q2)
{
	return kronecker(q1, q2);
}