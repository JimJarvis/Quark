/**********************************************
* Quantum Register  *
**********************************************/

#include "qureg.h"
#include "quop.h"
#include "qugate.h"
using namespace Quop;
using namespace Qugate;

// Private comprehensive ctor
// The last two args are only relevant to sparse mode
Qureg::Qureg(bool _dense, int _nqubit, qubase initBase, size_t reservedSize, bool init) :
	dense(_dense),
	nqubit(_nqubit),
	amp(vector<CX>(dense ? 1 << nqubit : 0))
{
	if (dense)
		amp[initBase] = 1;
	else
	{
		amp.reserve(reservedSize);
		basis = vector<qubase>();
		basis.reserve(reservedSize);
		basemap = unordered_map<qubase, size_t>(reservedSize);
		if (init)
		{
			basis.push_back(initBase);
			amp.push_back(CX(1));
			basemap[initBase] = 0; // indexed at 0
		}
	}
}


template<bool checkExists>
void Qureg::add_base(qubase base, CX a)
{
	if (checkExists && contains_base(base))
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
		oss << "|" << PRINT_KET(get_base<true>(i)) << "> ";
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