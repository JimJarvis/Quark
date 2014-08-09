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

void Qureg::add_base(qubase base, CX a)
{
	basemap[base] = amp.size();
	basis.push_back(base);
	amp.push_back(a);
}

// Remove near-zero amplitudes
Qureg& Qureg::purge()
{
	if (dense) return *this; // do nothing
	vector<CX> purgedAmp;
	purgedAmp.reserve(amp.capacity());
	vector<qubase> purgedBasis;
	purgedBasis.reserve(basis.capacity());
	CX a;
	size_t s = 0; // new size
	for (qubase& base : basis)
	{
		a = (*this)[base];
		if (abs(a) > TOL)
		{
			purgedAmp.push_back(a);
			purgedBasis.push_back(base);
			basemap[base] = s++;
		}
		else // remove from basemap
			basemap.erase(base);
	}
	// update with new vectors
	amp = purgedAmp;
	basis = purgedBasis;
	return *this;
}

#define BIT_PRINT
#ifdef BIT_PRINT
#define PRINT_KET(ket) bits2str<3>(ket)
#else
#define PRINT_KET(ket) (ket)
#endif // BIT_PRINT

string Qureg::to_string(bool nonZeroOnly, bool bigEndian)
{
	ostringstream oss;
	oss << setprecision(3) << "Qureg[";
	size_t actualPrints = 0;
	for (size_t i = 0; i < size() ; ++i)
	{
		CX a = dense ? 
			amp[bigEndian ? to_bigend(i) : i] : 
			(*this)[get_base(i)];
		float magnitude = abs(a);
		if (nonZeroOnly && magnitude < TOL)
			continue;

		if (actualPrints != 0)
		{
			oss << ", ";
			if (actualPrints % 4 == 0)
				oss << "\n";
		}

		oss << "|" << PRINT_KET(
			(dense ? i : 
				(bigEndian ? 
					to_bigend(get_base(i)) : get_base(i)))
			) << "> ";
		oss << a.real() << "+"
			<< a.imag() << "i"
			<< " (" << magnitude << ")";

		++ actualPrints;
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