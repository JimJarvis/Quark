/**********************************************
* Quantum Register  *
**********************************************/

#include "qureg.h"

#define BIT_PRINT

#ifdef BIT_PRINT
#define PRINT_KET(ket) bits2str<5>(ket)
#else
#define PRINT_KET(ket) (ket)
#endif // BIT_PRINT


Qureg::operator string()
{
	ostringstream oss;
	oss << setprecision(3) << "Qureg[";
	for (int i = 0; i < size ; ++i)
	{
		oss << "|" << PRINT_KET(isSparse() ? basis[i] : i) << "> ";
		CX a = amp[i];
		oss << a.real() << "+"
			<< a.imag() << "i"
			<< " (" << abs(a) << ")";
		if (i != size - 1)
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
	if (isDense())
	{
		size = 1 << nqubit;
		amp.resize(size, CX(0));
	}
	return *this;
}