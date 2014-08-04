/**********************************************
* Quantum Register  *
**********************************************/

#include "qureg.h"
#include "quop.h"
#include "qugate.h"
using namespace Quop;
using namespace Qugate;

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
		oss << "|" << PRINT_KET(isSparse() ? basis[i] : i) << "> ";
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
	if (isDense())
		amp.resize(1 << nqubit, CX(0));
	return *this;
}

///////***** Quop *****///////
Qureg operator*(Q1, Q2)
{
	return kronecker(q1, q2);
}