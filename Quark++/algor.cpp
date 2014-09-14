/**********************************************
* Quantum Algorithms  *
**********************************************/

#include "algor.h"
#include "qureg.h"
#include "qumat.h"
#include "qugate.h"
using namespace Qumat;
using namespace Qugate;

uint64_t deutsch_josza_parity(int nbit, uint64_t secret_u, bool dense)
{
	if (secret_u > (1 << nbit))
		throw QuantumException("secret_u should not exceed 2^nbit");

	// return parity
	oracle_function oracle = [=](uint64_t x)
	{
		uint64_t dot = x & secret_u;
		uint64_t ans = 0;
		for (int bi = 0; bi < nbit ; ++bi)
			if (dot & (1 << bi))
				ans = !ans;
		return ans;
	};

	// output bit init to 1, all others 0
	Qureg q = dense ? 
		Qureg::create<true>(nbit + 1, qubase(1)) :
		Qureg::create<false>(nbit + 1, 1 << (nbit + 1), qubase(1));

	hadamard(q);

	apply_oracle(q, oracle, nbit);

	for (int qi = 0; qi < nbit; ++qi)
		hadamard(q, qi);

	// Randomly choose a measurement scheme
	if (rand_int(0, 2) == 0)
		// full measurement: discard the last output bit
		return measure(q) >> 1;
	else
	{
		// partial measurement: discard the last output bit
		uint64_t result = 0; 
		for (int qi = 0; qi < nbit; ++qi)
			// going from msb to lsb
			result |= measure(q, qi) << (nbit-1 - qi);
		return result;
	}
}