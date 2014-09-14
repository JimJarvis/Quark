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

	// Check consistency of measurement schemes
	// result from full measurement
	uint64_t result_full = measure_top(q, nbit, false);
	// result from partial measurement
	uint64_t result_partial = measure_top(q, nbit);

	if (result_full != result_partial)
		throw QuantumException("Full measurement doesn't agree with partial measurement");

	return result_full;
}