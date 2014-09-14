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
		/*uint64_t dot = x & secret_u;
		uint64_t ans = 0;
		for (int bi = 0; bi < nbit ; ++bi)
			if (dot & (1 << bi))
				ans = !ans;*/
		return bitwise_dot(x, secret_u);
	};

	// output bit init to 1, all others 0
	Qureg q = dense ? 
		Qureg::create<true>(nbit + 1, qubase(1)) :
		Qureg::create<false>(nbit + 1, 1 << (nbit + 1), qubase(1));

	hadamard(q);

	apply_oracle(q, oracle, nbit);

	hadamard_top(q, nbit);

	// Check consistency of measurement schemes
	// result from full measurement
	uint64_t result_full = measure_top(q, nbit, false);
	// result from partial measurement
	uint64_t result_partial = measure_top(q, nbit);

	if (result_full != result_partial)
		throw QuantumException("Full measurement doesn't agree with partial measurement");

	return result_full;
}

std::pair<Qureg, uint64_t> simon_period(int nbit, uint64_t period, bool dense)
{
	if (period > (1 << nbit))
		throw QuantumException("period should not exceed 2^nbit");

	// Let's make a table of this periodic function
	// map x to f(x)
	unordered_map<uint64_t, uint64_t> ftable(1 << (nbit - 1)); // half of 2^nbit 

	uint64_t rem = 1 << (nbit - 1); // remaining
	uint64_t fval = 0;
	uint64_t N = 1 << nbit;
	uint64_t xval = 0;
	// generate unique function values
	while (rem > 0)
	{
		if (rand_double() <= double(rem) / N)
		{
			while (contains(ftable, xval) || contains(ftable, xval ^ period))
				++ xval;
			//pr("xval " << xval << "; xval^ " << (xval ^ period) << " val = " << fval);
			ftable[xval] = fval;
			-- rem;
		}
		++ fval;
		-- N;
	}

	// capture ftable by ref and everything else by value
	oracle_function oracle = [=, &ftable](uint64_t x)
	{
		return ftable[contains(ftable, x) ? x : x ^ period];
	};

	// output bit init to 1, all others 0
	Qureg q = dense ? 
		Qureg::create<true>(nbit * 2, qubase(1)) :
		Qureg::create<false>(nbit * 2, 1 << nbit, qubase(1));

	hadamard_top(q, nbit);

	apply_oracle(q, oracle, nbit);

	hadamard_top(q, nbit);

	// non-destructive measurement
	// TODO post-classically solve for the period
	uint64_t result = measure_top(q, nbit, false);

	return pair<Qureg, uint64_t>(move(q), result);
}