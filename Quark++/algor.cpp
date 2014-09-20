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

Qureg qft_period(int nbit, uint64_t period, bool dense /* = true */)
{
	// output bit init to 1, all others 0
	Qureg q = dense ? 
		Qureg::create<true>(nbit * 2, qubase(0)) :
		Qureg::create<false>(nbit * 2, 1 << nbit, qubase(0));

	qft(q, 0, nbit);

	apply_oracle(q, [=](uint64_t x){ return x % period + 1; }, nbit);

	// This measurement shouldn't really matter
	for (int tar = nbit + 1; tar < nbit * 2; ++tar)
		measure(q, tar);

	qft(q, 0, nbit);

	return q;
}

int shor_factorize(int nbit, int prime1, int prime2, bool dense)
{
	// output bit init to 1, all others 0
	Qureg q0 = dense ?
		Qureg::create<true>(nbit * 2, qubase(0)) :
		Qureg::create<false>(nbit * 2, 1 << nbit, qubase(0));

	// top n qubits
	qft(q0, 0, nbit);

	int M = prime1 * prime2;
	// randomly pick a base
	for (int b : Range<>(2, M / 2))
	{
		if (gcd(b, M) != 1)
			continue;

		Qureg q = q0.clone();

		apply_oracle(q, shor_oracle(b, M), nbit);

		// This measurement shouldn't really matter
		for (int tar = nbit + 1; tar < nbit * 2; ++tar)
			measure(q, tar);

		qft(q, 0, nbit);

	}
	return 0;
}

void shor_factorize_verbose(int nbit, int prime1, int prime2, bool dense)
{
	// output bit init to 1, all others 0
	Qureg q0 = dense ? 
		Qureg::create<true>(nbit * 2, qubase(0)) :
		Qureg::create<false>(nbit * 2, 1 << nbit, qubase(0));

	// top n qubits
	qft(q0, 0, nbit);

	int M = prime1 * prime2;
	// randomly pick a base
	for (int b : Range<>(2, M/2))
	{
		if (gcd(b, M) != 1)
			continue;

		Qureg q = q0.clone();

		apply_oracle(q, shor_oracle(b, M), nbit);

		// This measurement shouldn't really matter
		for (int tar = nbit + 1; tar < nbit * 2; ++tar)
			measure(q, tar);

		qft(q, 0, nbit);

		// verification
		int period = smallest_period(b, M);
		auto sorted = q.sorted_non_zero_states();

		pr("-------");
		pr("Smallest period: " << b << " ^ " << period << " = 1 mod " << M);
		float expectProb = 1.0 / period;
		pr("Measurement (m) should be multiple of " << (1 << nbit)*1.0 / period);
		pr("Expected prob 1/r = " << expectProb);
		pr("m\tm*r/N\t\tprob");
		for (int i = 0; i < sorted.size(); ++i)
		{
			int base = sorted[i].first >> nbit;
			float prob = sorted[i].second;

			// The measurement is unlikely to happen
			if (prob < expectProb * .7)  break;

			// measured * r / N  should be very close to the nearest integer
			pr(base << setprecision(5) << "\t" << (1.0*base*period / (1 << nbit)) << "\t\t" << prob);
		}
	}
}

///////************** Helpers **************///////
int gcd(int a, int b)
{
	int c;
	while (a != 0)
	{
		c = a;
		a = b % a;
		b = c;
	}
	return b;
}

uint64_t exp_mod(uint64_t b, uint64_t e, uint64_t m)
{
	int remainder;
	uint64_t x = 1;

	while (e != 0)
	{
		remainder = e % 2;
		e >>= 1;

		if (remainder == 1)
			x = (x * b) % m;
		b = (b * b) % m; // New base equal b^2 % m
	}
	return x;
}

oracle_function shor_oracle(int b, int M)
{
	return
		[=](uint64_t x)
		{
			return exp_mod(b, x, M);
		};
}

int smallest_period(int b, int M)
{
	int i = 1;
	while (exp_mod(b, i, M) != 1) ++i;
	return i;
}