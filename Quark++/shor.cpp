// shor.cpp : Defines the entry point for the console application.
//

#include "qureg.h"
#include "qumat.h"
#include "qugate.h"
using namespace Qumat;
using namespace Qugate;


///////************** Helpers **************///////
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

int smallest_period(int b, int M)
{
	int i = 1;
	while (exp_mod(b, i, M) != 1) ++i;
	return i;
}

uint64_t long_pow(uint64_t a, int p)
{
	if (p == 1)
		return a;
	uint64_t partial = long_pow(a, p / 2);
	if (p % 2 == 1)
		return partial * partial * a;
	else
		return partial * partial;
}

int log2_int(int x)
{
	int ans = 0;
	while (x > 0)
	{
		x = x >> 1;
		++ans;
	}
	return ans;
}

/*
*	If size == 0, continue until 0
*/
vector<int64_t> to_continued_fraction(Frac frac, int size)
{
	vector<int64_t> cfrac;
	if (size > 0)
		cfrac.reserve(size);
	int i = 0;
	while (size < 1 || i < size)
	{
		cfrac.push_back(int64_t(frac));
		frac -= cfrac[i];
		if (frac.num == 0) break;
		frac = Frac(frac.denom, frac.num);
		++ i;
	}
	return cfrac;
}

int M = 17 * 13;
int nbit = log2_int(M) + 1;

// This is the user defined function that should be passed as an argument
uint64_t shor_oracle(uint64_t x)
{
	return exp_mod(nbit, x, M);
}

int main(int argc, char **argv)
{
	Qureg q0 = Qureg::create<true>(nbit * 2, 0);

	// top n qubits. This can be a user-library function
	qft(q0, 0, nbit);

	int b, i;
	while (1)
	{
		b = rand_int(2, M);

		if (gcd(b, M) != 1) continue;

		// Shor's circuit goes here
		Qureg q = q0.clone();

		// We can compile "shor_oracle" as a string and then get rid of the double quotes
		apply_oracle(q, shor_oracle, nbit);

		qft(q, 0, nbit);

		int mTrial = 0; // measurement trial
		int measured;

		while (mTrial++ < 10)
			// If 0, measure again
			// q ? nbit
		if ((measured = measure_top(q, nbit, false)) != 0)
		{
			// Use continued fraction approximation of {N / measured = r / k}
			vector<int64_t> cfrac = to_continued_fraction(Frac(1 << nbit, measured), 0);
			// We reduce the continued fraction more and more to get a simpler approximate fraction
			for (int size = cfrac.size(); size >= 1; --size)
			{
				// the actual period can be a multiple of p
				int p = to_frac(cfrac, size).num;
				int P = p;
				// 64-bit long long limit
				while (P < 128 && P < M)
				{
					if (P % 2 == 0
						&& exp_mod(b, P, M) == 1)
					{
						// further check  b^(p/2) != +/-1 mod M
						int check = exp_mod(b, P / 2, M);
						if (check != 1 && check != M - 1)
						{
							//pr("Almost there b = " << b << "; p = " << p << "; P = " << P 
							//   << "; cfrac = " << to_frac(cfrac, size) << " VS " << to_frac(cfrac));
							// We almost found it. Might have some numerical overflow here:
							uint64_t b_P_1 = long_pow(b, P / 2) - 1;
							int prime = gcd(M, b_P_1);
							// due to overflow, prime might become +/-1
							if (prime != 1 && prime != -1)
							{
								pr("Found period r = " << P);
								pr("b ^ r = " << b << " ^ " << P << " = 1 mod " << M);
								pr("b ^ (r/2) = " << b << " ^ " << P / 2 << " = " << check << " mod " << M);
								int prime2 = gcd(M, b_P_1 + 2); // b^(P/2) + 1
								pr("gcd(" << M << ", " << b_P_1 << ") = " << prime);
								pr("gcd(" << M << ", " << b_P_1 + 2 << ") = " << prime2);
								pr("\nFactorize " << M << " = " << prime << " * " << (prime2 == 1 ? M / prime : prime2));
								return 0;
							}
						}
					}
					// Try the next multiple of p, the next candidate of a possible period hit
					P += p;
				}
			}
		}
	}
	// We failed!!! ;(
	pr("FAIL");

	return 0;
}

