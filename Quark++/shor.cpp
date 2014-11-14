// shor.cpp : Defines the entry point for the console application.
//

#include "qureg.h"
#include "qumat.h"
#include "qugate.h"
using namespace Qumat;
using namespace Qugate;


///////************** Helpers **************///////
int64_t exp_mod(int64_t b, int64_t e, int64_t m)
{
	int64_t remainder;
	int64_t x = 1;

	while (e != 0)
	{
		remainder = e % 2;
		e = e >> 1;

		if (remainder == 1)
			x = (x * b) % m;
		b = (b * b) % m; // New base equal b^2 % m
	}
	return x;
}

int64_t smallest_period(int64_t b, int64_t M)
{
	int64_t i = 1;
	while (exp_mod(b, i, M) != 1) ++i;
	return i;
}

int64_t long_pow(int64_t a, int64_t p)
{
	if (p == 1)
		return a;
	int64_t partial = long_pow(a, p / 2);
	if (p % 2 == 1)
		return partial * partial * a;
	else
		return partial * partial;
}

int64_t log2_int(int64_t x)
{
	int64_t ans = 0;
	while (x > 0)
	{
		x = x >> 1;
		ans ++;
	}
	return ans;
}

/*
*	If size == 0, continue until 0
*/
vector<int64_t> to_continued_fraction(Frac frac, int64_t size)
{
	vector<int64_t> cfrac;
	int64_t i = 0;
	while (size < 1 || i < size)
	{
		cfrac.push_back(frac.num / frac.denom);
		frac -= cfrac[i];
		if (frac.num == 0) break;
		frac = ~frac;
		i ++;
	}
	return cfrac;
}

Frac to_fraction(const vector<int64_t>& cfrac, int64_t size)
{
	if (size < 1)
		size = cfrac.size();
	Frac ans(1, cfrac[size - 1]);
	for (int i = size - 2; i >= 1; --i)
	{
		ans += cfrac[i];
		ans = ~ans;
	}
	return ans + cfrac[0];
}

int64_t M = 17 * 13;
int64_t nbit = log2_int(M) + 1;

// This is the user defined function that should be passed as an argument
int64_t shor_oracle(int64_t x)
{
	return exp_mod(nbit, x, M);
}

int main(int argc, char **argv)
{
	Qureg q0 = Qureg::create<true>(nbit * 2, 0);

	// top n qubits. This can be a user-library function
	qft(q0, 0, nbit);

	int64_t b, i;
	while (true)
	{
		// built-in function: random int
		b = rand_int(2, M);

		if (gcd(b, M) != 1) continue;

		// Shor's circuit goes here
		Qureg q = q0.clone();

		// We can compile "shor_oracle" as a string and then get rid of the double quotes
		apply_oracle(q, shor_oracle, nbit);

		qft(q, 0, nbit);

		int64_t mTrial = 0; // measurement trial
		int64_t measured;

		while (mTrial < 10)
		{
			mTrial++;
			// If 0, measure again
			// q ? nbit
			measured = measure_top(q, nbit, false);
			if (measured != 0)
			{
				// Use continued fraction approximation of {N / measured = r / k}
				vector<int64_t> cfrac = to_continued_fraction(Frac(1 << nbit, measured), 0);
				// We reduce the continued fraction more and more to get a simpler approximate fraction
				for (int64_t size = cfrac.size(); size >= 1; --size)
				{
					// the actual period can be a multiple of p
					int64_t p = to_fraction(cfrac, size).num;
					int64_t P = p;
					// 64-bit long long limit
					while (P < 128 && P < M)
					{
						if (P % 2 == 0
							&& exp_mod(b, P, M) == 1)
						{
							// further check  b^(p/2) != +/-1 mod M
							int64_t check = exp_mod(b, P / 2, M);
							if (check != 1 && check != M - 1)
							{
								//pr("Almost there b = " << b << "; p = " << p << "; P = " << P 
								//   << "; cfrac = " << to_frac(cfrac, size) << " VS " << to_frac(cfrac));
								// We almost found it. Might have some numerical overflow here:
								int64_t b_P_1 = long_pow(b, P / 2) - 1;
								int64_t prime = gcd(M, b_P_1);
								// due to overflow, prime might become +/-1
								if (prime != 1 && prime != -1)
								{
									pr("Found period r = " << P);
									pr("b ^ r = " << b << " ^ " << P << " = 1 mod " << M);
									pr("b ^ (r/2) = " << b << " ^ " << P / 2 << " = " << check << " mod " << M);
									int64_t prime2 = gcd(M, b_P_1 + 2); // b^(P/2) + 1
									pr("gcd(" << M << ", " << b_P_1 << ") = " << prime);
									pr("gcd(" << M << ", " << b_P_1 + 2 << ") = " << prime2);
									int64_t other_prime;
									if (prime2 == 1)
										other_prime = M / prime;
									else
										other_prime = prime2;
									pr("\nFactorize " << M << " = " << prime << " * " << other_prime);
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
	}
	// We failed!!! ;(
	pr("FAIL");

	return 0;
}

