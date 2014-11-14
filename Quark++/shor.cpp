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

int main(int argc, char **argv)
{
	int M = 17 * 13;
	int nbit = log2_int(M) + 1;

	Qureg q0 = Qureg::create<true>(nbit * 2, 0);

	// top n qubits. This can be a user-library function
	qft(q0, 0, nbit);

	// Don't try b's already tried
	unordered_set<int> bTriedSet;
	// If the hash set is more than 2/3 filled, we try all the rest values sequentially
	bool hashMode = true;
	vector<int> remainings; // b values we haven't tried yet
	int b, i;
	while (1)
	{
		if (hashMode && bTriedSet.size() < M * 2.0 / 3)
		{
			// the random picked base must be co-prime with M
			b = rand_int(2, M);

			if (contains(bTriedSet, b))
				continue;
			else
				bTriedSet.insert(b);
		}
		else
		{
			// init the vector once
			if (hashMode)
			{
				hashMode = false;
				i = 0;
				for (int btmp = 2; btmp < M; ++btmp)
				if (!contains(bTriedSet, btmp))
					remainings.push_back(btmp);
			}
			if (i < remainings.size())
				b = remainings[i++];
			else
				break; // we tried every b and failed. Break the infinite loop
		}

		if (gcd(b, M) != 1) continue;

		// Shor's circuit goes here
		Qureg q = q0.clone();

		apply_oracle(q, shor_oracle(b, M), nbit);

		qft(q, 0, nbit);

		int mTrial = 0; // measurement trial
		int measured;

		while (mTrial++ < 10)
			// If 0, measure again
		if ((measured = measure_top(q, nbit, false)) != 0)
		{
			// Use continued fraction approximation of {N / measured = r / k}
			ContFrac cfrac = to_cont_frac(Frac(1 << nbit, measured));
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

