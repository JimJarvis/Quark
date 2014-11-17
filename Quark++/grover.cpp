#include "qureg.h"
#include "qumat.h"
#include "qugate.h"
using namespace Qumat;
using namespace Qugate;

int64_t nbit = 7;
int64_t key = 73; // secret key

int64_t grover_oracle(int64_t x)
{
	return x == key;
}

int main(int argc, char ** argv)
{
	// the last output bit init to 1
	Qureg q = Qureg::create<true>(nbit + 1, 1);

	int64_t N = 1 << nbit;
	int64_t sqrtN = floor(sqrt(N));

	// Init by superposition
	hadamard(q);

	int64_t ans = 0;
	vector<float> probAtKey;

	// optimal: PI / 4 * sqrt(N) times
	int64_t iter;
	for (iter = 0; iter < sqrtN * 2; ++iter)
	{
		// phase inversion
		apply_oracle(q, grover_oracle, nbit);

		// mean inversion
		hadamard_top(q, nbit);
		grover_diffuse(q);
		hadamard_top(q, nbit);

		if (iter == sqrtN)
			// measure until we get the solution
			while (grover_oracle(ans) == 0)
				ans = measure_top(q, nbit, false);

		probAtKey.push_back( q.prefix_prob(nbit, key) );
	}

	pr("Found key: " << ans);
	pr("Probability = " << probAtKey);
}

