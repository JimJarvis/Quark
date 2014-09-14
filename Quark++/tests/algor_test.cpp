#include "tests.h"

TEST(Algor, DeutschJosza)
{
	for (int nbit : Range<>(1, 7))
		for (uint64_t secret_u : Range<uint64_t>(1 << nbit))
			// dense and sparse mode
			for (int dense : Range<>(2))
			{
				uint64_t result = deutsch_josza_parity(nbit, secret_u, dense);
				ASSERT_EQ(secret_u, result)
					<< "Inconsistency: secret_u = " << secret_u << "  result = " << result 
					<< " at nbit " << nbit << "; dense " << bool(dense);
			}
}

TEST(Algor, Simon)
{
	for (int nbit : Range<>(2, 6))
		for (uint64_t period : Range<uint64_t>(1, 1 << nbit))
			for (int dense : Range<>(2))
			{
				auto simon_result = simon_period(nbit, period, dense);
				Qureg q = move(simon_result.first);
				VectorXcf state(q);
				// period * y = 0 and each y should be of equal probability
				unordered_map<uint64_t, float> results;
				for (uint64_t base : Range<uint64_t>(state.size()))
				{
					uint64_t y = base >> nbit; // shift out the output bits
					float prob = norm(state[base]);
					if (contains(results, y))
						results[y] += prob;
					else
						if (prob > TOL)
							results[y] = prob;
				}

				// should all be the same
				float equal_prob = 0;
				for (auto& entry : results)
				{
					ASSERT_EQ(bitwise_dot(entry.first, period), 0) << "period * y == 0";
					if (equal_prob == 0)
						equal_prob = entry.second;
					else
						ASSERT_EQ(equal_prob, entry.second);
				}
			}
}