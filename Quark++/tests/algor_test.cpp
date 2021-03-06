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

				uint64_t half_size = 1 << (nbit - 1);
				ASSERT_EQ(results.size(), half_size) << "Should have exactly half non-zero states";
				float equal_prob = 0;
				for (auto& entry : results)
				{
					ASSERT_EQ(bitwise_dot(entry.first, period), 0) << "period * y == 0";
					float expected_prob = 1.0 / half_size;
					ASSERT_NEAR(entry.second, expected_prob, TOL) << "probability should be " << expected_prob;
				}
			}
}

TEST(Algor, QftPeriod)
{
	for (int nbit : Range<>(4, 8))
		// N > 2 * r^2
		for (int period = 2; period < sqrt(1<<(nbit - 2)) ; ++period)
			for (int dense : Range<>(2))
			{
				Qureg q = qft_period(nbit, period, dense);
				auto sorted = q.sorted_non_zero_states();

				// Should be multiple of:  N / r  = 2^nbit / period
				// The good measured values will approach the following expected probability
				float expectProb = 1.0 / period;
				for (int i = 0; i < sorted.size(); ++i)
				{
					float prob = sorted[i].second;

					// This means the measurement will not be likely
					if (prob < expectProb * 0.7)
						break;

					qubase base = sorted[i].first >> nbit;

					// {base(measured) * r / N} should be as close to an integer as possible
					float k = (1.0 * base * period) / (1 << nbit);
					ASSERT_NEAR(k, round(k), 0.01)
						<< setprecision(3) << "Measured base = " << base << "\n1/period = " << 1.0 / period
						<< "\nprob = " << prob << "\nrank in the sorted bases = " << i;
				}
			}
}

TEST(Algor, Shor)
{
	// Each entry is a 3-tuple: (nbit, prime1, prime2)
	vector<vector<int>> trials =
	{
		{5, 3, 5},
		{6, 13, 5},
		{6, 11, 7},
		{7, 19, 3},
		{8, 17, 11},
		{9, 29, 13},
		//{9, 23, 19},
		{11, 31, 37}
		//{12, 41, 83}, // = 3403
	};

	for (auto& entry : trials)
	{
		int nbit = entry[0];
		int prime1 = entry[1];
		int prime2 = entry[2];
		int M = prime1 * prime2;

		pr("---------------");
		auto ans = shor_factorize(nbit, M, false);

		ASSERT_TRUE((prime1 == ans.first || prime1 == ans.second) && ans.first * ans.second == M)
			<< "M = " << M << " != " << ans.first << " * " << ans.second;

		pr("\nFactorize " << M << " = " << ans.first << " * " << ans.second);
	}
}

TEST(Algor, Grover)
{
	for (int nbit : Range<>(3, 8))
		for (int trial : Range<>(10))
			for (int dense : Range<>(2))
			{
				uint64_t key = rand_int(0, 1 << nbit);

				auto result = grover_search(nbit, key, dense);

				uint64_t foundKey = result.first;
				auto probHistory = result.second;

				ASSERT_EQ(key, foundKey) << " Key not found!";

				// The maximum probability should not be much less than 0.5
				// and should be near the sqrt(N) iteration
				int sqrtN = floor(sqrt(1 << nbit));

				int maxIdx = max_index(result.second);

				ASSERT_GT(probHistory[maxIdx], 0.40) 
					<< " Max probability should be at least 0.4";

				ASSERT_NEAR(maxIdx, sqrtN, 3) 
					<< "Max probability should occur near iteration " << sqrtN;
			}
}

TEST(Algor, Teleport)
{
	for (int trial : Range<>(10))
		for (int dense : Range<>(2))
		{
			// Create a random state in A register
			Qureg qa = dense ?
				Qureg::create<true>(1, 0) :
				Qureg::create<false>(1, 1<<1, 0);
			// random ops
			hadamard(qa);
			rot_X(qa, rand_float(PI / 12, PI / 2), 0);
			phase_shift(qa, rand_float(PI / 12, PI / 2), 0);

			VectorXcf qaAmp(qa);

			vector<CX> teleAmp = teleport(qa, dense);

			for (int i : Range<>(2))
			{
				ASSERT_CX_EQ(qaAmp(i), teleAmp[i], "Teleportation fails", 5e-7);
			}
		}
}