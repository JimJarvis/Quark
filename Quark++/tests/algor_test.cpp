#include "tests.h"

TEST(Algor, Deutsch)
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