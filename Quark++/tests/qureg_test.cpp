#include "tests.h"

TEST(Hello, World)
{
	Qureg q = Qureg::create<true>(3);
	pr(q);
}