#include "tests.h"

TEST(Measure, All)
{
	srand(time(0));
	Qureg q = Qureg::create<false>(6, 10, qubase(19));
	hadamard(q);
	pvec(q.amp);
	vector<int> hist(64); // histogram

	for (int i : Range<>(6400))
		++hist[measure(q)];

	pvec(hist);

}