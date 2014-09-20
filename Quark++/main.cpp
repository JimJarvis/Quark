#include "qureg.h"
#include "qumat.h"
#include "qugate.h"
#include "algor.h"
using namespace Qumat;
using namespace Qugate;

int main(int argc, char **argv)
{
	ptitle("SHOR's ALGORITHM SIMULATION");
	pr("");
	bool mode = true;
	if (argc == 3)
		mode = str2int(argv[2]);

	int M = str2int(argv[1]);
	auto ans = shor_factorize(std::ceil(log2(M)), M, mode);
	pr("Factorize " << M << " = " << ans.first << " * " << ans.second);

	return 0;
}