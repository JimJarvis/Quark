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
	while (1)
	{
		int M;
		cout << "\nEnter an integer to factorize: ";
		cin >> M;

		auto ans = shor_factorize(std::ceil(log2(M)), M, true);
		pr("\nFactorize " << M << " = " << ans.first << " * " << ans.second);
	}
	return 0;
}