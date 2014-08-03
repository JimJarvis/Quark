#include "qureg.h"
#include "vld.h"
using namespace Testing;

int main(int argc, char **argv)
{
	Qureg qureg(5);
	qureg = Qureg(3);
	qureg = Qureg(3, 5);

	pr(qureg);

	return 0;
}