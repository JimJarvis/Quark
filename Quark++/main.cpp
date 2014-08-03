#include "qureg.h"
#include "vld.h"
using namespace Testing;

int main(int argc, char **argv)
{
	Qureg qureg(5);
	qureg = Qureg(5, qubase(14));

	pvec(qureg.amp);
	pvec(qureg.basis);

	return 0;
}