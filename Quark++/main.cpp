#include "qureg.h"
#include "quop.h"
#include "qugate.h"
#include <Eigen/Dense>
#include "vld.h"
using namespace Testing;
using namespace Quop;
using namespace Qugate;
using Eigen::Matrix2cf;

void eigen_demo()
{
	Matrix2cf m;
	m(0, 0) = CX(2,3);
	m(0, 1) = CX(3,4);
	m(1, 0) = CX(1,1);
	m(1, 1) = CX(2,3);
	pr(m);
	pr(m.adjoint());
}

void ctor()
{
	// Dense init
	Qureg qureg1(3);
	// Sparse init with only 1 base at start
	Qureg qureg2(3, qubase(8));
	// Sparse init with N
	Qureg qureg3(3, 8);

	qureg1 += 3;
	qureg2 += 3;
	qureg3 += 3;

	ptitle("Dense reg1");
	pr(qureg1);
	pause();
	ptitle("Sparse reg single");
	pr(qureg2);
	pause();
	ptitle("Sparse reg N");
	pr(qureg3);
}

int main(int argc, char **argv)
{
	// Dense init
	Qureg qureg1(2);
	// Sparse init with only 1 base at start
	Qureg qureg2(2, qubase(3));
	// Sparse init with N
	Qureg qureg3(3, 8);

	auto& amp = qureg1.amp;
	amp[0] = CX(0.2, 1);
	amp[1] = CX(0.3, -1);
	amp[2] = CX(0.1, .5);
	amp[3] = CX(0.5, -.5);
	pr(kronecker(qureg1, qureg2));

	return 0;
}