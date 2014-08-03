#include "qureg.h"
#include "quop.h"
#include "qugate.h"
#include <Eigen/Dense>
#include "vld.h"
using namespace Testing;
using namespace Quop;
using namespace Qugate;
using Eigen::Matrix2cf;

// For testing
// Dense init
Qureg qureg1(2);
// Sparse init with only 1 base at start
Qureg qureg2(2, qubase(2));
// Sparse init with N
Qureg qureg3(2, 4);
// Dense init
Qureg qureg4(2);

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

void init()
{
	auto& amp1 = qureg1.amp;
	amp1[0] = CX(1, 0);
	amp1[1] = CX(2, 0);
	amp1[2] = CX(3, 0);
	amp1[3] = CX(4, 0);

	auto& amp4 = qureg4.amp;
	amp4[0] = CX(3, 0);
	amp4[1] = CX(4, 0);
	amp4[2] = CX(2, 0);
	amp4[3] = CX(1, 0);

	auto& amp3 = qureg3.amp;
	amp3[0] = CX(1, 0);
	amp3[1] = CX(2, 0);
	amp3[2] = CX(3, 0);
	amp3[3] = CX(4, 0);
	auto& basis3 = qureg3.basis;
	basis3[3] = qubase(0);
	basis3[2] = qubase(1);
	basis3[0] = qubase(2);
	basis3[1] = qubase(3);
}

int main(int argc, char **argv)
{
	init();

	pr(kronecker(qureg1, qureg4));

	return 0;
}