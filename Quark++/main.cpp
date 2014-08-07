#include "qureg.h"
#include "quop.h"
#include "qugate.h"
#include "vld.h"
using namespace Testing;
using namespace Quop;
using namespace Qugate;

// For testing
// Dense init
Qureg qureg1 = Qureg::create<true>(2);
// Sparse init with only 1 base at start
Qureg qureg3 = Qureg::create<false>(2, 3, qubase(2));
// Dense init
Qureg qureg4 = Qureg::create<true>(2);

void eigen_demo()
{
	Matrix2cf m;
	m << CX(2,3), CX(3,4), 
		CX(1,1), CX(2,3);
	//m(0, 0) = CX(2,3);
	//m(0, 1) = CX(3,4);
	//m(1, 0) = CX(1,1);
	//m(1, 1) = CX(2,3);
	pr(m);
	pr(m.adjoint());
	pr(m(0, 0) << m(0, 1) << m(1, 0) << m(1, 1));
}

void ctor()
{
	// Dense init
	Qureg qureg1 = Qureg::create<true>(3);
	// Sparse init with only 1 base at start
	Qureg qureg2 = Qureg::create<false>(3, 2, qubase(8));

	qureg1 += 3;
	qureg2 += 3;

	ptitle("Dense reg1");
	pr(qureg1);
	pause();
	ptitle("Sparse reg single");
	pr(qureg2);
	pause();
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

	qureg3.add_base(qubase(2), CX(1));
	qureg3.add_base(qubase(3), CX(2));
	qureg3.add_base(qubase(1), CX(3));
	qureg3.add_base(qubase(0), CX(4));
}

void dense_hadamard()
{
	int nqubit = 3;
	//int qi = 1;
	for (int qi = 0; qi < 1<<nqubit ; ++qi)
	{
		Qureg q = Qureg::create<false>(nqubit, 1, qubase(qi));
		hadamard(q);
		pr(q);
	}
}

int main(int argc, char **argv)
{
	init();
	//pr(qureg1 * qureg4);
	//eigen_demo();
	dense_hadamard();

	return 0;
}