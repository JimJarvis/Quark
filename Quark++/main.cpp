#include "qureg.h"
#include "quop.h"
#include "qugate.h"
#include "vld.h"
using namespace Testing;
using namespace Quop;
using namespace Qugate;

// For testing
// Dense init
Qureg qureg1(2);
// Sparse init with only 1 base at start
Qureg qureg3(2, 3, true, qubase(2));
// Dense init
Qureg qureg4(2);

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
	Qureg qureg1(3);
	// Sparse init with only 1 base at start
	Qureg qureg2(3, 2, true, qubase(8));

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
	for (int qi = 0; qi < 1<<nqubit ; ++qi)
	{
		Qureg q(nqubit, qubase(qi), false);
		hadamard(q);
		pr(q);
	}
}

int main(int argc, char **argv)
{
	init();
	pr(qureg1 * qureg4);
	//eigen_demo();
	//dense_hadamard();

	unordered_map<int, int> dud;
	dud[2] = 3;
	pr((dud.find(2))->second);

	return 0;
}