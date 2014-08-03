#include "qureg.h"
#include <Eigen/Dense>
#include "vld.h"
using namespace Testing;
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

int main(int argc, char **argv)
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

	vector<int> vec;
	vec.reserve(30);
	pr(vec.size());

	return 0;
}