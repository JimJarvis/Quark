#include "qureg.h"
#include <Eigen/Dense>
#include "vld.h"
using namespace Testing;
using Eigen::Matrix2cf;

int main(int argc, char **argv)
{
	Qureg qureg(5);
	qureg = Qureg(3);
	qureg = Qureg(3, 5);

	pr(qureg);

	Matrix2cf m;
	m(0, 0) = CX(2,3);
	m(0, 1) = CX(3,4);
	m(1, 0) = CX(1,1);
	m(1, 1) = CX(2,3);
	pr(m);
	pr(m.adjoint());

	return 0;
}