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
	Qureg qureg(5);
	qureg = Qureg(3);
	qureg = Qureg(3, qubase(8));

	pr(qureg);

	vector<int> vec;
	vec.reserve(30);
	pr(vec.size());

	return 0;
}