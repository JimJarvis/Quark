#include "qureg.h"
#include "quop.h"
#include "qugate.h"
#include "vld.h"
using namespace Testing;
using namespace Quop;
using namespace Qugate;
using namespace Eigen;

// For testing
// Dense init
Qureg qureg1 = Qureg::create<true>(2);
// Sparse init with no base at start
Qureg qureg3 = Qureg::create<false>(2, 4);
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

// Create new Qureg with dummy amplitude values
Qureg dummy_amp(int nqubit, bool dense)
{
	Qureg q = dense ?
		Qureg::create<true>(nqubit) :
		Qureg::create<false>(nqubit, 1);

	for (qubase base = 0; base < 1<<nqubit; ++base)
		if (dense)
			q.set_base_d(base, (base+1) * 10);
		else
			q.add_base(base, (base+1) * 10);
	return q;
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
	pause;
	ptitle("Sparse reg single");
	pr(qureg2);
	pause;
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

template<bool dense>
void test_hadamard()
{
	int nqubit = 3;
	Qureg q;
	for (size_t qi = 0; qi < 1 << nqubit; ++qi)
	{
		if (dense)
			q = Qureg::create<true>(nqubit, qi);
		else
			q = Qureg::create<false>(nqubit, 1, qi);
		hadamard(q);
		pr(q.sort());
	}
}

void test_cnot()
{
	int nqubit = 3;
	Qureg q;
	for (int qi = 0; qi < 1<<nqubit ; ++qi)
	{
		q = Qureg::create<true>(nqubit, qubase(qi));
		cnot(q, 1, 2);
		pr(q);
		q = Qureg::create<false>(nqubit, 1, qubase(qi));
		cnot(q, 1, 2);
		pr(q.sort().purge());
	}
}


int main(int argc, char **argv)
{
	init();
	pr((qureg1 * qureg3).sort());
	//eigen_demo();
	pr("dense hadamard");
	test_hadamard<true>();
	pr("sparse hadamard");
	test_hadamard<false>();
	//test_cnot();

	Qureg qq = dummy_amp(2, true);
	cnot(qq, 0, 1);
	pr(qq);

	qq = dummy_amp(2, false);
	Matrix2cf m;
	m <<
		0, 1,
		1, 0;
	generic_control(qq, m, 0, 1);
	pr(qq);

	ptitle("start vectoriong");
	vector<int> a;
	for (int i : Range<int, false>(10))
	{
		a.push_back(i * 3);
	}
	pr(vec2str(a));

	a.clear();
	for (int i : Range<int>(10))
	{
		a.push_back(i * 3);
	}
	pr(vec2str(a));

	for (int i : VecRange<int>(a, 3, 7))
	//for (int i : VecRange<int, false>(a, 1))
	{
		//i = 3;

		a.push_back(i * 10);
	}
	pr(vec2str(a));

	MatrixXcf mm(8, 8);
	MatrixXcf zz = MatrixXcf::Zero(6, 2);
	// toffoli gate
	mm << MatrixXcf::Identity(6, 6), zz, zz.transpose(), MatrixXcf::Identity(2,2).colwise().reverse();
	ptitle("toffoli");
	pr(mm);

	mm = MatrixXcf::Zero(3, 3);
	RowVectorXcf avec2(1); avec2 << 1;
	VectorXcf avec = VectorXcf::Zero(3);
	avec.rowwise() += avec2;
	pr(avec);
	mm.row(0) = avec.transpose();
	pr(mm);

	pr(hadamard_mat(3));

	VectorXcf v1(3);
	v1 << 3, -2, 4;
	MatrixXcf m1(2, 3);
	m1 << 10, 20, 30,
		40, 50, 60;
	pr(kronecker_mat(m1, v1));

	return 0;
}