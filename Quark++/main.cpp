#include "qureg.h"
#include "qumat.h"
#include "qugate.h"
#include "vld.h"
using namespace Qumat;
using namespace Qugate;

// For testing
// Dense init
Qureg qureg1 = Qureg::create<true>(2);
// Sparse init with no base at start
Qureg qureg3 = Qureg::create<false>(2, 4);
// Dense init
Qureg qureg4 = Qureg::create<true>(2);

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

	mm = MatrixXcf::Zero(3, 3);
	RowVectorXcf avec2(1); avec2 << 1;
	VectorXcf avec = VectorXcf::Zero(3);
	avec.rowwise() += avec2;
	pr(avec);
	mm.row(0) = avec.transpose();
	pr(mm);

	MatrixXcf v1(3, 2);
	v1 << 3, -2, 
		4, -1, 
		1, 0;
	MatrixXcf m1(2, 3);
	m1 << 10, 20, 30,
		40, 50, 60;
	pr(kronecker_mat(v1, m1));

	ptitle("Toffoli");
	pr(toffoli_mat());
	pr(cnot_mat());
	pr(cnot_mat());

	return 0;
}