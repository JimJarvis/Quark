#include "qureg.h"
#include "qumat.h"
#include "qugate.h"
#include "algor.h"
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

void test_toffoli()
{
	int nqubit = 3;
	Qureg q;
	for (int qi = 0; qi < 1 << nqubit; ++qi)
	{
		q = Qureg::create<true>(nqubit, qubase(qi));
		toffoli(q, 0, 1, 2);
		pr(q);
		q = Qureg::create<false>(nqubit, 1, qubase(qi));
		toffoli(q, 0, 1, 2);
		pr(q.sort().purge());
	}
}

void test_rand_unique()
{
	vector<int> histogram(100, 0);
	vector<int> holder(5, 0);
	for (int i : Range<>(1000000))
	for (size_t j : rand_unique<int>(holder, 5, 100))
		++histogram[j];
	pr(vec2str(histogram));
}

int main(int argc, char **argv)
{
	init();
	pr((qureg1 & qureg3).sort());
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
	MatrixXcf m1(2, 2);
	m1 << 10, 20,
		40, 50;

	Qureg q = Qureg::create<false>(3, 1, qubase(0));
	//q = Qureg::create<true>(3, qubase(0));
	hadamard(q);

	rand_seed(10);
	pr("printed " << measure(q, 0));
	pr(q);

	ptitle("ALGORITHMS");
	//for (int i = 0; i < 15 ; ++i)
	//	deutsch_josza_parity(4, i);
	auto pair_result = simon_period(4, 1, true);
	pr("Simon's " << pair_result.first << "\nPeriod = " << pair_result.second);

	pr(exp_mod(65489, 251025, 894603));

	//pr(smallest_period(2, 35));
	//pr(smallest_period(3, 35));
	//pr(smallest_period(4, 35));
	//shor_factorize_verbose(7, 7, 11);

	ptitle("Continued Frac");

	ContFrac vv = { 2, 3, -6, -1, 2, 7 };

	pr(to_frac(vv));

	pvec(to_cont_frac(to_frac(vv), vv.size()));

	pr(to_frac(to_cont_frac(Frac(8192, 2567), 5)));

	return 0;
}