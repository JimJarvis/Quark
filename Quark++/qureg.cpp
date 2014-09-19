/**********************************************
* Quantum Register  *
**********************************************/

#include "qureg.h"
#include "qumat.h"
#include "qugate.h"
using namespace Qumat;
using namespace Qugate;

// Private comprehensive ctor
// The last two args are only relevant to sparse mode
Qureg::Qureg(bool _dense, int _nqubit, qubase initBase, size_t reservedSize, bool init) :
	dense(_dense),
	nqubit(_nqubit),
	amp(vector<CX>(dense ? 1 << nqubit : 0))
{
	if (dense)
		amp[initBase] = 1;
	else
	{
		amp.reserve(reservedSize);
		basis = vector<qubase>();
		basis.reserve(reservedSize);
		basemap = unordered_map<qubase, size_t>(reservedSize);
		if (init)
		{
			basis.push_back(initBase);
			amp.push_back(CX(1));
			basemap[initBase] = 0; // indexed at 0
		}
	}
}

// Remove near-zero amplitudes
Qureg& Qureg::purge()
{
	if (dense) return *this; // do nothing
	vector<CX> purgedAmp;
	purgedAmp.reserve(amp.capacity());
	vector<qubase> purgedBasis;
	purgedBasis.reserve(basis.capacity());
	CX a;
	size_t s = 0; // new size
	for (qubase& base : basis)
	{
		a = (*this)[base];
		if (norm(a) > TOL)
		{
			purgedAmp.push_back(a);
			purgedBasis.push_back(base);
			basemap[base] = s++;
		}
		else // remove from basemap
			basemap.erase(base);
	}
	// update with new vectors
	amp = move(purgedAmp);
	basis = move(purgedBasis);
	return *this;
}

Qureg Qureg::clone()
{
	Qureg qc;
	qc.nqubit = this->nqubit;
	qc.dense = this->dense;
	qc.amp = this->amp;
	qc.basis = this->basis;
	qc.basemap = this->basemap;
	return qc;
}

#define BIT_PRINT
#ifdef BIT_PRINT
#define PRINT_KET(ket) bits2str(ket, nqubit)
#else
#define PRINT_KET(ket) (ket)
#endif // BIT_PRINT

string Qureg::to_string(bool nonZeroOnly)
{
	ostringstream oss;
	oss << setprecision(3) << "Qureg[";
	size_t actualPrints = 0;
	for (size_t i = 0; i < size(); ++i)
	{
		CX a = dense ? amp[i] : (*this)[get_base(i)];
		float prob = norm(a);
		if (nonZeroOnly && prob < TOL)
			continue;

		if (actualPrints != 0)
		{
			oss << ", ";
			if (actualPrints % 4 == 0)
				oss << "\n";
		}
		oss << "|" << PRINT_KET(get_base_d_s(i)) << "> ";
		oss << a.real() << "+"
			<< a.imag() << "i"
			<< " (" << prob << ")";

		++actualPrints;
	}
	oss << "]";
	return oss.str();
}

Qureg& Qureg::operator+=(int scratchQubits)
{
	if (dense)
	{
		amp.resize(1 << (nqubit + scratchQubits), CX(0));
		// Move old amplitudes over
		// Do it in the reverse order to avoid conflict
		for (size_t i = (1<<nqubit) - 1; i >= 0 ; --i)
		{
			amp[i << scratchQubits] = amp[i];
			amp[i] = 0;
		}
	}
	else
		// Update sparse and hashmap
		for (qubase& base : basis)
		{
			size_t idx = basemap[base];
			basemap.erase(base);
			base <<= scratchQubits;
			basemap[base] = idx;
		}
	nqubit += scratchQubits;
	return *this;
}

Qureg::operator VectorXcf()
{
	VectorXcf vec(1 << nqubit);
	if (dense)
		for (qubase base = 0; base < 1<<nqubit ; ++base)
			vec(base) = amp[base];
	else
		for (qubase base = 0; base < 1<<nqubit ; ++base)
			vec(base) = contains_base(base) ? (*this)[base] : 0;
	return vec;
}

vector<qubase> Qureg::non_zero_states()
{
	vector<qubase> nonZeros;
	if (dense)
		for (qubase base = 0; base < 1<<nqubit ; ++base)
			if (norm(amp[base]) > TOL)
				nonZeros.push_back(base);
	else
		for (qubase& base : basis)
			if (norm((*this)[base]) > TOL)
				nonZeros.push_back(base);
	return nonZeros;
}

vector<pair<qubase, float>> Qureg::sorted_non_zero_states()
{
	typedef pair<qubase, float> qentry;
	auto cmp = [](const qentry& x1, const qentry& x2)
	{
		return x1.second < x2.second;
	};
	std::priority_queue<qentry, vector<qentry>, decltype(cmp)> que(cmp);

	if (dense)
		for (qubase base = 0; base < 1 << nqubit; ++base)
		{
			float prob = norm(amp[base]);
			if (prob > TOL)
				que.push(qentry(base, prob));
		}
	else
		for (qubase& base : basis)
		{
			float prob = norm((*this)[base]);
			if (prob > TOL)
				que.push(qentry(base, prob));
		}
	
	vector<qentry> sortedNonZeros;
	while (!que.empty())
	{
		sortedNonZeros.push_back(que.top());
		que.pop();
	}
	return sortedNonZeros;
}

qubase measure(Q)
{
	float prob = rand_float();
	if (q.dense)
		DENSE_ITER(base)
		{
			prob -= norm(q.amp[base]);
			if (prob <= 0)
				return base;
		}
	else
		for (qubase& base : q.base_iter())
		{
			prob -= norm(q[base]);
			if (prob <= 0)
				return base;
		}

	// should never be here
	throw QuantumException("Measurement fails, probability doesn't sum up to 1.");
}

int measure(Q, int tar)
{
	int result = 0;
	float prob = 0; // probability of being 1
	qubase t = q.to_qubase(tar);

	// get probability of this target bit collapsing to 0 or 1
	if (q.dense)
	{
		auto& amp = q.amp;
		DENSE_ITER(base)
			if (base & t) // target bit 1
				prob += norm(amp[base]);
		if (prob > rand_float())
			result = 1;
		float newNorm = result == 1 ? sqrt(prob) : sqrt(1 - prob);
		// eliminate all states that don't agree
		DENSE_ITER(base)
			if (((base & t) != 0) == result)
				amp[base] /= newNorm;
			else
				amp[base] = CX(0);
	}
	else // sparse
	{
		for (qubase& base : q.base_iter())
			if (base & t) // target bit 1
				prob += norm(q[base]);
		if (prob > rand_float())
			result = 1;
		float newNorm = result == 1 ? sqrt(prob) : sqrt(1 - prob);

		vector<CX> newAmp;
		newAmp.reserve(q.amp.capacity() / 2);
		vector<qubase> newBasis;
		newBasis.reserve(q.basis.capacity() / 2);
		size_t s = 0; // new size
		for (qubase& base : q.base_iter())
			if (((base & t) != 0) == result)
			{
				newAmp.push_back(q[base] / newNorm);
				newBasis.push_back(base);
				q.basemap[base] = s++;
			}
			else // eliminate that base
				q.basemap.erase(base);
		q.amp = move(newAmp);
		q.basis = move(newBasis);
	}
	return result;
}

uint64_t measure_top(Q, int topQubits, bool destructive)
{
	if (destructive)
	{
		// partial measurement: discard the last output bits
		uint64_t result = 0;
		for (int qi = 0; qi < topQubits; ++qi)
			// going from msb to lsb
			result |= measure(q, qi) << (topQubits - 1 - qi);
		return result;
	}
	else
		return measure(q) >> (q.nqubit - topQubits);
}


// Apply to n most significant bits
void apply_oracle(Q, const oracle_function& oracle, int inputQubits)
{
	if (inputQubits >= q.nqubit)
		throw QuantumException("inputQubits should not exceed total qubits");
	int outputQubits = q.nqubit - inputQubits;
	auto& amp = q.amp;
	// needs a temporary amp that holds outputQubits' original amplitude
	uint64_t outputSize = 1 << outputQubits;
	if (q.dense)
	{
		vector<CX> tmpAmp(outputSize);
		for (uint64_t input = 0; input < 1 << inputQubits; ++input)
		{
			uint64_t ans = oracle(input);
			std::copy_n(&amp[input << outputQubits], outputSize, &tmpAmp[0]);
			for (uint64_t output = 0; output < outputSize; ++output)
			{
				qubase newBase = (input << outputQubits) | (output ^ ans);
				amp[newBase] = tmpAmp[output];
			}
		}
	}
	else
	{
		uint64_t outputMask = outputSize - 1; // all ones in the lower bits
		vector<CX> newAmp;
		newAmp.reserve(q.amp.capacity());
		vector<qubase> newBasis;
		newBasis.reserve(q.basis.capacity());
		decltype(q.basemap) newBasemap(q.basemap.size());
		size_t s = 0; // new size

		for (qubase& base : q.base_iter())
		{
			uint64_t input = base >> outputQubits; // most sig bits
			uint64_t output = base & outputMask;
			uint64_t ans = oracle(input);
			if (norm(q[base]) < TOL)
				continue;  // purge along the way
			newAmp.push_back(q[base]);
			qubase newBase = input << outputQubits | (output ^ ans);
			newBasis.push_back(newBase);
			newBasemap[newBase] = s++;
		}
		q.amp = move(newAmp);
		q.basis = move(newBasis);
		q.basemap = move(newBasemap);
	}
}