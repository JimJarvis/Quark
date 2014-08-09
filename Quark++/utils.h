#ifndef utils_h__
#define utils_h__

#include <iostream>
#include <string>
#include <complex>
#include <cmath>
#include <vector>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <algorithm>
#include <memory>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <bitset>
#include <unordered_map>
#include <Eigen/Dense>

using namespace std;
using Eigen::Matrix2cf;
using Eigen::Matrix4cf;
using Eigen::VectorXcf;
using Eigen::MatrixXcf;

typedef float REAL;
typedef complex<REAL> CX;
typedef unsigned long long qubase;

// Amplitude tolerance: smaller than this will be considered 0
#define TOL 1e-7

inline string int2str(int a)
{
	ostringstream oss;
	oss << a;
	return oss.str();
}

template<typename T>
string vec2str(vector<T> vec)
{
	ostringstream oss;
	oss << "[";
	for (T& ele : vec)
		oss << ele << ", ";
	string s = oss.str();
	return s.substr(0, s.size() - 2) + "]";
}

//******** Bit ops
/*
 *	Convert to a bit string
 * May also specify minimal number of bits to print (template parameter)
 */
template<int minbit>
string bits2str(uint64_t b)
{
	return bitset<minbit>(b).to_string();
}

// Convert to bit string
template<>
inline string bits2str<0>(qubase b)
{
	string s = bitset<32>(b).to_string();

	// kill leading zeros
	int i = 0;
	while (i < s.size() && s[i] == '0') { ++i; }
	return i != s.size() ? s.substr(i) : "0";
}
inline string bits2str(qubase b) { return bits2str<0>(b); }

template <typename T, T m, int k>
static inline T swapbits(T p)
{
	T q = ((p >> k) ^ p)&m;
	return p^q ^ (q << k);
}
inline uint64_t bit_reverse(uint64_t n)
{
	static const uint64_t m0 = 0x5555555555555555ULL;
	static const uint64_t m1 = 0x0300c0303030c303ULL;
	static const uint64_t m2 = 0x00c0300c03f0003fULL;
	static const uint64_t m3 = 0x00000ffc00003fffULL;
	n = ((n >> 1)&m0) | (n&m0) << 1;
	n = swapbits<uint64_t, m1, 4>(n);
	n = swapbits<uint64_t, m2, 8>(n);
	n = swapbits<uint64_t, m3, 20>(n);
	n = (n >> 34) | (n << 30);
	return n;
}

///////************** For-range loop iterables **************///////
/*
 * inclusive/exclusive: 
 *	Forward: [begin, end)
 * Reverse: (begin, end]
 * In reverse mode, if begin is 0, traverse from the last element
 */
template<typename T, bool forward = true>
struct VecRange
{
	vector<T>& vec;
	size_t beginIdx, endIdx;

	// customized iterator
	struct iter
	{
		size_t i;
		vector<T>& vec;
		iter(vector<T>& _vec, size_t _i) : vec(_vec), i(_i) {}

		T& operator *() { return vec[i]; }

		// prefix
		size_t& operator++() { return forward ? ++i : --i; }

		bool operator==(const iter& other) const
		{
			return this->i == other.i;
		}
		bool operator!=(const iter& other) const
		{
			return this->i != other.i;
		}
	};

	VecRange(vector<T>& _vec, size_t begin, size_t end) :
		vec(_vec), 
		beginIdx(forward ? begin : begin - 1), 
		endIdx(forward ? end : end -1) { }
	VecRange(vector<T>& _vec, size_t begin = 0) :
		vec(_vec),
		beginIdx(forward ? begin : 
				(begin == 0 ? vec.size() : begin) - 1),
				endIdx(forward ? vec.size() : -1) { }

	iter begin() { return iter(vec, beginIdx); }
	iter end() { return iter(vec, endIdx); }
};

/*
 * inclusive/exclusive: 
 *	Forward: [begin, end)
 * Reverse: (begin, end]
 */
template<typename IntType = int, bool forward = true>
struct Range
{
	IntType beginIdx, endIdx;

	// customized iterator
	struct iter
	{
		IntType i;
		iter(IntType _i) : i(_i) {}

		IntType& operator *() { return i; }
		// prefix
		IntType& operator++() { return forward ? ++i : --i; }
		bool operator==(const iter& other) const
		{
			return this->i == other.i;
		}
		bool operator!=(const iter& other) const
		{
			return this->i != other.i;
		}
	};

	Range(IntType begin, IntType end) :
		beginIdx(forward ? begin : begin - 1), 
		endIdx(forward ? end : end -1) { }
	Range(IntType end) :
		beginIdx(forward ? 0 : end - 1), 
		endIdx(forward ? end : -1) { }

	iter begin() { return iter(beginIdx); }
	iter end() { return iter(endIdx); }
};

///////************** Exceptions **************///////
// thrown when requested file can't be opened
class QuantumException : public exception
{
public:
	QuantumException(string fp) : errmsg("Quantum error: ")
	{
		errmsg += fp;
	}
#ifdef _MSC_VER  // VC++ doesn't yet support noexcept()
	virtual const char* what() const throw()
	{
		return errmsg.c_str();
	}
#else  // C++11 noexcept operator. Required by gcc
	virtual const char* what() noexcept(true)
	{
		return errmsg.c_str();
	}
#endif
private:
	string errmsg;
};


///////************** Debugging **************///////
namespace Testing
{
	inline void ptitle(string title = "") 
	{ cout << "！！！！！！！！！！" << title << "！！！！！！！！！！！" << endl; }
#define pr(X) cout << X << endl
#define pause std::cin.get()

	template<typename T>
	inline void pvec(vector<T> vec) { pr(vec2str(vec)); }
}

#endif // utils_h__