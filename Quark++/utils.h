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

string int2str(int a);
/*
 *	Convert to a bit string
 * May also specify minimal number of bits to print (template parameter)
 */
template<int minbit>
string bits2str(qubase b)
{
	return bitset<minbit>(b).to_string();
}
template <>
string bits2str<0>(qubase b);
string bits2str(qubase b);

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

	template<typename T>
	inline void pvec(vector<T> vec) { pr(vec2str(vec)); }

	inline void pause() { cin.get(); }
}

#endif // utils_h__