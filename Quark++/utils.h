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
#include <Eigen/Dense>

using namespace std;
using Eigen::Matrix2cf;
using Eigen::Matrix4cf;

typedef float REAL;
typedef complex<REAL> CX;
typedef unsigned long long qubase;

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

///////***** Testing *****///////
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