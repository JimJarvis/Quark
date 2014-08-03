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

using namespace std;

typedef float REAL;
typedef complex<REAL> CX;
typedef unsigned long long qubase;

template<typename T>
extern string vec2str(vector<T> vec);
extern string int2str(int a);
/*
 *	Convert to a bit string
 * May also specify minimal number of bits to print (template parameter)
 */
template<int minbit>
extern string bits2str(qubase b)
{
	return bitset<minbit>(b).to_string();
}
template <>
extern string bits2str<0>(qubase b);
extern string bits2str(qubase b);

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