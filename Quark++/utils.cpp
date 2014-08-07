#include "utils.h"

// Convert to bit string
template<>
string bits2str<0>(qubase b)
{
	string s = bitset<32>(b).to_string();

	// kill leading zeros
	int i = 0;
	while (i < s.size() && s[i] == '0') { ++i; }
	return i != s.size() ? s.substr(i) : "0";
}
string bits2str(qubase b) { return bits2str<0>(b); }

string int2str(int a)
{
	ostringstream oss;
	oss << a;
	return oss.str();
}