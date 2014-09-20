#ifndef frac_h__
#define frac_h__

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <cmath>
using namespace std;

inline int64_t gcd(int64_t a, int64_t b)
{
	int64_t c;
	while (a != 0)
	{
		c = a;
		a = b % a;
		b = c;
	}
	return b;
}

class Frac
{
public:
	int64_t num, denom;

	Frac()
	{
		num = 0;
		denom = 1;
	}

	Frac(int64_t n)
	{
		num = n;
		denom = 1;
	}

	Frac(int64_t n, int64_t d)
	{
		if (d == 0)
		{
			cerr << "Denominator may not be 0." << endl;
			exit(0);
		}
		else if (n == 0)
		{
			num = 0;
			denom = 1;
		}
		else
		{
			int sign = 1;
			if (n < 0)
			{
				sign *= -1;
				n *= -1;
			}
			if (d < 0)
			{
				sign *= -1;
				d *= -1;
			}

			long long tmp = gcd(n, d);
			num = n / tmp*sign;
			denom = d / tmp;
		}
	}

	operator int() { return (num) / denom; }
	operator float() { return ((float)num) / denom; }
	operator double() { return ((double)num) / denom; }
};

Frac operator+(const Frac& lhs, const Frac& rhs)
{
	return Frac(lhs.num*rhs.denom
				 + rhs.num*lhs.denom,
				 lhs.denom*rhs.denom);
}

Frac operator+=(Frac& lhs, const Frac& rhs)
{
	Frac tmp(lhs.num*rhs.denom
				 + rhs.num*lhs.denom,
				 lhs.denom*rhs.denom);
	lhs = tmp;
	return lhs;
}

Frac operator-(const Frac& lhs, const Frac& rhs)
{
	return Frac(lhs.num*rhs.denom
				 - rhs.num*lhs.denom,
				 lhs.denom*rhs.denom);
}

Frac operator-=(Frac& lhs, const Frac& rhs)
{
	Frac tmp(lhs.num*rhs.denom
				 - rhs.num*lhs.denom,
				 lhs.denom*rhs.denom);
	lhs = tmp;
	return lhs;
}

Frac operator*(const Frac& lhs, const Frac& rhs)
{
	return Frac(lhs.num*rhs.num,
				 lhs.denom*rhs.denom);
}

Frac operator*=(Frac& lhs, const Frac& rhs)
{
	Frac tmp(lhs.num*rhs.num,
				 lhs.denom*rhs.denom);
	lhs = tmp;
	return lhs;
}

Frac operator*(int lhs, const Frac& rhs)
{
	return Frac(lhs*rhs.num, rhs.denom);
}

Frac operator*(const Frac& rhs, int lhs)
{
	return Frac(lhs*rhs.num, rhs.denom);
}

Frac operator/(const Frac& lhs, const Frac& rhs)
{
	return Frac(lhs.num*rhs.denom,
				 lhs.denom*rhs.num);
}

// reciprocal
Frac operator~(const Frac& f)
{
	return Frac(f.denom, f.num);
}

std::ostream& operator<<(std::ostream &strm, const Frac &a)
{
	if (a.denom == 1)
		strm << a.num;
	else
		strm << a.num << "/" << a.denom;
	return strm;
}

///////************** Continued fraction **************///////
typedef vector<int> ContFrac;

inline Frac to_frac(ContFrac cfrac, int size = 0)
{
	if (size < 1)
		size = cfrac.size();
	Frac ans(1, cfrac[size - 1]);
	for (int i = size - 2; i >= 1; --i)
		ans = ~(ans + Frac(cfrac[i]));
	return ans + Frac(cfrac[0]);
}

#endif // frac_h__
