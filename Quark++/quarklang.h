/**********************************************/
/*********** Specific for quarklang compiler ***********/
/**********************************************/
#ifndef quarklang_h__
#define quarklang_h__
#include "utils.h"

template <typename T>
vector<T> concat_vector(vector<T> vec1, vector<T> vec2)
{
	vector<T> ans = vec1;
	ans.insert(ans.end(), vec2.begin(), vec2.end());
	return ans;
}


#endif // quarklang_h__

