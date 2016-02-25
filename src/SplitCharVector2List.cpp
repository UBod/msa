#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

#include "SplitCharVector2List.h"

RcppExport SEXP SplitCharVector2List(SEXP xR)
{
    vector<string> x = as< vector<string> >(xR);
    int i, j, n = x.size();
    List out;

    for (i = 0; i < n; i++)
    {
	int len = x[i].length();
	vector<string> tmp;
	
	for (j = 0; j < len; j++)
	    tmp.push_back(x[i].substr(j, 1));

	out.push_back(tmp);
    }

    return(out);
}
