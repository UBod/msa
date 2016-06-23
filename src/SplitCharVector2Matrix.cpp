#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

#include "SplitCharVector2Matrix.h"

RcppExport SEXP SplitCharVector2Matrix(SEXP xR, SEXP replR)
{
    vector<string> x = as< vector<string> >(xR);
    string repl = as< string >(replR);
    int i, j, n = x.size(), m = x[0].length();
    CharacterMatrix out(n, m);

    for (i = 0; i < n; i++)
    {
	m = x[i].length();
	
	for (j = 0; j < m; j++)
	{
	    string tmp = x[i].substr(j, 1);

	    if (tmp.compare("-") == 0)
		tmp = repl;

	    out(i, j) = tmp;
	}
    }

    return(out);
}
