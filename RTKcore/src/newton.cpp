//////////////////////////////////////////////////////////////////////////////
//  ニュートン法によりekを算出
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <math.h>
#include "global_extern.h"

using namespace std;

double newton(double mk, double e, int i)
{
	double ek;
	double a1,a2,a3;
	double tol = 1.0e-15;
	int iter, lmax = 4;
	i=i;

	ek = 0.0;
	if(fabs(mk) > 0.0)
	{
		a1 = mk + e*sin(mk);
		iter = 1;
		while(iter <= lmax)
		{
			a2 = a1 - e*sin(a1) - mk;
			a3 = a1 - a2/2.0;
			a3 = 1.0 - e*cos(a3);
			ek = a1 - a2/a3;
			a2 = a1 - ek;
			if(fabs(a2) < tol){break ;}
			++iter;
			if(iter <= lmax)
			{
				a1 = ek;
			}
			else{
				i=i;
//				fprintf(fp[3],"Newton Error   PRN=%d\n",i);
				//cout << "newton_error!!   " << i << endl;
			}
		}
	}
	return(ek);
}
