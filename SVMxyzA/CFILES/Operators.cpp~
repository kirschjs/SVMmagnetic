#include "Operators.h"
#include "Input.h"
#include <vector>
#include <Eigen> 
using namespace Eigen;
using namespace std;

Operators::Operators(Input &input)
{
	npar    = input.npar;
	spin    = input.spin;
	isospin = input.isospin;

}
//=============================================================================
/* 
    calculate the matrix elements in spin-isospin space 
         < S T | Oij(kind) Perm | S T >
    | S T > - vector state defined in the constructor
    O(kind) - operator
    i,j     - interacting particles
*/
double Operators::O(int i, int j, VectorXi Perm, int kind)
{
	double x = 0;
	for (int icmp = 0; icmp < spin.ncmp; icmp++){
	  for (int jcmp = 0; jcmp < spin.ncmp; jcmp++){
	    x = x + spin.coef[icmp] * spin.coef[jcmp] 
	          * SpinOp(spin.cmp[icmp].sz,spin.cmp[jcmp].sz ,i ,j ,Perm ,kind);
	  }
	}

	double y = 0;
	for (int icmp = 0; icmp < isospin.ncmp; icmp++){
	  for (int jcmp = 0; jcmp < isospin.ncmp; jcmp++){
	    y = y + isospin.coef[icmp] * isospin.coef[jcmp] 
                  * IsospinOp(isospin.cmp[icmp].tz, isospin.cmp[jcmp].tz, i, j, Perm, kind);
	  }
	}

//	cout << "\t iop=" << kind << "\t me t = " << y << "\t me s = " << x << endl;
	return x*y;
}


double Operators::SpinOp(std::vector<int> sz1, std::vector<int> sz2, int i, int j, VectorXi Perm, int kind)
{

	for (int ipar = 0; ipar < npar; ipar++)
	{
		if ((ipar != i) && (ipar != j))
		{
			if (sz1[ipar] != sz2[Perm(ipar)]) return 0;
		}
	}

	/* iop: 0-Wigner, 1-Psigma*Ptau, 2-Psigma, 3-Ptau */
	if (kind == 0) 	return O_0(sz1[i], sz1[j], sz2[Perm(i)], sz2[Perm(j)]);
	if (kind == 1)	return O_1(sz1[i], sz1[j], sz2[Perm(i)], sz2[Perm(j)]);
	if (kind == 2)	return O_1(sz1[i], sz1[j], sz2[Perm(i)], sz2[Perm(j)]);
	if (kind == 3)	return O_0(sz1[i], sz1[j], sz2[Perm(i)], sz2[Perm(j)]);
}


double Operators::IsospinOp(std::vector<int> tz1, std::vector<int> tz2, int i, int j, VectorXi Perm, int kind)
{

	for (int ipar = 0; ipar < npar; ipar++)
	{
		if ((ipar != i) && (ipar != j))
		{
			if (tz1[ipar] != tz2[Perm(ipar)]) return 0;
		}
	}

	/* iop: 0-Wigner, 1-Psigma*Ptau, 2-Psigma, 3-Ptau */
	if (kind == 0) 	return O_0(tz1[i], tz1[j], tz2[Perm(i)], tz2[Perm(j)]);
	if (kind == 1)	return O_1(tz1[i], tz1[j], tz2[Perm(i)], tz2[Perm(j)]);
	if (kind == 2)	return O_0(tz1[i], tz1[j], tz2[Perm(i)], tz2[Perm(j)]);
	if (kind == 3)	return O_1(tz1[i], tz1[j], tz2[Perm(i)], tz2[Perm(j)]);
}


double Operators::O_0(int i1, int i2, int i3, int i4) //identity operator
{
	//	1 for spin up
	//	2 for spin down
	double x = 0;
	if ((i1 == i3) && (i2 == i4)) x = 1.0;
	return x;
}
double Operators::O_1(int i1, int i2, int i3, int i4) //permutation operator
{
	//	1 for spin up
	//	2 for spin down
	double x = 0;
	if ((i1 == i4) && (i2 == i3)) x = 1;
	return x;
}

Operators::~Operators()
{
}
