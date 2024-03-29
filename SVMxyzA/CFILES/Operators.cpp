#include "Operators.h"
#include "Input.h"
#include "BasisState.h"
#include "MatrixElement.h"
#include <vector>
#include <Eigen> 
using namespace Eigen;
using namespace std;

Operators::Operators(Input &input)
{
	npar    = input.npar;
	spin    = input.spin;
	isospin = input.isospin;
	nop2b   = input.nop2b;
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

//=============================================================================
/* 
    calculate the matrix elements in spin-isospin space 
         < S T | Oij(kind) Perm | S T >
    | S T > - vector state defined in the constructor
    O(iop) - operator
    i,j     - interacting particles
*/

vector<ME2bST>  Operators::OST_2bme(SpinIsospinState &state1, SpinIsospinState &state2, VectorXi &Perm){
    std::vector<ME2bST> op2b(nop2b);

    for (int iop2b = 0; iop2b < nop2b; iop2b++){
	op2b[iop2b].me.resize(npar,npar);
	op2b[iop2b].me = MatrixXd::Zero(npar,npar);
    }
    
    bool     keypr= false;
    double   spin_me,isospin_me;
    VectorXi sz1(npar),sz2(npar),tz1(npar),tz2(npar);
    for (int icmp = 0; icmp < state1.ncmp; icmp++){
	for (int jcmp = 0; jcmp < state2.ncmp; jcmp++){
	    sz1=state1.sz.col(icmp); sz2=state2.sz.col(jcmp);
	    tz1=state1.tz.col(icmp); tz2=state2.tz.col(jcmp);
	    if (keypr) printf("OST_2bme 1: %4d C = %9.4f |%3d %3d >\n",icmp,state1.coef(icmp),sz1[0],sz1[1]);
	    if (keypr) printf("OST_2bme 2: %4d C = %9.4f |%3d %3d >\n",jcmp,state2.coef(jcmp),sz2[0],sz2[1]);
	    double fctr = state1.coef(icmp) * state2.coef(jcmp);
	    for (int iop2b = 0; iop2b < nop2b; iop2b++){       
		for (int ip = 0; ip < npar; ip++){
		    for (int jp = ip + 1; jp < npar; jp++){
			spin_me    = Spin2bME(   sz1,sz2 ,ip ,jp ,Perm ,iop2b);
			isospin_me = Isospin2bME(tz1,tz2 ,ip ,jp ,Perm ,iop2b);
			op2b[iop2b].me(ip,jp) = op2b[iop2b].me(ip,jp) +
			    fctr*spin_me*isospin_me;
			if (keypr) 
			    printf("\t  iop= %3d  fctr= %9.4f  me= %9.4f\n",iop2b,fctr,op2b[iop2b].me(ip,jp));
		    }
		}
	    }
	}
    }

    if (keypr) {
	for (int iop2b = 0; iop2b < nop2b; iop2b++){       
	    for (int ip = 0; ip < npar; ip++){
		for (int jp = ip + 1; jp < npar; jp++){
		    spin_me    = Spin2bME(   sz1,sz2 ,ip ,jp ,Perm ,iop2b);
		    isospin_me = Isospin2bME(tz1,tz2 ,ip ,jp ,Perm ,iop2b);
		    printf("\t  OST_2bme: iop= %3d    me= %9.4f\n",iop2b,op2b[iop2b].me(ip,jp));
		}
	    }
	}
    }

    return op2b;

}

//=============================================================================
double Operators::Spin2bME(VectorXi &sz1, VectorXi &sz2, int i, int j, VectorXi &Perm, int iop){

    for (int ipar = 0; ipar < npar; ipar++){
	if ((ipar != i) && (ipar != j)){
	    if (sz1(ipar) != sz2(Perm(ipar))) return 0;
	}
    }

    /* iop: 0-Wigner, 1-Psigma*Ptau, 2-Psigma, 3-Ptau */
         if (iop == 0) return O_0(sz1(i), sz1(j), sz2(Perm(i)), sz2(Perm(j)) );
    else if (iop == 1) return O_1(sz1(i), sz1(j), sz2(Perm(i)), sz2(Perm(j)) );
    else if (iop == 2) return O_1(sz1(i), sz1(j), sz2(Perm(i)), sz2(Perm(j)) );
    else if (iop == 3) return O_0(sz1(i), sz1(j), sz2(Perm(i)), sz2(Perm(j)) );
    else {cout << "\n\tERROR undefined iop (Spin2bME)"<<endl;exit(0);}
}
//=============================================================================
double Operators::Isospin2bME(VectorXi &tz1, VectorXi &tz2, int i, int j, VectorXi &Perm, int iop){

    for (int ipar = 0; ipar < npar; ipar++){
	if ((ipar != i) && (ipar != j)){
	    if (tz1(ipar) != tz2(Perm(ipar))) return 0;
	}
    }

    /* iop: 0-Wigner, 1-Psigma*Ptau, 2-Psigma, 3-Ptau */
         if (iop == 0) return O_0(tz1(i), tz1(j), tz2(Perm(i)), tz2(Perm(j)) );
    else if (iop == 1) return O_1(tz1(i), tz1(j), tz2(Perm(i)), tz2(Perm(j)) );
    else if (iop == 2) return O_0(tz1(i), tz1(j), tz2(Perm(i)), tz2(Perm(j)) );
    else if (iop == 3) return O_1(tz1(i), tz1(j), tz2(Perm(i)), tz2(Perm(j)) );
    else {cout << "\n\tERROR undefined iop (Isospin2bME)"<<endl;exit(0);}
}
//=============================================================================
double Operators::SpinOp(std::vector<int> sz1, std::vector<int> sz2, int i, int j, VectorXi Perm, int kind)
{

	for (int ipar = 0; ipar < npar; ipar++)	{
	    if ((ipar != i) && (ipar != j)){
		if (sz1[ipar] != sz2[Perm(ipar)]) return 0;
	    }
	}

	/* iop: 0-Wigner, 1-Psigma*Ptau, 2-Psigma, 3-Ptau */
	if (kind == 0) 	return O_0(sz1[i], sz1[j], sz2[Perm(i)], sz2[Perm(j)]);
	if (kind == 1)	return O_1(sz1[i], sz1[j], sz2[Perm(i)], sz2[Perm(j)]);
	if (kind == 2)	return O_1(sz1[i], sz1[j], sz2[Perm(i)], sz2[Perm(j)]);
	if (kind == 3)	return O_0(sz1[i], sz1[j], sz2[Perm(i)], sz2[Perm(j)]);
}

//=============================================================================
double Operators::IsospinOp(std::vector<int> tz1, std::vector<int> tz2, int i, int j, VectorXi Perm, int kind)
{

	for (int ipar = 0; ipar < npar; ipar++)	{
	    if ((ipar != i) && (ipar != j)){
		if (tz1[ipar] != tz2[Perm(ipar)]) return 0;
	    }
	}

	/* iop: 0-Wigner, 1-Psigma*Ptau, 2-Psigma, 3-Ptau */
	if (kind == 0) 	return O_0(tz1[i], tz1[j], tz2[Perm(i)], tz2[Perm(j)]);
	if (kind == 1)	return O_1(tz1[i], tz1[j], tz2[Perm(i)], tz2[Perm(j)]);
	if (kind == 2)	return O_0(tz1[i], tz1[j], tz2[Perm(i)], tz2[Perm(j)]);
	if (kind == 3)	return O_1(tz1[i], tz1[j], tz2[Perm(i)], tz2[Perm(j)]);
}

//=============================================================================
double Operators::O_0(int i1, int i2, int i3, int i4) //identity operator
{
	//	1 for spin up
	//	2 for spin down
	double x = 0;
	if ((i1 == i3) && (i2 == i4)) x = 1.0;
	return x;
}
//=============================================================================
double Operators::O_1(int i1, int i2, int i3, int i4) //permutation operator
{
	//	1 for spin up
	//	2 for spin down
	double x = 0;
	if ((i1 == i4) && (i2 == i3)) x = 1;
	return x;
}
//=============================================================================
Operators::~Operators(){}
//=============================================================================
