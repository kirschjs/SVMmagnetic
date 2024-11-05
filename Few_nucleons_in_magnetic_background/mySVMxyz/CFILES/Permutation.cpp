#include "Permutation.h"
#include <cmath>
#include <Eigen> 
#include <vector>
#include <iostream>
using namespace Eigen;
using namespace std;


Permutation::Permutation(int N) : npar(N)
{
	VectorXi nperms_i = VectorXi::Zero(npar);
	VectorXi perm = VectorXi::Zero(npar);

	int nperms = 1;
	for (int i = 1; i < npar + 1; i++) { nperms = nperms*i; }
	permutation.reserve(nperms);
	perm_matrix.reserve(nperms);
	parity.reserve(nperms);

	/* create permutation list */
	for (int i = 0; i < npar; i++) { perm(i) = i; }
	int np = 1;
	permutation.push_back(perm);
	parity.push_back(1);
	nperms_i(0) = 1;
	for (int ipar = 1; ipar < npar; ipar++) {
	  for (int iperm = 0; iperm < nperms_i(ipar - 1); iperm++) {
	    VectorXi refperm = permutation[iperm];
	    int      refparity = parity[iperm];
	    for (int jpar = 0; jpar < ipar; jpar++) {
	      np = np + 1;
	      perm = refperm;
	      perm(jpar) = ipar; perm(ipar) = refperm(jpar);
	      permutation.push_back(perm); parity.push_back(-1 * refparity);
	    }
	  }
	  nperms_i[ipar] = np;
	}
	if (np != nperms){cout << "Permutations: funny error" << endl; exit(0);}

	/* create paermutation matrices */
	for (int ip = 0; ip < nperms; ip++) {
		MatrixXd pmat = MatrixXd::Zero(npar, npar);
		perm = permutation[ip];
		for (int i = 0; i<npar; i++) { pmat(i, perm(i)) = 1.0; }
		perm_matrix.push_back(pmat.transpose());
	}


}

Permutation::~Permutation(){}
/*
//=============================================================================
void Permutation::printperm(int iperm)
{
  VectorXi &perm = permutation[iperm];
  fprintf(printfile, "\t permutation %5d   parity = %2d   ["
	  ,iperm,parity[iperm]);
  for (int i=0; i < npar; i++)
    {fprintf(printfile,"%2d",perm(i));}
  fprintf(printfile, "]\n");
}
//=============================================================================
void Permutation::printpmat(int iperm)
{
  MatrixXd &pmat = perm_matrix[iperm];
  std::string cperm ;
  for (int i=0; i < npar; i++){
    fprintf(printfile,"\t\t\t\t\t") ;
    for (int j=0; j < npar; j++){
      double x = pmat(i,j);
      fprintf(printfile,"%6.1f ",x) ;
    }
    fprintf(printfile,"\n") ;
  }
}
//=============================================================================
*/



