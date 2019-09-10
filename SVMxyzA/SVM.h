#ifndef SVM_H
#define SVM_H
#include "MatrixElement.h"
#include "Rand.h"
#include "BasisState.h"
#include "Input.h"
#include <Eigen> 
#include <vector>
using namespace Eigen;
using namespace std;

class SVM
{
private:
	Rand &rr;
	MatrixElement me;
	int N, seed, mm0, kk0, mnb;
	double bmin, bmax;
        int iBoxInf;

public:
  SVM(Rand &r, Input &input);
  ~SVM();
  int CheckOverlap(vector<BasisState> &Basis);
  BasisState FirstNewState();
  double NewEnergy(vector<BasisState> &Basis, MatrixXd C, VectorXd D, double E, double EE);
  BasisState NewState(vector<BasisState> &Basis, MatrixXd C, VectorXd D, double E, double EE);
  MatrixXd Dmatrix();
  MatrixXd A(MatrixXd d);
  MatrixXd CalcNormMatrix(vector<BasisState> &Basis);
  MatrixXd CalcHamiltonianMatrix(vector<BasisState> &Basis);
  MatrixXd NormMatrix(vector<BasisState> &Basis);
  MatrixXd HamiltonianMatrix(vector<BasisState> &Basis);
  MatrixXd Hmatrix,Nmatrix;
  void UpdateNorm(vector<BasisState> &Basis);
  void UpdateHamiltonian(vector<BasisState> &Basis);
};


#endif 







