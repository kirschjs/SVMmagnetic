
#ifndef SVM_H
#define SVM_H
#include "MatrixElement.h"
#include "Rand.h"
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
	int CheckOverlap(vector<vector<MatrixXd>> Basis);
	vector<MatrixXd> FirstNewState();
	double NewEnergy(vector<vector<MatrixXd>> Basis, MatrixXd C, VectorXd D, double E, double EE);
	vector<MatrixXd> NewState(vector<vector<MatrixXd>> Basis, MatrixXd C, VectorXd D, double E, double EE);
	MatrixXd Dmatrix();
	MatrixXd A(MatrixXd d);
	MatrixXd CalcNormMatrix(vector<vector<MatrixXd>> Basis);
	MatrixXd CalcHamiltonianMatrix(vector<vector<MatrixXd>> Basis);
	MatrixXd NormMatrix(vector<vector<MatrixXd>> Basis);
	MatrixXd HamiltonianMatrix(vector<vector<MatrixXd>> Basis);
	MatrixXd Hmatrix,Nmatrix;
	void UpdateNorm(vector<vector<MatrixXd>> Basis);
	void UpdateHamiltonian(vector<vector<MatrixXd>> Basis);
};


#endif 







