#ifndef MatrixElement_H
#define MatrixElement_H
#include "Input.h"
#include <vector>
#include <Eigen> 
using namespace Eigen;

class MatrixElement
{
private:
	const double pi = 3.1415926;
	double h2m, L, apot3b, vpot3b, eB, momega;
	std::vector<VectorXd> bbox; //vector of cells
	std::vector<MatrixXd> PM;
	std::vector<MatrixXd> NonZeroPM;
	std::vector<VectorXi> PV;
	std::vector<int> parity;
	std::vector<VectorXd> TBP;
	int npar, npt, nop, dmax, nPairs, nPerm, nConf,ibf;
	std::vector<double> stme;
	MatrixXd vpot, apot;
        MatrixXd mass;
	void PreparePotential(Input &input);
	void PrepareSpinIsospinME(Input &input);
        int iBoxInf;

public:
	MatrixElement(Input &input);
	~MatrixElement();
	double overlap(std::vector<MatrixXd> state1, std::vector<MatrixXd> state2);
	double energy(std::vector<MatrixXd> state1, std::vector<MatrixXd> state2);
};
#endif 
