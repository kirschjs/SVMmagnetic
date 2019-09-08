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
	double h2m, L, apot3b, vpot3b, eB, eBspin, momega;
	int npar, npt, nop, dmax, nPairs, nPerm, nConf,ibf;
	
	std::vector<VectorXd> bbox; //vector of cells
	std::vector<MatrixXd> PM;
	std::vector<MatrixXd> NonZeroPM;
	std::vector<VectorXi> PV;
	std::vector<int> parity;
	std::vector<VectorXd> TBP;
		
	std::vector<double> stme;
    std::vector<double> magnetic_me;
    std::vector<double> magnetic_spin_me;
    std::vector<double> magnetic_charge_me;
	
	MatrixXd vpot, apot;
    MatrixXd mass;

	void PreparePotential(Input &input);
	void PrepareSpinIsospinME(Input &input);
    
    void PrepareMagneticSpinME(Input &input);
    void PrepareMagneticME(Input &input);
    double One_magnetic_me(int ipar, VectorXi Perm, int i_spin_coupling);

    std::vector<VectorXd> magnetic_me_pair;
    void PrepareMagneticMePair(Input &input);
    double One_magnetic_me_pair(int i, int j, VectorXi Perm, int index);

    int NchargPar;
    int CountNchargPar();

    IsospinState isospin;
    SpinState spin;
//==============
    double KinEnergy_cof, MagneticEnergy_cof, harmonic_cof, MagneticSpin_cof;
//==============

public:
	MatrixElement(Input &input);
	~MatrixElement();
	double overlap(std::vector<MatrixXd> state1, std::vector<MatrixXd> state2);
	double energy(std::vector<MatrixXd> state1, std::vector<MatrixXd> state2);
};
#endif 
