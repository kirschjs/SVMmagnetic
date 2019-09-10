#ifndef MatrixElement_H
#define MatrixElement_H
#include "Input.h"
#include "BasisState.h"
#include <vector>
#include <Eigen> 
using namespace Eigen;

//=============================================================================
/*
    actp    - active particles
    op      - spin-isospin operator
    me      - actual matrix elements 

    iperm   - index of permutation
    ipair   - index of pair
    ipar    - index of particle
    iop     - index of ST operator
    ijst    - index of left/right spin-isospin A-body basis states.  
              ijst = ist*nst+jst.   
 
    me1b = stme1b[ijst].perm[iperm].op[iop].me(ipar)
    me2b = stme2b[ijst].perm[iperm].op[iop].me(ipar,jpar)
*/
struct ME1bST   {VectorXd me;};
struct ME2bST   {VectorXd me;};
struct MEofSTOp {std::vector<ME1bST>   op1b; std::vector<ME2bST>   op2b;};
struct MESTperm {std::vector<MEofSTOp> perm;};
//=============================================================================

class MatrixElement
{
private:
    const double pi = 4.0*atan(1.0);

    double h2m, L, apot3b, vpot3b, eB, eBspin, momega;
    int npar, npt, nop, dmax, nPairs, nPerm, nConf,ibf;
    int nop2b,nop1b;
    
    std::vector<VectorXd> bbox; //vector of cells
    std::vector<MatrixXd> PM;
    std::vector<MatrixXd> NonZeroPM;
    std::vector<VectorXi> PV;
    std::vector<int> parity;
    std::vector<VectorXd> TBP;
    
//    std::vector<double> stme;
    std::vector<double> magnetic_me;
    std::vector<double> magnetic_spin_me;
    std::vector<double> magnetic_charge_me;
	
    std::vector<MESTperm> stmeop;

    MatrixXd vpot, apot;
    MatrixXd mass;

    void PreparePotential(Input &input);
    void PrepareSpinIsospinME(Input &input);
    
    void PrepareMagneticSpinME(Input &input);
    void PrepareMagneticME(Input &input);
    double One_magnetic_me(int ipar, VectorXi Perm, int i_spin_coupling);
    double SpinOp(std::vector<int> sz1, std::vector<int> sz2, int npar, int ipar, 
		  VectorXi Perm, int i_spin_coupling);
    double IsospinOp(std::vector<int> tz1, std::vector<int> tz2, int npar, int ipar, 
		     VectorXi Perm, int i_spin_coupling);
    double SpinOp_pair(std::vector<int> sz1, std::vector<int> sz2, int npar, int i, int j, 
		       VectorXi Perm);
    double IsospinOp_pair(std::vector<int> tz1, std::vector<int> tz2, int npar, int i, int j, 
			  VectorXi Perm, int index);

    int factorial(int npar);
    int sign(double x);

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
    double overlap(BasisState &state1, BasisState &state2);
    double energy(BasisState &state1, BasisState &state2);
};
#endif 
