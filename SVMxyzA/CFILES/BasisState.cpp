#include "BasisState.h"
#include <iostream>
#include <Eigen> 
using namespace Eigen;
using namespace std;

/* define class: BasisState */
//=============================================================================
BasisState::BasisState(){notdefined = false;}
//=============================================================================
BasisState::BasisState( 
		       MatrixXd Ax_, MatrixXd Ay_, MatrixXd Az_, 
		       int ts_): 
    Ax(Ax_),Ay(Ay_),Az(Az_),ts(ts_)
{notdefined = false;};
//=============================================================================
void BasisState::set(MatrixXd Ax_, MatrixXd Ay_, MatrixXd Az_, int ts_){
    Ax = Ax_ ; Ay = Ay_ ; Az = Az_;
    ts = ts_;
    notdefined = false;
};
//============================================================================= 
void BasisState::print(){ 
    IOFormat OctaveFmt(StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
    cout << "Ax=\n" << Ax.format(OctaveFmt) << "\n" ;
    cout << "Ay=\n" << Ay.format(OctaveFmt) << "\n" ;
    cout << "Az=\n" << Az.format(OctaveFmt) << "\n" ;
    cout << "isospin,spin state=" << ts << "\n";
};
//============================================================================= 
BasisState::~BasisState(){};
//=============================================================================

/* define class: SpinIsospinState */
//=============================================================================
SpinIsospinState::SpinIsospinState(){}
//=============================================================================
SpinIsospinState::SpinIsospinState(int ncmp_, VectorXd coef_, 
				   MatrixXi sz_, MatrixXi tz_):
    ncmp(ncmp_),coef(coef_),sz(sz_),tz(tz_)
{};
//=============================================================================
void SpinIsospinState::set(int ncmp_, VectorXd coef_, 
			   MatrixXi sz_, MatrixXi tz_){
    ncmp = ncmp_; 
    coef = coef_;
    sz   = sz_;
    tz   = tz_;
};
//============================================================================= 
void SpinIsospinState::print(){ 
    IOFormat OctaveFmt(StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
    cout << "SpinIsospinState Coefs=\n" << coef.format(OctaveFmt) << "\n" ;
    cout << "SpinIsospinState sz=\n" << sz.format(OctaveFmt) << "\n" ;
    cout << "SpinIsospinState tz=\n" << tz.format(OctaveFmt) << "\n" ;
};
//============================================================================= 
SpinIsospinState::~SpinIsospinState(){};
