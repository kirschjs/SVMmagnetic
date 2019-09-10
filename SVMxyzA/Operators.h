#pragma once
#ifndef Operators_H
#define Operators_H
#include "Input.h"
#include "BasisState.h"
#include "MatrixElement.h"
#include <Eigen> 
using namespace Eigen;

class Operators
{
private:
    int npar,nop2b;
    SpinState spin;
    IsospinState isospin;

public:
    Operators(Input &input);
    double O(int i, int j, VectorXi Perm, int kind);
    double SpinOp(   std::vector<int> sz1, std::vector<int> sz2, int i, int j, VectorXi Perm, int kind);
    double IsospinOp(std::vector<int> tz1, std::vector<int> tz2, int i, int j, VectorXi Perm, int kind);
    double Spin2bME(   VectorXi &sz1, VectorXi &sz2, int i, int j, VectorXi &Perm, int iop);
    double Isospin2bME(VectorXi &sz1, VectorXi &sz2, int i, int j, VectorXi &Perm, int iop);
    double O_0(int i1, int i2, int i3, int i4);
    double O_1(int i1, int i2, int i3, int i4);
    std::vector<ME2bST> OST_2bme(SpinIsospinState &state1, SpinIsospinState &state2, VectorXi &Perm);
    ~Operators();
};
#endif 
