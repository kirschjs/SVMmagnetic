#pragma once
#include <Eigen> 
using namespace Eigen;

class BasisState
{
public:
    BasisState();
    BasisState(MatrixXd Ax, MatrixXd Ay, MatrixXd Az, int ts);
    void set(  MatrixXd Ax, MatrixXd Ay, MatrixXd Az, int ts);
    void print();
    MatrixXd Ax, Ay, Az; 
    int ts;
    bool notdefined = false;
    ~BasisState();
};

class SpinIsospinState
{
public:
    SpinIsospinState();
    SpinIsospinState(int ncmp, VectorXd coef, MatrixXi sz, MatrixXi tz);
    void set(        int ncmp, VectorXd coef, MatrixXi sz, MatrixXi tz);
    void print();
    int ncmp;
    VectorXd coef;
    MatrixXi sz, tz;
    ~SpinIsospinState();
};

