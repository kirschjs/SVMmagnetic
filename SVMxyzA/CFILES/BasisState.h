#pragma once
#include <Eigen> 
using namespace Eigen;

class BasisState
{
public:
    BasisState();
    BasisState(MatrixXd Ax, MatrixXd Ay, MatrixXd Az, VectorXd ts_vector);
    void set(  MatrixXd Ax, MatrixXd Ay, MatrixXd Az, VectorXd ts_vector);
    void print();
    MatrixXd Ax, Ay, Az; 
    VectorXd ts_vector;
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

