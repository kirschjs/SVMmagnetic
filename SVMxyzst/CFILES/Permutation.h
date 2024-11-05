#pragma once
#ifndef Permutation_H
#define Permutation_H
#include <vector>
#include <Eigen> 
using namespace Eigen;

class Permutation
{
private:
	int npar;
public:
	Permutation(int N);
	std::vector<MatrixXd> perm_matrix;
	std::vector<VectorXi> permutation;
	std::vector<int> parity;
	~Permutation();
};
#endif 
