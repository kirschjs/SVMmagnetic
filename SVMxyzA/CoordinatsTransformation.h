#ifndef CoordinatsTransformation_H
#define CoordinatsTransformation_H
#pragma once
#include <Eigen>
using namespace Eigen;

class CoordinatsTransformation
{
private:
        int N;
public:
	CoordinatsTransformation(int npar);
	VectorXd SingleParticle(int i, int j);
	~CoordinatsTransformation();
};

#endif
