#include "CoordinatsTransformation.h"
#include <Eigen>
using namespace Eigen;


CoordinatsTransformation::CoordinatsTransformation(int npar)
{
     N=npar;
}
VectorXd CoordinatsTransformation::SingleParticle(int i, int j)
{
	VectorXd C = VectorXd::Zero(N);

	C(i) = 1.0;
	C(j) = -1.0;

	return C;
}

CoordinatsTransformation::~CoordinatsTransformation()
{
}
