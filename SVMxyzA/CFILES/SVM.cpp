#include "SVM.h"
#include "Input.h"
#include "MatrixElement.h"
#include "Rand.h"
#include "BasisState.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen> 
using namespace Eigen;

//=============================================================================
SVM::SVM(Rand &r, Input &input) :rr(r), me(input)
{
	N = input.npar;
	bmin = input.rndmin;
	bmax = input.rndmax;
	mm0 = input.mm0;
	kk0 = input.kk0;
	int ndb = input.maxbasis;
        iBoxInf=input.keycontinue;
	Hmatrix = MatrixXd::Zero(ndb,ndb);
	Nmatrix = MatrixXd::Zero(ndb,ndb);
}
//=============================================================================
int SVM::CheckOverlap(vector<BasisState> &Basis)
{
  int itr = Basis.size() - 1;
  float vnorm,vdotv;
  vnorm = Nmatrix(itr,itr);
  if (vnorm < 1e-8) return 0;
  if (itr > 0){
    for (int i = 0; i < itr; i++){
      vdotv = Nmatrix(itr,i)/sqrt(Nmatrix(i,i)*Nmatrix(itr,itr)) ;
      if (vdotv > 0.99){
	return 0;
      }
    }
  }
  return 1;
}
//=============================================================================
double EigenValuesEquation(int itr, VectorXd D, VectorXd q, double aa, double xx)
{
	double vv = 1;
	double ww = 1;
	double yy = 0;
	double zz = 0;

	for (int n1 = 0; n1 < itr; n1++)
	{
		vv = vv * (D(n1) - xx);
		ww = 1;
		for (int n2 = 0; n2 < itr; n2++)
		{
			if (n2 != n1) ww = ww*(D(n2) - xx);
		}
		yy = yy + q(n1)*q(n1)*ww;

	}
	zz = (aa - xx)*vv - yy;

	return zz;
}
//=============================================================================
double SVM::NewEnergy(vector<BasisState> &Basis, MatrixXd C, VectorXd D, double E, double EE)
{
	int itr = Basis.size() - 1;
	VectorXd c = VectorXd::Zero(itr + 1);
	VectorXd Overlap = VectorXd::Zero(itr);
	VectorXd q = VectorXd::Zero(itr);
	double aa = 0;
	double NN = 0;


	//create the vector of i+1 state ortogonal to all provius ortogonal eigenvectors.

	for (int k1 = 0; k1 < itr; k1++)
	{
		for (int k2 = 0; k2 < itr; k2++)
		{
			Overlap(k1) = Overlap(k1) + C(k2, k1)* Nmatrix(itr,k2);
		}
	}

	c(itr) = 1;
	for (int k1 = 0; k1 < itr; k1++)
	{
		for (int k2 = 0; k2 < itr; k2++)
		{
			c(k1) = c(k1) - Overlap(k2)*C(k1, k2);
		}
	}

	//normelaize the i+1 vector  ==================================================

	for (int k1 = 0; k1 < itr + 1; k1++)
	{
		for (int k2 = 0; k2 < itr + 1; k2++)
		{
			NN = NN + c(k1) * c(k2) * Nmatrix(k1,k2);
		}
	}
	for (int k1 = 0; k1 < itr + 1; k1++)
	{
		c(k1) = c(k1) / sqrt(NN);
	}


	//creat the matrix element of i+1 vector with all prvius ortogonal eigenvectors.
	for (int k1 = 0; k1 < itr; k1++)
	{
	  for (int k2 = 0; k2 < itr; k2++)
	    {
		for (int k3 = 0; k3 < itr + 1; k3++)
		  {
		    q(k1) = q(k1) + Hmatrix(k2,k3)*c(k3) *C(k2, k1);
		  }
	    }
	}
	for (int k1 = 0; k1 < itr + 1; k1++)
	{
		for (int k2 = 0; k2 < itr + 1; k2++)
		{
			aa = aa + c(k1) * c(k2) * Hmatrix(k1,k2);
		}
	}

	//solving the equation for the new eigenvalue===============================================  
	int count = 0;
	double e1 = E;
	double e2 = E - abs(0.5*(E - EE));
	double e3 = E;
	double Ee1 = EigenValuesEquation(itr, D, q, aa, e1);
	//std::cout << "intial e1= " << e1 << ".   intial e2= " << e2<<endl;
	while (count < 101)
	{
		
		if (Ee1*EigenValuesEquation(itr, D, q, aa, e2) < 0)  break;
		else
		{
			e1 = e2;
			e2 = e2 - abs(0.5*(E - EE));
			count++;
		}
                
	}
	//std::cout  << "counter= " << count << std::endl;
	//if (count > 100)  std::cout << "finding root lees then the last fail " << std::endl;
        if (count <= 100)
	{
		count = 0;
		//	std::cout << "e3= ";
		while (abs((e1 - e2) / e2) > abs(1e-5*(E - EE) / EE))
		{
			e3 = (e1 + e2) / 2;
			//	std::cout << e3 << "  ";
			if (EigenValuesEquation(itr, D, q, aa, e3)*EigenValuesEquation(itr, D, q, aa, e2) < 0)  e1 = e3;
			else e2 = e3;
			count++;
			if (count > 100)  break;
		}
		//	std::cout << std::endl;
             
	}
	return e3;
	//==============================================================================

}

//=============================================================================
MatrixXd SVM::Dmatrix()
{
	MatrixXd d = MatrixXd::Zero(N, N);
	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			d(i, j) = bmin + (bmax - bmin)*rr.doub();
			d(j, i) = d(i, j);
		}
	}
	return d;
}
//=============================================================================
MatrixXd SVM::A(MatrixXd d)
{
	MatrixXd A = MatrixXd::Zero(N, N);

	for (int i = 0; i < N; i++){
	  for (int j = i; j < N; j++){
	    if (i == j){
	      for (int k = 0; k < N; k++)
                {
		if (i != k) A(i, j) = A(i, j) + 2 * pow(d(i, k), -2);
	        }
              if(N==2) A(i, j) = A(i, j) + 0.1*rr.doub() * pow(bmin + (bmax - bmin)*rr.doub(), -2);
               
	    }
	    else{
	      A(i, j) = -2 * pow(d(i, j), -2);
	      A(j, i) = A(i, j);
	    }
	  }
	}
	for (int i = 0; i < N; i++)
	{
            A(i, i) = A(i, i) + 1.e-6;	
	}
	return A;
}


//=============================================================================
BasisState SVM::FirstNewState(){
    VectorXd ts(1); ts(0)=1.;
    BasisState NewState(A(Dmatrix()),A(Dmatrix()),A(Dmatrix()),ts);

    double e_overlap=0;
    int count=0;
    while (count<10){
	count++;
	VectorXd ts(1); ts(0)=1.;
	NewState.set(A(Dmatrix()),A(Dmatrix()),A(Dmatrix()),ts);
	e_overlap = me.overlap(NewState, NewState);
	if(e_overlap>1e-8) break;
    }

    if(count >= 8){NewState.notdefined = true; return NewState;}

    double MinE = me.energy(NewState, NewState) / e_overlap;
    double NewE;
    BasisState State;
    State = NewState;
    count = 0;
    while (count < 1000){
	count++;
	VectorXd ts(1); ts(0)=1.;
	NewState.set(A(Dmatrix()),A(Dmatrix()),A(Dmatrix()),ts);
	e_overlap = me.overlap(NewState, NewState);
	if (e_overlap < 1e-8) continue;
	NewE = me.energy(NewState, NewState) / e_overlap;
	if (NewE < MinE){ MinE = NewE; State = NewState;}
    }
    return State;
}

//=============================================================================
BasisState SVM::NewState(vector<BasisState> &Basis, MatrixXd C, VectorXd D, double E, double EE)
{

        BasisState NewState;

	MatrixXd dx = Dmatrix();  //x axis
	MatrixXd dy = Dmatrix();  //y axis
	MatrixXd dz = Dmatrix();  //z axis
	VectorXd ts(1); ts(0)=1.;
	NewState.set(A(dx),A(dy),A(dz),ts);

	Basis.push_back(NewState);
	UpdateNorm(Basis);
	UpdateHamiltonian(Basis);
        
        int Bsize=Basis.size()-1;

	BasisState State; 

        MatrixXd mindx=dx;
        MatrixXd mindy=dy;
        MatrixXd mindz=dz;
 

        int ix,jx,kkx;
        int iy,jy,kky;
        int iz,jz,kkz;
        int count1,count2,count3, count4;
        int xx=0;
        count4=0; 
        double minE, NewE;

	while(count4<=mm0){
	    ix=0; jx=1; kkx=0;
	    iy=0; jy=1; kky=0;
	    iz=0; jz=1; kkz=0;
	    count1=0; count2=0; count3=0; 
	    minE=E; NewE=E;       

	    while (count1 < 3*mm0*kk0*N*(N - 1)/2){
		if (CheckOverlap(Basis) == 1){
		    NewE = NewEnergy(Basis, C, D, E, EE);
		    if (NewE < minE) {
		      minE = NewE;
		      xx=1;
		      State=NewState;

		      mindx=dx;
		      mindy=dy;
		      mindz=dz;
		    }
		    count2 = 0;
		}
                count1++;
		//============================
		if(count1%kk0==0){ 
		  dx=mindx;
		  dy=mindy;
		  dz=mindz;
		}
		//============================
		count2++;
		if (count2 > 200){
		  State.Ax=NewState.Ax;
		  State.notdefined = true;
		  break;
		}
		//============ Ax ===================
		count3++;
                if(count3<=kk0*N*(N - 1)/2){                   
		    dx(ix,jx)= bmin + (bmax - bmin)*rr.doub();
		    dx(jx, ix) = dx(ix, jx);
		    kkx++;
		    if (kkx == kk0){
		      kkx = 0;
		      jx++;
		      if (jx == N){
			ix++;
			if (ix == N - 1) ix = 0;
			jx = ix + 1;
		      }
		    }
		    NewState.Ax = A(dx);
                }
		//============ Ay ===================
                else if(count3<=(2*kk0*N*(N - 1)/2)){                   
		    dy(iy,jy)= bmin + (bmax - bmin)*rr.doub();
		    dy(jy, iy) = dy(iy, jy);
		    kky++;
		    if (kky == kk0){
			kky = 0;
			jy++;
			if (jy == N){
			  iy++;
			  if (iy == N - 1) iy = 0;
			  jy = iy + 1;
			}
		    }
		    NewState.Ay = A(dy);
                }
		//============ Az ===================
                else{                   
		  dz(iz,jz)= bmin + (bmax - bmin)*rr.doub();
		  dz(jz, iz) = dz(iz, jz);
		  kkz++;
		  if (kkz == kk0){ 
		    kkz = 0;
		    jz++;
		    if (jz == N){
		      iz++;
		      if (iz == N - 1) iz = 0;
		      jz = iz + 1;
		    }
		  }
		  NewState.Az = A(dz);
                }
                //================end z===================
                if(count3==3*kk0*N*(N - 1)/2) count3=0;
          
		Basis[Bsize] = NewState;
		UpdateNorm(Basis);
		UpdateHamiltonian(Basis);
	}
        if(xx==0){
         count4++;
         dx = Dmatrix();
         dy = Dmatrix();
         dz = Dmatrix();
	 VectorXd ts(1); ts(0)=1.;
         NewState.set(A(dx),A(dy),A(dz),ts);
	 Basis[Bsize]=NewState;
	 UpdateNorm(Basis);
	 UpdateHamiltonian(Basis);
         State.Ax=NewState.Ax; State.notdefined = true;
         }
         if(xx==1) break;
} 
	return State;
}
//=============================================================================
MatrixXd SVM::NormMatrix(vector<BasisState> &Basis)
{
	int itr = Basis.size();
	MatrixXd Norm = MatrixXd::Zero(itr, itr);
	for (int n1 = 0; n1 < itr; n1++)
	{
		for (int n2 = n1; n2 < itr; n2++)
		{
			Norm(n1, n2) = Nmatrix(n1,n2);
			Norm(n2, n1) = Norm(n1, n2);
		}    
	}
	return Norm;
}
//=============================================================================
MatrixXd SVM::HamiltonianMatrix(vector<BasisState> &Basis)
{
	int itr = Basis.size();
	MatrixXd H = MatrixXd::Zero(itr, itr);
	for (int n1 = 0; n1 < itr; n1++)
	{
		for (int n2 = n1; n2 < itr; n2++)
		{
			H(n1, n2) = Hmatrix(n1,n2);
			H(n2, n1) = H(n1,n2);
		}
	}
	return H;
}
//=============================================================================
void SVM::UpdateNorm(vector<BasisState> &Basis)
{
//  cout << Basis.size() << "\n";
  int ncur = Basis.size()-1;
  for (int n1 = 0; n1 <= ncur; n1++)
    {
//      cout << Basis.size() << " " << n1 << "\n";
      Nmatrix(n1, ncur) = me.overlap(Basis[n1], Basis[ncur]);
      Nmatrix(ncur, n1) = Nmatrix(n1, ncur);
    }   
}
//=============================================================================
void SVM::UpdateHamiltonian(vector<BasisState> &Basis)
{
  int ncur = Basis.size()-1;
  for (int n1 = 0; n1 <= ncur; n1++)
    {
      Hmatrix(n1, ncur) = me.energy(Basis[n1], Basis[ncur]);
      Hmatrix(ncur, n1) = Hmatrix(n1, ncur);
    }    
}
//=============================================================================
SVM::~SVM()
{
}
//=============================================================================

