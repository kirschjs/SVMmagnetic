#include "SVM.h"
#include "Input.h"
#include "MatrixElement.h"
#include "Rand.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen> 
using namespace Eigen;

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
int SVM::CheckOverlap(vector<vector<MatrixXd>> Basis)
{
  int itr = Basis.size() - 1;
  float vnorm,vdotv;
  vnorm = Nmatrix(itr,itr);
//  printf("\t itr = %4d     vnrom = %12.6f \n",itr,vnorm);
  if (vnorm < 1e-8) return 0;
  if (itr > 0){
    for (int i = 0; i < itr; i++){
      vdotv = Nmatrix(itr,i)/sqrt(Nmatrix(i,i)*Nmatrix(itr,itr)) ;
//      printf("\t\t i = %4d     vdotv = %12.6f \n",i,vdotv);
      if (vdotv > 0.99){
	return 0;
      }
    }
  }
  return 1;
}
//=============================================================================
double EigenValusEquation(int itr, VectorXd D, VectorXd q, double aa, double xx)
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

double SVM::NewEnergy(vector<vector<MatrixXd>> Basis, MatrixXd C, VectorXd D, double E, double EE)
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
	double Ee1 = EigenValusEquation(itr, D, q, aa, e1);
	//std::cout << "intial e1= " << e1 << ".   intial e2= " << e2<<endl;
	while (count < 101)
	{
		
		if (Ee1*EigenValusEquation(itr, D, q, aa, e2) < 0)  break;
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
			if (EigenValusEquation(itr, D, q, aa, e3)*EigenValusEquation(itr, D, q, aa, e2) < 0)  e1 = e3;
			else e2 = e3;
			count++;
			if (count > 100)  break;
		}
		//	std::cout << std::endl;
             
	}
	return e3;
	//==============================================================================

}


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

MatrixXd SVM::diagDmatrix()
{
	MatrixXd dd = MatrixXd::Zero(N,N);
	for (int i = 0; i < N; i++)
	{
		dd(i,i) = bmin + (bmax - bmin)*rr.doub();
	}
	return dd;
}

MatrixXd SVM::A(MatrixXd d)
{
	MatrixXd A = MatrixXd::Zero(N, N);

	for (int i = 0; i < N; i++){
	  for (int j = i; j < N; j++){
	    if (i == j){
	      for (int k = 0; k < N; k++){
		if (i != k) A(i, j) = A(i, j) + 2 * pow(d(i, k), -2);
	      }
	    }
	    else{
	      A(i, j) = -2 * pow(d(i, j), -2);
	      A(j, i) = A(i, j);
	    }
	  }
	}
	return A;
}

MatrixXd SVM::B(MatrixXd dd)
{
	MatrixXd B = MatrixXd::Zero(N, N);
	for (int i = 0; i < N; i++)
	{
            B(i, i) =1.0 / (dd(i,i)*dd(i,i));	
	    //if(iBoxInf==1) B(i, i) = 1.e-6; 
	}
	return B;
}

MatrixXd SVM::s(MatrixXd dd)
{
	MatrixXd s = MatrixXd::Zero(N, N);
	for (int i = 0; i < N; i++)
	{
                s(i, i) =-0.5*bmax+dd(i,i);
	}
	return s;
}

//=============================================================================
vector<MatrixXd> SVM::FirstNewState()
{
  vector<MatrixXd> NewState;
  NewState.push_back(A(Dmatrix()));
  NewState.push_back(B(diagDmatrix()));
  NewState.push_back(s(diagDmatrix()));
  NewState.push_back(A(Dmatrix()));
  NewState.push_back(B(diagDmatrix()));
  NewState.push_back(s(diagDmatrix()));
  NewState.push_back(A(Dmatrix()));
  NewState.push_back(B(diagDmatrix()));
  NewState.push_back(s(diagDmatrix()));
  double e_overlap;
  int count=0;
  while (count<10)
    {
      count++;
      NewState[0] = A(Dmatrix());
      NewState[3] = A(Dmatrix());
      NewState[6] = A(Dmatrix());
      e_overlap = me.overlap(NewState, NewState);
      if(e_overlap>1e-8) break;
    }
  if (count >= 8)
    {
      NewState[0](0, 0) = 2000;
      return NewState;
    }
  else
    {
      double MinE = me.energy(NewState, NewState) / e_overlap;
      double NewE;
      vector<MatrixXd> State(NewState.size());
      State = NewState;
      int count = 0;
      while (count < 1000)  //1000
	{
	  count++;
      NewState[0] = A(Dmatrix());
      NewState[3] = A(Dmatrix());
      NewState[6] = A(Dmatrix());
	  e_overlap = me.overlap(NewState, NewState);
	  if (e_overlap < 1e-8) continue;
	  NewE = me.energy(NewState, NewState) / e_overlap;
	  if (NewE < MinE)
	    {
	      MinE = NewE;
	      State[0] = NewState[0];
          State[3] = NewState[3];
          State[6] = NewState[6];
	    }
	}
	
     // Hmatrix(0,0)=me.energy(State, State);
     // Nmatrix(0,0)=me.overlap(State, State);
      return State;
    }
}

//=============================================================================
vector<MatrixXd> SVM::NewState(vector<vector<MatrixXd>> Basis, MatrixXd C, VectorXd D, double E, double EE)
{
        vector<MatrixXd> NewState;

	MatrixXd dx = Dmatrix();
	NewState.push_back(A(dx));
	MatrixXd ddx = diagDmatrix();
	NewState.push_back(B(ddx));
	MatrixXd ssx = diagDmatrix();
	NewState.push_back(s(ssx));

	MatrixXd dy = Dmatrix();
	NewState.push_back(A(dy));
	MatrixXd ddy = diagDmatrix();
	NewState.push_back(B(ddy));
	MatrixXd ssy = diagDmatrix();
	NewState.push_back(s(ssy));

	MatrixXd dz = Dmatrix();
	NewState.push_back(A(dz));
	MatrixXd ddz = diagDmatrix();
	NewState.push_back(B(ddz));
	MatrixXd ssz = diagDmatrix();
	NewState.push_back(s(ssz));

	Basis.push_back(NewState);
	UpdateNorm(Basis);
	UpdateHamiltonian(Basis);
        
        int Bsize=Basis.size()-1;

	int size = NewState.size();
	vector<MatrixXd> State(size); 

        MatrixXd mindx=dx;
        MatrixXd minddx=ddx;
        MatrixXd minssx=ssx; 

        MatrixXd mindy=dy;
        MatrixXd minddy=ddy;
        MatrixXd minssy=ssy; 

        MatrixXd mindz=dz;
        MatrixXd minddz=ddz;
        MatrixXd minssz=ssz;    

        int ix, jx, iix, jjx, kkx, llx, nnx;
        int iy, jy, iiy, jjy, kky, lly, nny;
        int iz, jz, iiz, jjz, kkz, llz, nnz;
        int count0, count1, count2, count3, count4;
        int xx=0;
        int yy=0;
        count4=0; 
        count0=0;
        double minE, NewE;

/*
        minE=E; NewE=E;
	while (count0 <= mm0)
	{
		if (CheckOverlap(Basis) == 1)
		{
			NewE = NewEnergy(Basis, C, D, E, EE);
			if (NewE < minE) 
                        {
                                yy=1;
				minE = NewE;
	                        State = NewState;                                         			
		        }
                }
	        NewState[0]=A(Dmatrix());
		Basis[Bsize] = NewState;
		UpdateNorm(Basis);
		UpdateHamiltonian(Basis);
                count0++;
	}
        if(yy==1) 
        {
        Basis[Bsize] = State;
	UpdateNorm(Basis);
	UpdateHamiltonian(Basis);
        }
  */    


while(count4<=mm0)
{
        ix=0; jx=1; iix=0; jjx=0; kkx=0; llx=0; nnx=0; 
        iy=0; jy=1; iiy=0; jjy=0; kky=0; lly=0; nny=0;
        iz=0; jz=1; iiz=0; jjz=0; kkz=0; llz=0; nnz=0;
        count1=0; count2=0; count3=0; 
        minE=E; NewE=E;       

	while (count1 < 3*mm0*kk0*(N*(N - 1)/2+N+N))
	{
		if (CheckOverlap(Basis) == 1)
		{
			NewE = NewEnergy(Basis, C, D, E, EE);
			if (NewE < minE) {
				minE = NewE;
                                xx=1;
                                State=NewState;

                                mindx=dx;
                                minddx=ddx;
                                minssx=ssx;

                                mindy=dy;
                                minddy=ddy;
                                minssy=ssy;

                                mindz=dz;
                                minddz=ddz;
                                minssz=ssz;
			}
			count2 = 0;
			//std::cout << "i= " << i << "   NewE= " << NewE << "   minE= " << minE << std::endl;
		}
            count1++;
             //============================

                    if(count1%kk0==0)  
                        { 
                           dx=mindx;
                           ddx=minddx;
                           ssx=minssx;

                           dy=mindy;
                           ddy=minddy;
                           ssy=minssy;

                           dz=mindz;
                           ddz=minddz;
                           ssz=minssz;
                        }

             //============================
		count2++;
                count3++;
		if (count2 > 200) {
                        State[0]=NewState[0];
			State[0](0, 0) = 2000;
			break;
		}
		
//==============================================================	

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
                if(count3<=kk0*N*(N - 1)/2)
                {                   
		      dx(ix,jx)= bmin + (bmax - bmin)*rr.doub();
		      dx(jx, ix) = dx(ix, jx);
		      kkx++;
		      if (kkx == kk0) 
                      {
			  kkx = 0;
			  jx++;
			  if (jx == N)
                          {
				ix++;
				if (ix == N - 1) ix = 0;
				jx = ix + 1;
			  }
		      }
		      NewState[0]= A(dx);
                }
                else if(count3<=kk0*(N*(N - 1)/2+N))
                {                   
		      ddx(iix,iix)= bmin + (bmax - bmin)*rr.doub();
		      llx++;
		      if (llx == kk0) 
                      {
			  llx = 0;
	                  iix++;
                          if (iix == N) 
                          {
                             iix=0;
                          }
		      }
		      NewState[1]= B(ddx);
                }
                else if(count3<=kk0*(N*(N - 1)/2+N+N))
                {
                      ssx(jjx,jjx)= bmin + (bmax - bmin)*rr.doub();
		      nnx++;
		      if (nnx == kk0) 
                      {
			  nnx = 0;
	                  jjx++;
                          if (jjx == N) 
                          {
                             jjx=0;
                             //count3=0;
                          }
		      }
		      NewState[2] = s(ssx);
                 }
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
                else if(count3<=kk0*(2*N*(N - 1)/2+N+N))
                {                   
		      dy(iy,jy)= bmin + (bmax - bmin)*rr.doub();
		      dy(jy, iy) = dy(iy, jy);
		      kky++;
		      if (kky == kk0) 
                      {
			  kky = 0;
			  jy++;
			  if (jy == N)
                          {
				iy++;
				if (iy == N - 1) iy= 0;
				jy = iy + 1;
			  }
		      }
		      NewState[3]= A(dy);
                }
                else if(count3<=kk0*(2*N*(N - 1)/2+2*N+N))
                {                   
		      ddy(iiy,iiy)= bmin + (bmax - bmin)*rr.doub();
		      lly++;
		      if (lly == kk0) 
                      {
			  lly = 0;
	                  iiy++;
                          if (iiy == N) 
                          {
                             iiy=0;
                          }
		      }
		      NewState[4]= B(ddy);
                }
                else  if(count3<=kk0*(2*N*(N - 1)/2+2*N+2*N))
                {
                      ssy(jjy,jjy)= bmin + (bmax - bmin)*rr.doub();
		      nny++;
		      if (nny == kk0) 
                      {
			  nny = 0;
	                  jjy++;
                          if (jjy == N) 
                          {
                             jjy=0;
                             //count3=0;
                          }
		      }
		      NewState[5] = s(ssy);
                 }
//yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
//zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
                else  if(count3<=kk0*(3*N*(N - 1)/2+2*N+2*N))
                {                   
		      dz(iz,jz)= bmin + (bmax - bmin)*rr.doub();
		      dz(jz, iz) = dz(iz, jz);
		      kkz++;
		      if (kkz == kk0) 
                      {
			  kkz = 0;
			  jz++;
			  if (jz == N)
                          {
				iz++;
				if (iz == N - 1) iz = 0;
				jz = iz + 1;
			  }
		      }
		      NewState[6]= A(dz);
                }
                else if(count3<=kk0*(3*N*(N - 1)/2+3*N+2*N))
                {                   
		      ddz(iiz,iiz)= bmin + (bmax - bmin)*rr.doub();
		      llz++;
		      if (llz == kk0) 
                      {
			  llz = 0;
	                  iiz++;
                          if (iiz == N) 
                          {
                             iiz=0;
                          }
		      }
		      NewState[7]= B(ddz);
                }
                else  if(count3<=kk0*(3*N*(N - 1)/2+3*N+3*N))
                {
                      ssz(jjz,jjz)= bmin + (bmax - bmin)*rr.doub();
		      nnz++;
		      if (nnz == kk0) 
                      {
			  nnz = 0;
	                  jjz++;
                          if (jjz == N) 
                          {
                             jjz=0;
                             //count3=0;
                          }
		      }
		      NewState[8] = s(ssz);
                 }
//zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
                 if(count3==3*kk0*(N*(N - 1)/2+N+N)) count3=0;


		Basis[Bsize] = NewState;
		UpdateNorm(Basis);
		UpdateHamiltonian(Basis);
	}
        if(xx==0){
         count4++;
         dx = Dmatrix();
         dy = Dmatrix();
         dz = Dmatrix();

         NewState[0]=A(dx);
         NewState[3]=A(dy);
         NewState[6]=A(dz);
         //dd = diagDmatrix();
	 //NewState[1]=B(dd);
	 //ss = diagDmatrix();
	 //NewState[2]=s(ss);
        Basis[Bsize] = NewState;
	 UpdateNorm(Basis);
	 UpdateHamiltonian(Basis);
         State[0]=NewState[0]; State[0](0, 0) = 2000;
         }
         if(xx==1) break;
} //end while(count4<mm0)
//cout<<"count4=  "<<count4<<endl;
	return State;
}

//=============================================================================
MatrixXd SVM::NormMatrix(vector<vector<MatrixXd>> Basis)
{
	int itr = Basis.size();
	MatrixXd Norm = MatrixXd::Zero(itr, itr);
	for (int n1 = 0; n1 < itr; n1++)
	{
		for (int n2 = n1; n2 < itr; n2++)
		{
			Norm(n1, n2) = me.overlap(Basis[n1], Basis[n2]);
			Norm(n2, n1) = Norm(n1, n2);
		}    
	}
	return Norm;
}

//=============================================================================
MatrixXd SVM::HamiltonianMatrix(vector<vector<MatrixXd>> Basis)
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
void SVM::UpdateNorm(vector<vector<MatrixXd>> Basis)
{
  int ncur = Basis.size()-1;
  for (int n1 = 0; n1 <= ncur; n1++)
    {
      Nmatrix(n1, ncur) = me.overlap(Basis[n1], Basis[ncur]);
      Nmatrix(ncur, n1) = Nmatrix(n1, ncur);
    }   
}
//=============================================================================
void SVM::UpdateHamiltonian(vector<vector<MatrixXd>> Basis)
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

