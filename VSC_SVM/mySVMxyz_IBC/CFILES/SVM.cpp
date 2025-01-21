#include "SVM.h"
#include "Input.h"
#include "MatrixElement.h"
#include "Rand.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>  // OpenMP for parallelization
#include <Eigen> 


using namespace Eigen;

SVM::SVM(Rand &r, Input &input) :rr(r), me(input)
{
	N = input.npar;                      // number of particles
	bmin = input.rndmin;				 // width-smapling interval bounds
	bmax = input.rndmax;
	mm0 = input.mm0;                     //
	kk0 = input.kk0;                     //
	int ndb = input.maxbasis;            // variable denoted "mnb" in input file 
    iBoxInf=input.keycontinue;
	Hmatrix = MatrixXd::Zero(ndb,ndb);
	Nmatrix = MatrixXd::Zero(ndb,ndb);
}

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

// Eq.(15)
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

double SVM::NewEnergy(vector<vector<MatrixXd>> Basis, MatrixXd C, VectorXd D, double E, double EE)
{
	int itr = Basis.size() - 1;
	VectorXd c = VectorXd::Zero(itr + 1);
	VectorXd Overlap = VectorXd::Zero(itr);
	VectorXd q = VectorXd::Zero(itr);
	double aa = 0;
	double NN = 0;

	/* Gram-Schmidt orthogonalization adapted from Eq.(12)
	   sum of overlaps over the N-1 Eigenvectors times the Eigenvector
	   ECCE: the |Psi> in (12) are eigenvectors which we need to expand in our basis -> C(.,.) */
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

	// normelaize the i+1 vector
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

	// create the matrix element of the i+1 vector with all previous orthogonal Eigenvectors.
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

	// solve the secular equation and thus obtain the new eigenvalue
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

// instantiate a symmetric NxN matrix with elements randomly chosen 
// from the interval [bmin,bmax]
// d_ij controls the distance between particle i and j 
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

// take the D matrix and transform it into an A matrix, i.e., 
// the corresponding quadratic form between single-particle coordinates
// sum_{i<j}^N (x_i-x_j)^2/(2dij^2)-sum_i epsilon x_i^2 = -1/2 sum a_ij x_i x_j
MatrixXd SVM::A(MatrixXd d)
{
	MatrixXd A = MatrixXd::Zero(N, N);

	for (int i = 0; i < N; i++){
	  
	  for (int j = i; j < N; j++){
	    
		if (i == j){
	      for (int k = 0; k < N; k++){
			if (i != k) A(i, j) = A(i, j) + 2 * pow(d(i, k), -2);
	      }
		  if(N==2) A(i, j) = A(i, j) + 4.83*rr.doub() * pow(bmin + (bmax - bmin)*rr.doub(), -1);
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

vector<MatrixXd> SVM::FirstNewState()
{
	vector<MatrixXd> NewState;
	
	NewState.push_back(A(Dmatrix()));  //x
	NewState.push_back(A(Dmatrix()));  //y
	NewState.push_back(A(Dmatrix()));  //z
	
	// instantiate state with non-zero norm
	double e_overlap=0;
	int count=0;
	while (count<10)
	{
	   	count++;
	   	NewState[0] = A(Dmatrix());
	   	NewState[1] = A(Dmatrix());
	   	NewState[2] = A(Dmatrix());
	   	e_overlap = me.overlap(NewState, NewState);
	   	if(e_overlap>1e-8) break;
	}
	// return if too many trails were needed
	if(count >= 8)
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
	    count = 0;
	    while (count < 1000)
		{
	  		count++;
	  		NewState[0] = A(Dmatrix());
	        NewState[1] = A(Dmatrix());
	        NewState[2] = A(Dmatrix());
	  		e_overlap = me.overlap(NewState, NewState);
			// abort if norm too small
	  		if (e_overlap < 1e-8) continue;
	  		NewE = me.energy(NewState, NewState) / e_overlap;
			// replace reference state
	  		if (NewE < MinE)
	    	{
	      		MinE = NewE;
	      		State = NewState;
	    	}
		}
		return State;
	}
}

// algorithm to add a basis state (see https://inspirehep.net/literature/398252)
vector<MatrixXd> SVM::NewState(vector<vector<MatrixXd>> Basis, MatrixXd C, VectorXd D, double E, double EE)
{

    vector<MatrixXd> NewState;

	// obtain the candidate basis state
	MatrixXd dx = Dmatrix(); //x axis
	NewState.push_back(A(dx));

	MatrixXd dy = Dmatrix();  //y axis
	NewState.push_back(A(dy));

	MatrixXd dz = Dmatrix();  //z axis
	NewState.push_back(A(dz));

	// add candidate to basis
	// ECCE: check if the update does compute only matrix elements involving the new state
	Basis.push_back(NewState);
	UpdateNorm(Basis);
	UpdateHamiltonian(Basis);
        
    int Bsize=Basis.size()-1;

	int size = NewState.size();
	vector<MatrixXd> State(size); 

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

	while(count4<=mm0)
	{
        ix=0; jx=1; kkx=0;
        iy=0; jy=1; kky=0;
        iz=0; jz=1; kkz=0;
        count1=0; count2=0; count3=0; 
        minE=E; NewE=E;       

		while (count1 < 3*mm0*kk0*N*(N - 1)/2)
		{     
			if (CheckOverlap(Basis) == 1)
			{
				// obtain the ground state from the secular equation (15)
				NewE = NewEnergy(Basis, C, D, E, EE);
				// does the additional state lower the energy?
				if (NewE < minE) 
				{
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
             
            if(count1%kk0==0)  
            { 
               dx=mindx;
               dy=mindy;
               dz=mindz;
            }
			count2++;
			if (count2 > 200) 
			{
                State[0]=NewState[0];
				State[0](0, 0) = 2000;
				break;
			}

            count3++;
			//================= X ====================
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
            //================end x===================
            //================= Y ====================
            else if(count3<=(2*kk0*N*(N - 1)/2))
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
						if (iy == N - 1) iy = 0;
						jy = iy + 1;
			  		}
		      	}
		      	NewState[1]= A(dy);
            }
            //================end y===================
            //================= Z ====================
            else
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
		      	NewState[2]= A(dz);
            }
            //================end z===================
            if(count3==3*kk0*N*(N - 1)/2) count3=0;
          
			Basis[Bsize] = NewState;
			UpdateNorm(Basis);
			UpdateHamiltonian(Basis);
		}
        if(xx==0)
		{
        	count4++;
        	dx = Dmatrix();
        	dy = Dmatrix();
        	dz = Dmatrix();
        	NewState[0]=A(dx);
        	NewState[1]=A(dy);
        	NewState[2]=A(dz);
        	Basis[Bsize] = NewState;
	 		UpdateNorm(Basis);
	 		UpdateHamiltonian(Basis);
         	State[0]=NewState[0]; State[0](0, 0) = 2000;
        }
        if(xx==1) break;
	} //end while(count4<5)
	//cout<<"count4=  "<<count4<<endl;
	return State;
}

//=============================non linear optimization=======================

// vector<MatrixXd> SVM::NewState(vector<vector<MatrixXd>> Basis, MatrixXd C, VectorXd D, double E, double EE) {
//     vector<MatrixXd> NewState;

//     // Initialize D-matrices
//     MatrixXd dx = Dmatrix(); // x-axis
//     MatrixXd dy = Dmatrix(); // y-axis
//     MatrixXd dz = Dmatrix(); // z-axis

//     // Add initial state to Basis
//     NewState.push_back(A(dx));
//     NewState.push_back(A(dy));
//     NewState.push_back(A(dz));

//     Basis.push_back(NewState);
//     UpdateNorm(Basis);
//     UpdateHamiltonian(Basis);

//     int Bsize = Basis.size() - 1;

//     // Parameters for optimization
//     double learning_rate = 0.05; // Initial step size
//     double decay_rate = 0.99;    // Step size decay
//     int max_iterations = 1000;  // Maximum number of iterations
//     double tolerance = 1e-8;    // Convergence threshold

//     double minEnergy = E;
//     double currentEnergy = E;

//     MatrixXd best_dx = dx;
//     MatrixXd best_dy = dy;
//     MatrixXd best_dz = dz;

//     // Gradient matrices
//     MatrixXd grad_dx(dx.rows(), dx.cols());
//     MatrixXd grad_dy(dy.rows(), dy.cols());
//     MatrixXd grad_dz(dz.rows(), dz.cols());

//     for (int iter = 0; iter < max_iterations; ++iter) {
//         // Update NewState with current dx, dy, dz
//         NewState[0] = A(dx);
//         NewState[1] = A(dy);
//         NewState[2] = A(dz);

//         Basis[Bsize] = NewState;
//         UpdateNorm(Basis);
//         UpdateHamiltonian(Basis);

//         currentEnergy = NewEnergy(Basis, C, D, E, EE);

//         // Update minimum energy and save best states
//         if (currentEnergy < minEnergy) {
//             minEnergy = currentEnergy;
//             best_dx = dx;
//             best_dy = dy;
//             best_dz = dz;

//             if (std::abs(currentEnergy - minEnergy) < tolerance) {
//                 break;
//             }
//         }

//         // Gradient approximation using finite difference
//         double delta = 1e-5;

//         for (int i = 0; i < dx.rows(); ++i) {
//             for (int j = 0; j < dx.cols(); ++j) {
//                 dx(i, j) += delta;
//                 Basis[Bsize][0] = A(dx);
//                 UpdateNorm(Basis);
//                 UpdateHamiltonian(Basis);
//                 double energy_dx = NewEnergy(Basis, C, D, E, EE);
//                 dx(i, j) -= delta;

//                 grad_dx(i, j) = (energy_dx - currentEnergy) / delta;
//             }
//         }

//         for (int i = 0; i < dy.rows(); ++i) {
//             for (int j = 0; j < dy.cols(); ++j) {
//                 dy(i, j) += delta;
//                 Basis[Bsize][1] = A(dy);
//                 UpdateNorm(Basis);
//                 UpdateHamiltonian(Basis);
//                 double energy_dy = NewEnergy(Basis, C, D, E, EE);
//                 dy(i, j) -= delta;

//                 grad_dy(i, j) = (energy_dy - currentEnergy) / delta;
//             }
//         }

//         for (int i = 0; i < dz.rows(); ++i) {
//             for (int j = 0; j < dz.cols(); ++j) {
//                 dz(i, j) += delta;
//                 Basis[Bsize][2] = A(dz);
//                 UpdateNorm(Basis);
//                 UpdateHamiltonian(Basis);
//                 double energy_dz = NewEnergy(Basis, C, D, E, EE);
//                 dz(i, j) -= delta;

//                 grad_dz(i, j) = (energy_dz - currentEnergy) / delta;
//             }
//         }

//         // Nonlinear gradient transform.
//         grad_dx = grad_dx.array().tanh();
//         grad_dy = grad_dy.array().tanh();
//         grad_dz = grad_dz.array().tanh();

        
//         dx -= learning_rate * grad_dx;
//         dy -= learning_rate * grad_dy;
//         dz -= learning_rate * grad_dz;

      
//         learning_rate *= decay_rate;

//         if (iter % 100 == 0) {
//             std::cout << "Iteration: " << iter << " | Energy: " << currentEnergy << std::endl;
//         }
//     }

    
//     vector<MatrixXd> State;
//     State.push_back(A(best_dx));
//     State.push_back(A(best_dy));
//     State.push_back(A(best_dz));

//     return State;
// }




//===========================end here =======================================

// calculate and return the norm matrix of a basis
MatrixXd SVM::NormMatrix(vector<vector<MatrixXd>> Basis)
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

// calculate and return the Hamilton matrix of a basis
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

// calculate the norm matrix of a basis and equate it to the
// SVM member variable Nmatrix
void SVM::UpdateNorm(vector<vector<MatrixXd>> Basis)
{
	int ncur = Basis.size()-1;
  	for (int n1 = 0; n1 <= ncur; n1++)
    {
    	Nmatrix(n1, ncur) = me.overlap(Basis[n1], Basis[ncur]);
		Nmatrix(ncur, n1) = Nmatrix(n1, ncur);
    }   
}

// calc. the hamilton matrix and set it equal to the
// designated SVM member variable
void SVM::UpdateHamiltonian(vector<vector<MatrixXd>> Basis)
{
	int ncur = Basis.size()-1;
  	for (int n1 = 0; n1 <= ncur; n1++)
    {
    	Hmatrix(n1, ncur) = me.energy(Basis[n1], Basis[ncur]);
    	Hmatrix(ncur, n1) = Hmatrix(n1, ncur);
    }    
}

SVM::~SVM()
{
}