
#include "Rand.h"
#include "Input.h"
#include "SVM.h"
#include "MatrixElement.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

#include <Eigen> 

#include <iomanip>
#include <string>
#include <stdlib.h> 
#include <ctime> 
#include <cstdio>
#include <cstdlib>
#include <stdio.h>

using namespace Eigen;
using namespace std;

void Erorr(int i);
int main(int argc, char* argv[])
{
clock_t begin = clock();
	string jobname;
	if (argc < 2) { jobname = "2H_EFT"; }
	else { jobname = argv[1]; } 

/* Read and print input data */
	Input input("./input/"+jobname+".inp");
	input.print();
//============================
        ifstream  srcc("./input/"+jobname+".inp");
        cout<< srcc.rdbuf();


        ifstream  src("./input/"+jobname+".inp");
        ofstream  dst("./output/"+jobname+".txt");
        dst << src.rdbuf();
        dst.close();
//============================
	Rand rand(input.irand);
	printf("\n\t Initialize SVM \n");
	SVM svm(rand, input);
	vector<MatrixXd> NewState;
	vector<vector<MatrixXd>> Basis;

	MatrixXd Norm;
	MatrixXd H;
	MatrixXd C;
	VectorXd D;
	double E;
	double EE;
        int n_accuracy=1;
        double dE;
	GeneralizedSelfAdjointEigenSolver<MatrixXd> ges;

//============================================================
int FileKey=input.file_key;
int itr;
fstream  basis_list;
int NP=input.npar;

if(FileKey==1)
{
    itr=1;
    NewState = svm.FirstNewState();
    if (NewState[0](0, 0) == 2000) {Erorr(1); return 0;}
    Basis.push_back(NewState);
    svm.UpdateNorm(Basis);
    svm.UpdateHamiltonian(Basis);
}

if(FileKey==2)
{
    string line;
    double x;
    MatrixXd Mx = MatrixXd::Zero(NP, NP);
    vector<MatrixXd> State(3);
    int ib=0, jb=0, mb=0, nb=0;
    basis_list.open("./output/Basis_"+jobname+".txt", ios::in);
       	while (getline(basis_list, line))
	{
                  x=stod(line);
                  Mx(mb,nb)=x;
                  ib++;
                  mb++; 
                  if(mb==NP){mb=0; nb++;}
                  if(nb==NP){nb=0; State[jb]=Mx; jb++;}
                  if(jb==3) jb=0;
                  if(ib==3*NP*NP) {ib=0; Basis.push_back(State); svm.UpdateNorm(Basis);  svm.UpdateHamiltonian(Basis); }              
	}
       basis_list.close();
       itr= Basis.size()+1;

       vector<vector<MatrixXd>> r_Basis;
       r_Basis=Basis;
       r_Basis.pop_back();   
       	    Norm = svm.NormMatrix(r_Basis);
	    H    = svm.HamiltonianMatrix(r_Basis);
	    ges.compute(H, Norm);
	    C = ges.eigenvectors();
	    D = ges.eigenvalues();
	    EE = D.minCoeff(); 
}

//=========================================================

	
	printf("\t Start SVM iters\n\n");
	while (itr < input.maxbasis)
	  {
	    Norm = svm.NormMatrix(Basis);
	    H    = svm.HamiltonianMatrix(Basis);
	    ges.compute(H, Norm);
	    C = ges.eigenvectors();
	    D = ges.eigenvalues();
	    E = D.minCoeff();   
        if (itr == 1)  EE = E + abs(E / 2);
        dE=abs((EE - E) / E);

	    printf("\t iter = %4d     E = %14.8f    dE = %14.8f \n",itr,E,dE);
//=======================================================================================
            dst.open("./output/"+jobname+".txt", ios::app);
            dst<<"      itr= "<<itr<<"        E= "<<fixed<<E<<"  dE= "<<fixed<<dE<<endl;
            dst.close();
//=======================================================================================
//=======================================================================================
    basis_list.open("./output/Basis_"+jobname+".txt",  ios::out | ios::trunc);
    int BS=Basis.size();
for(int ibasis=0; ibasis<BS; ibasis++)
{
    for(int i=0; i<3; i++)
       {
       for(int j=0; j<NP; j++)
          {
          for(int k=0; k<NP; k++)
             {
                basis_list.precision(17);
                basis_list<<fixed<<Basis[ibasis][i](j,k)<<endl;
             }
          }
       }
}
basis_list.close();      
//========================================================================================
            if(dE > pow(10, -5)) n_accuracy=1;
            if(dE < pow(10, -5)) n_accuracy++;
           // if(n_accuracy==10) break;
            NewState = svm.NewState(Basis, C, D, E, EE);
            if (NewState[0](0, 0) == 2000) {Erorr(1); return 0;}
            Basis.push_back(NewState);
            svm.UpdateNorm(Basis);
	    svm.UpdateHamiltonian(Basis);
	  
            if(itr%10==0)
               {
                  dst.open("./output/"+jobname+".txt", ios::app);
                  dst<<"  more eigenvalues=  ";
                  cout<<"   more eigenvalues=  ";
                  for(int ii=1; ii<itr-1; ii++)
                  //for(int ii=1; ii<4; ii++)
                  {
                       dst<<D(ii)<<"  ";
                       cout<<D(ii)<<"  ";
                  }
                  dst<<endl;
                  dst.close();
                  cout<<endl;
               }
            EE = E; 
            itr = itr + 1; 

	  }

clock_t end = clock();
double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
            dst.open("./output/"+jobname+".txt", ios::app);
            dst<< "time=  " << elapsed_secs << endl;
            dst.close();
cout << "time=  " << elapsed_secs << endl;

	return 0;
}

void Erorr(int i)
{
    	if (i==1) cout << "finding new state with appropriate overloop failed" << endl << endl;
}



