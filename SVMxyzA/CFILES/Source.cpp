
#include "Rand.h"
#include "Input.h"
#include "SVM.h"
#include "MatrixElement.h"
#include "BasisState.h"

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


int main(int argc, char* argv[])
{
    clock_t begin = clock();
    string jobname;
    if (argc < 2) { jobname = "2H_EFT"; }
    else { jobname = argv[1]; } 

/* Read and print input data */
    std::size_t found;
    found=jobname.find(".inp");
    if (found != std::string::npos) {jobname = jobname.substr(0,found);}
    found=jobname.find("input/");
    if (found != std::string::npos) {jobname = jobname.substr(found+6);}
    cout << " jobname = " << jobname << "\n";
    Input input("./input/"+jobname+".inp");
    input.print();
//===========
    ifstream  srcc("./input/"+jobname+".inp");
    cout<< srcc.rdbuf();


    ifstream  src("./input/"+jobname+".inp");
    ofstream  dst("./output/"+jobname+".txt");
    dst << src.rdbuf();
    dst.close();

//===========
//===========
        //ofstream  sta("stats.txt");
       // sta.close();
//===========


/* Initialize the random numbr generator */
    Rand rand(input.irand);
/* Initialize SVM  */
    printf("\n\t Initialize SVM \n");
    SVM svm(rand, input);
    BasisState NewState;
    vector<BasisState> Basis;
	
    NewState = svm.FirstNewState();
//    NewState.print();

    if (NewState.notdefined) {  
	cout << NewState.notdefined << "\n";
	cout << "finding new state with appropriate overloop failed" << endl << endl;
	return 0; }
/* start SVM iterations */
    Basis.push_back(NewState);
    svm.UpdateNorm(Basis);
    svm.UpdateHamiltonian(Basis);
    MatrixXd Norm;
    MatrixXd H;
    MatrixXd C;
    VectorXd D;
    double E;
    double EE;
    int n_accuracy=1;
    vector<double> dE;
    GeneralizedSelfAdjointEigenSolver<MatrixXd> ges;
    printf("\t Start SVM iters\n\n");
	int itr = 1;
	while (itr < input.maxbasis)
	  {
	    Norm = svm.NormMatrix(Basis);
	    H    = svm.HamiltonianMatrix(Basis);
//	    cout << "H\n" << H    << "\n";
//	    cout << "N\n" << Norm << "\n";
	    ges.compute(H, Norm);
	    C = ges.eigenvectors();
	    D = ges.eigenvalues();
	    E = D.minCoeff();
            if (itr == 1)  EE = E + abs(E / 2);
            dE.push_back(abs((EE - E) / E));
	    printf("\t iter = %4d     E = %14.8f    dE = %14.8f \n",itr,E,dE[itr-1]);
//==================
            dst.open("./output/"+jobname+".txt", ios::app);
            dst<<"      itr= "<<itr<<"        E= "<<fixed<<E<<endl;
            dst.close();
//==================



            if(dE[itr-1] > pow(10, -5)) n_accuracy=1;
            if(dE[itr-1] < pow(10, -5)) n_accuracy++;
            //if(n_accuracy==10) break;
            NewState = svm.NewState(Basis, C, D, E, EE);
	    if (NewState.notdefined){
	      cout << "Failed to find a new state with lower energy" << endl << endl;
	      break;}	   
	    Basis[itr] = NewState;
	    svm.UpdateNorm(Basis);
	    svm.UpdateHamiltonian(Basis);
            EE = E;  
            if(itr%5==0)
               {
                  dst.open("./output/"+jobname+".txt", ios::app);
                  dst<<"  more eigenvalues=  ";
                  cout<<"   more eigenvalues=  ";
                 // for(int ii=1; ii<itr-1; ii++)
                  int nbrevprint = itr;
                  if(itr > 30) nbrevprint = 30;
                   for(int ii=1; ii<nbrevprint; ii++)
                  {
                       dst<<D(ii)<<"  ";
                       cout<<D(ii)<<"  ";
                  }
                  dst<<endl;
                  dst.close();
                  cout<<endl;
               }
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




