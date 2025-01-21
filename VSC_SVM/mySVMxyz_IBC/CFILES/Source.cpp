#include "/home/sourav/svm/SVMmagnetic/VSC_SVM/mySVMxyz_IBC/CFILES/Input.h"
#include "/home/sourav/svm/SVMmagnetic/VSC_SVM/mySVMxyz_IBC/CFILES/Rand.h"
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


int main(int argc, char* argv[])
{
  	clock_t begin = clock();
	
	string jobname;
	
	// in debugging mode, load the default input
	if (argc < 2) { jobname = "default_2body"; }
	// run with specified input
	else { jobname = argv[1]; } 
	Input input("/home/sourav/svm/SVMmagnetic/VSC_SVM/mySVMxyz_IBC/input/"+jobname+".inp");
	input.print();

    ifstream  srcc("./input/"+jobname+".inp");
    cout<< srcc.rdbuf();

    ifstream  src("./input/"+jobname+".inp");
    ofstream  dst("./output/"+jobname+".txt");
    dst << src.rdbuf();
    dst.close();


	/* Initialize the random numbr generator */
	Rand rand(input.irand);
	/* Initialize SVM  */
	printf("\n\t Initialize SVM \n");
	SVM svm(rand, input);
	
	vector<MatrixXd> NewState;
	vector<vector<MatrixXd>> Basis;
	NewState = svm.FirstNewState();

	if (NewState[0](0, 0) == 2000)
	{
		cout << "finding new state with appropriate overloop failed" << endl << endl;
		return 0;
	}

	/* start SVM iterations */
	
	// after this step, the basis has dimension 1
	Basis.push_back(NewState);
	svm.UpdateNorm(Basis);
	svm.UpdateHamiltonian(Basis);

	MatrixXd Norm;
	MatrixXd H;
	MatrixXd C;
	VectorXd D;
	double E;
	double EE;
	double bmin,bmax ,b;
    int n_accuracy=1;
    vector<double> dE;
	GeneralizedSelfAdjointEigenSolver<MatrixXd> ges;
	
	printf("\t Start SVM iters\n\n");
	int itr = 1;
	while (itr < input.maxbasis)
	{
    	Norm = svm.NormMatrix(Basis);
	    H    = svm.HamiltonianMatrix(Basis);
	    ges.compute(H, Norm);
	    C = ges.eigenvectors();
	    D = ges.eigenvalues();
	    E = D.minCoeff();

        if (itr == 1)  EE = E + abs(E / 2);
        dE.push_back(abs((EE - E) / E));
	    printf("\t iter = %4d     E = %14.8f    dE = %14.8f \n",itr,E,dE[itr-1]);

        dst.open("./output/"+jobname+".txt", ios::app);
        dst<<"      itr= "<<itr<<"        E= "<<fixed<<E<<endl;
        dst.close();

        if(dE[itr-1] > pow(10, -5)) n_accuracy=1;
        if(dE[itr-1] < pow(10, -5)) n_accuracy++;
        //if(n_accuracy==10) break;

		// algorithm as in https://inspirehep.net/literature/398252
        NewState = svm.NewState(Basis, C, D, E, EE  );
	    if (NewState[0](0, 0) == 2000)
	    {
			cout << "finding new state with lower energy failed" << endl << endl;
			break;
	    }	   
	    Basis.push_back(NewState);
	    svm.UpdateNorm(Basis);
	    svm.UpdateHamiltonian(Basis);
        EE = E;  
        if(itr%5==0)
        {
			int outtmp = 5;
			if(itr > 35) outtmp = 40;
            dst.open("./output/"+jobname+".txt", ios::app);
            dst<<"  more eigenvalues=  ";
            cout<<"   more eigenvalues=  ";
            // for(int ii=1; ii<itr-1; ii++)
            for(int ii=0; ii<outtmp; ii++)
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




