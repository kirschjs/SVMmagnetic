#include "MatrixElement.h"
#include "Input.h"
#include "Operators.h"
#include "CoordinatsTransformation.h"
#include "Permutation.h"
#include <cmath>
#include <Eigen> 
#include <vector>
using namespace Eigen;
using namespace std;

vector<VectorXd> SumOverCells(int dmax, int npar);
//vector<VectorXd> TwoBodyPairs(int npar);
vector<VectorXd>  TwoBodySingles (int npar);

int factorial(int npar);

MatrixElement::MatrixElement(Input &input)
{
        cout << "\t Initialize MatrixElement\n";

	npar = input.npar;
	h2m = input.h2m;

        eB=input.eB;
        //momega=input.momega;
       // momega=0.0000001;
	   momega = 0.5 * eB;
	dmax = input.dmax;
	nop = input.nop;
	npt = input.npt;
	ibf = input.bosefermi;  // ibf=2 for bosons.  ibf=1 for fermions
        apot3b = input.apot3b;
        vpot3b = input.vpot3b;

	nPairs = npar*(npar - 1) / 2;
	nPerm = factorial(npar);
	nConf = pow(2 * dmax + 1, npar); //number of cells configurations
	cout << "\t npar = " << npar << endl;
	Permutation P(npar); 
	PM = P.perm_matrix;
	PV = P.permutation;
	if(ibf==1) parity = P.parity; //fermions
	if (ibf == 2) for (int iperm = 0; iperm < nPerm; iperm++) parity.push_back(1); //bosons

	mass=MatrixXd::Zero(npar,npar);
        for (int imas = 0; imas < npar; imas++) mass(imas,imas)=input.mass[imas];

      //  TBP = TwoBodyPairs(npar);
          TBP = TwoBodySingles(npar);
	PrepareSpinIsospinME(input);
	PreparePotential(input);
                        
}

double MatrixElement::overlap(std::vector<MatrixXd> state1, std::vector<MatrixXd> state2)
{
	MatrixXd A1x = state1[0];
	MatrixXd A1y = state1[1];
	MatrixXd A1z = state1[2];


	MatrixXd A2_0x = state2[0];
	MatrixXd A2_0y = state2[1];
	MatrixXd A2_0z = state2[2];

	MatrixXd A2x;
	MatrixXd A2y;
	MatrixXd A2z;

double detx, dety, detz;
long double x;

        double overlap=0;

	for (int iperm = 0; iperm < nPerm; iperm++)   //sum over permutation
	{
	  A2x = PM[iperm].transpose() * A2_0x * PM[iperm];

	  A2y = PM[iperm].transpose() * A2_0y * PM[iperm];

	  A2z = PM[iperm].transpose() * A2_0z * PM[iperm];


          detx=(A1x+A2x).determinant();
          dety=(A1y+A2y).determinant();
          detz=(A1z+A2z).determinant();
      
          x=1/sqrt(detx*dety*detz);//*detz);
	  overlap = overlap + parity[iperm]* stme[iperm*nPairs*nop] *x;
	}
	return overlap;

}


//==============================================================================================
double MatrixElement::energy(std::vector<MatrixXd> state1, std::vector<MatrixXd> state2)
{

	MatrixXd A1x = state1[0];
	MatrixXd A1y = state1[1];
	MatrixXd A1z = state1[2];


	MatrixXd A2_0x = state2[0];
	MatrixXd A2_0y = state2[1];
	MatrixXd A2_0z = state2[2];


	MatrixXd A2x;
	MatrixXd A2y;
	MatrixXd A2z;
  
        MatrixXd InvAAx,sx;
        MatrixXd InvAAy,sy;
        MatrixXd InvAAz,sz;

        MatrixXd TTx;
        MatrixXd TTy;
        MatrixXd TTz;

        double detx, dety, detz;
        long double x, y, z;
  
  double PotEnergy = 0, KinEnergy = 0, MagneticEnergy=0, harmonic=0, PotEnergy3B = 0, PotEnergy3BP = 0;

                //===========3-body===================

                MatrixXd Bx(2,2);
                MatrixXd By(2,2);
                MatrixXd Bz(2,2);
                MatrixXd I = MatrixXd::Identity(2,2);
                CoordinatsTransformation v(npar);
                VectorXd Cik(npar);
                VectorXd Cjk(npar);
                int i1, j1, k1;
                //====================================

  for (int iperm = 0; iperm < nPerm; iperm++)   //sum over permutation
    {
      A2x = PM[iperm].transpose() * A2_0x * PM[iperm];
      A2y = PM[iperm].transpose() * A2_0y * PM[iperm];
      A2z = PM[iperm].transpose() * A2_0z * PM[iperm];

      InvAAx = (A1x + A2x).inverse();
      InvAAy = (A1y + A2y).inverse();
      InvAAz = (A1z + A2z).inverse();

      detx=(A1x+A2x).determinant();
      dety=(A1y+A2y).determinant();
      detz=(A1z+A2z).determinant();
      
      x=1/sqrt(detx*dety*detz);//*detz);


      // Kinetic Energy 
      TTx = A1x * InvAAx * A2x;
      TTy = A1y * InvAAy * A2y;
      TTz = A1z * InvAAz * A2z;
      KinEnergy = KinEnergy + parity[iperm] * stme[iperm*nPairs*nop] *x *(TTx.trace()+TTy.trace()+TTz.trace());//+TTz.trace());
//===============================================================
//double sx,sy,sz, xme; 
//=============MagneticEnergy====================================     
// for (int ipair = 0; ipair < nPairs; ipair++)   //sum over pairs
// 	{
//     sx = TBP[ipair].transpose() * InvAAx * TBP[ipair];
//     sy = TBP[ipair].transpose() * InvAAy * TBP[ipair];
//     sz = TBP[ipair].transpose() * InvAAz * TBP[ipair];
//     MagneticEnergy = MagneticEnergy + parity[iperm] * stme[iperm*nPairs*nop] *x* (sx+sy+sz);
//     harmonic =        harmonic +      parity[iperm] * stme[iperm*nPairs*nop] *x*  sz;
//     }
//================================================================
 
// for (int ipart = 0; ipart < npar; ipart++)   // sum over particles
//         {
            
//             VectorXd singleParticleX = v.SingleParticle(ipart);
//             VectorXd singleParticleY = v.SingleParticle(ipart);
//             VectorXd singleParticleZ = v.SingleParticle(ipart);

//             // Magnetic Energy calculation for individual particles
//             double sx = singleParticleX.transpose() * InvAAx * singleParticleX;
//             double sy = singleParticleY.transpose() * InvAAy * singleParticleY;
//             double sz = singleParticleZ.transpose() * InvAAz * singleParticleZ;
// 			MagneticEnergy = MagneticEnergy + parity[ipart] *x* (sx+sy+sz);
// //     harmonic =        harmonic +      parity[iperm] * stme[iperm*nPairs*nop] *x*  sz;
//      }
//=============================================================
  // for (int ipart = 0 ; ipart <npar ; ipart++)
 
   sx = InvAAx;
   sy = InvAAy;
   sz = InvAAz;
   MagneticEnergy += parity[iperm] * stme[iperm*nPairs*nop]*  (sx.trace()+sy.trace()+sz.trace())*x;//+sz.trace())  * x;
 //*(InvAAx.trace()+InvAAy.trace()+InvAAz.trace());


//======2-body===================================================
    //   for (int ipair = 0; ipair < nPairs; ipair++)   //sum over pairs
	// {
    //   sx = TBP[ipair].transpose() * InvAAx * TBP[ipair];
	//   sy = TBP[ipair].transpose() * InvAAy * TBP[ipair];
	//   sz = TBP[ipair].transpose() * InvAAz * TBP[ipair];
	//   for (int iop = 0; iop < nop; iop++)
	//     {
	//       xme = parity[iperm] * stme[(iperm*nPairs+ipair)*nop+iop];
	//       for (int ipt = 0; ipt < npt; ipt++)
	// 	{
    //               y=1/sqrt((2.0*apot(iop, ipt)*sx + 1)*(2.0*apot(iop, ipt)*sy + 1)*(2.0*apot(iop, ipt)*sz + 1));
	// 	  PotEnergy = PotEnergy + xme * vpot(iop, ipt)*x*y;
	// 	}
	//     }
	// }
//===============================================================

    //   if (npar>2){

	//   PotEnergy3BP = 0;
	//   for (int i = 0; i < npar; i++){
	//       for (int j = i + 1; j < npar; j++){
	// 	  for (int k = j + 1; k < npar; k++){
	// 	      for (int cyc = 0; cyc < 3; cyc++)
	// 		{

	// 		  if (cyc == 0) {i1 = i; j1 = j; k1 = k;}
	// 		  if (cyc == 1) {i1 = j; j1 = k; k1 = i;}
	// 		  if (cyc == 2) {i1 = k; j1 = i; k1 = j;}
	// 		  Cik = v.SingleParticle(i1, k1);
	// 		  Cjk = v.SingleParticle(j1, k1);


	// 		  Bx(0,0) = Cik.transpose()*InvAAx*Cik;
	// 		  Bx(0,1) = Cik.transpose()*InvAAx*Cjk;
	// 		  Bx(1,0) = Cjk.transpose()*InvAAx*Cik;
	// 		  Bx(1,1) = Cjk.transpose()*InvAAx*Cjk;

	// 		  By(0,0) = Cik.transpose()*InvAAy*Cik;
	// 		  By(0,1) = Cik.transpose()*InvAAy*Cjk;
	// 		  By(1,0) = Cjk.transpose()*InvAAy*Cik;
	// 		  By(1,1) = Cjk.transpose()*InvAAy*Cjk;

	// 		  Bz(0,0) = Cik.transpose()*InvAAz*Cik;
	// 		  Bz(0,1) = Cik.transpose()*InvAAz*Cjk;
	// 		  Bz(1,0) = Cjk.transpose()*InvAAz*Cik;
	// 		  Bz(1,1) = Cjk.transpose()*InvAAz*Cjk;

    //                       z=1/sqrt((2.0*apot3b*Bx + I).determinant()*(2.0*apot3b*By + I).determinant()*(2.0*apot3b*Bz + I).determinant());                          
	// 		  PotEnergy3BP = PotEnergy3BP + z*x;
			  
	// 		}
	// 	    }
	// 	}
	//     }
	//   PotEnergy3B = PotEnergy3B + parity[iperm] * stme[iperm*nPairs*nop] *PotEnergy3BP;
	//   //====================================
	// }
    }    

     // PotEnergy3B = vpot3b * PotEnergy3B;

      KinEnergy = 0.5 * h2m * KinEnergy;
     // MagneticEnergy = 0.25 * pow(npar, -1)*0.5 * h2m * eB *       eB *    MagneticEnergy; 
   //   harmonic =              pow(npar, -1)*0.5 * h2m * momega  * momega * harmonic;
      MagneticEnergy =  0.5* h2m* eB *  eB *    MagneticEnergy;
 // return PotEnergy+MagneticEnergy+harmonic+KinEnergy+PotEnergy3B;
  return MagneticEnergy+KinEnergy;
} 

//=============================================================================



//=============================================================================

// vector<VectorXd> TwoBodyPairs(int npar)
// {
// 	vector<VectorXd> Cij(npar*(npar - 1) / 2);
// 	CoordinatsTransformation v(npar);
// 	int ipair = 0;
// 	for (int i = 0; i < npar; i++)	{
// 		for (int j = i + 1; j < npar; j++){
// 			Cij[ipair] = v.SingleParticle(i, j);
// 			ipair++;
// 		}
// 	}
// 	return Cij;
// }

//=========================================
vector<VectorXd> TwoBodySingles(int npar)
{
    vector<VectorXd> Cij(npar); 
    CoordinatsTransformation v(npar);
    
    for (int i = 0; i < npar; i++) {
       
        Cij[i] = v.SingleParticle(i);
    }
    
    return Cij;
}


//=============================================================================
//=============================================================================
void MatrixElement::PreparePotential(Input &input)
{
	vpot = MatrixXd::Zero(nop, npt);
	apot = MatrixXd::Zero(nop, npt);
//	cout << "PreparePotential: nop= " << nop << "  nterms= " << npt << "\n";
	for (int ipt = 0; ipt < npt; ipt++)
	{
		for (int iop = 0; iop < nop; iop++)
		{
			vpot(iop, ipt) = input.potop[iop].vpot[ipt];
			apot(iop, ipt) = input.potop[iop].aquad[ipt];
//			printf("\t\t iop = %4d   iterm = %4d   vpot = %8.3f   "
//		    "aquad = %8.3f   \n",iop,ipt,vpot(iop,ipt),apot(iop, ipt));
		}
	}
}
//=============================================================================
void MatrixElement::PrepareSpinIsospinME(Input &input)
{
        Operators operators(input);
	stme.resize(nPairs*nPerm*nop);
        //int keypr = 1;

	int ipair = -1;
	for (int ip = 0; ip < npar; ip++){
	  for (int jp = ip + 1; jp < npar; jp++){
	    ipair++;
	//    fprintf(input.printfile,
	//	    "\t\t PrepareSpinIsospinME: ipair = %4d   ipar = %4d   jpar = %4d \n",ipair,ip,jp);
	    for (int iperm = 0; iperm < nPerm; iperm++){
	      for (int iop = 0; iop < nop; iop++){       
		stme[(iperm*nPairs+ipair)*nop+iop] = operators.O(ip, jp, PV[iperm], iop);                               
               
		/* print spin-isospin matrix element */
		//if (keypr == 1) {
		 // fprintf(input.printfile,"\t\t   iperm = %4d  ",iperm);
		 // for (int i = 0; i < npar; i++) fprintf(input.printfile,"%1d",PV[iperm][i]);
		//  fprintf(input.printfile,
		//	 "      iop = %4d  me = %9.5f \n",iop,stme[(iperm*nPairs+ipair)*nop+iop]);
		//}

	      }
	    }
	  }
	}
}
//=============================================================================
int factorial(int npar)
{
	int N;
	if (npar <= 1) return 1;
	N = npar * factorial(npar - 1);
	return N;
}



MatrixElement::~MatrixElement()
{
}

