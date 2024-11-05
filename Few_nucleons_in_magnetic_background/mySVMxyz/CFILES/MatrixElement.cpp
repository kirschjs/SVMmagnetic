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
vector<VectorXd> TwoBodyPairs(int npar);

int factorial(int npar);

MatrixElement::MatrixElement(Input &input)
{
        cout << "\t Initialize MatrixElement\n";
	npar = input.npar;
	h2m = input.h2m;
        eB = input.eB;
        momega = input.momega;
	dmax = input.dmax;
	nop = input.nop;
	npt = input.npt;
	ibf = input.bosefermi;  // ibf=2 for bosons.  ibf=1 for fermions
        L = input.BoxSize;
        apot3b = input.apot3b;
        vpot3b = input.vpot3b;
	nPairs = npar*(npar - 1) / 2;
	nPerm = factorial(npar);
	nConf = pow(2 * dmax + 1, npar); //number of cells configurations

	Permutation P(npar); 
	PM = P.perm_matrix;
	PV = P.permutation;
	if(ibf==1) parity = P.parity; //fermions
	if (ibf == 2) for (int iperm = 0; iperm < nPerm; iperm++) parity.push_back(1); //bosons
	
	mass=MatrixXd::Zero(npar,npar);
        for (int imas = 0; imas < npar; imas++) mass(imas,imas)=input.mass[imas];

        TBP = TwoBodyPairs(npar);

	PrepareSpinIsospinME(input);
	PreparePotential(input);

	bbox = SumOverCells(dmax, npar);

        iBoxInf=input.keycontinue;   //iBoxInf=1 for INFINITY BC. iBoxInf=2 for BOX BC. iBoxInf is "ico" in the input file.
                        
}


//=============================================================================

double MatrixElement::overlap(std::vector<MatrixXd> state1, std::vector<MatrixXd> state2)
{
    /*
	MatrixXd A1x = state1[0];
	MatrixXd B1x = state1[1];
	VectorXd s1x = state1[2].diagonal();

	MatrixXd A2_0x = state2[0];
	MatrixXd B2_0x = state2[1];
	VectorXd s2_0x = state2[2].diagonal();

	MatrixXd A2x;
	MatrixXd B2x;
    VectorXd s2x;
    
    
    MatrixXd A1y = state1[3];
	MatrixXd B1y = state1[4];
	VectorXd s1y = state1[5].diagonal();

	MatrixXd A2_0y = state2[3];
	MatrixXd B2_0y = state2[4];
	VectorXd s2_0y = state2[5].diagonal();

	MatrixXd A2y;
	MatrixXd B2y;
    VectorXd s2y;
    
    
    MatrixXd A1z = state1[6];
	MatrixXd B1z = state1[7];
	VectorXd s1z = state1[8].diagonal();

	MatrixXd A2_0z = state2[6];
	MatrixXd B2_0z = state2[7];
	VectorXd s2_0z = state2[8].diagonal();

	MatrixXd A2z;
	MatrixXd B2z;
    VectorXd s2z;
    */
        MatrixXd A1x = state1[0];
	MatrixXd B1x = state1[1];
	VectorXd s1x = state1[2].diagonal();

	MatrixXd A2_0x = state2[0];
	MatrixXd B2_0x = state2[1];
	VectorXd s2_0x = state2[2].diagonal();

	MatrixXd A2x;
	MatrixXd B2x;
        VectorXd s2x;
    
    
        MatrixXd A1y = state1[3];
	MatrixXd B1y = state1[4];
	VectorXd s1y = state1[5].diagonal();

	MatrixXd A2_0y = state2[3];
	MatrixXd B2_0y = state2[4];
	VectorXd s2_0y = state2[5].diagonal();

	MatrixXd A2y;
	MatrixXd B2y;
        VectorXd s2y;
    
    
        MatrixXd A1z = state1[6];
	MatrixXd B1z = state1[7];
	VectorXd s1z = state1[8].diagonal();

	MatrixXd A2_0z = state2[6];
	MatrixXd B2_0z = state2[7];
	VectorXd s2_0z = state2[8].diagonal();
    
    MatrixXd A2z;
    MatrixXd B2z;
    VectorXd s2z;

    double overlap = 0;

    MatrixXd AAx;
    MatrixXd InvAAx;
     	
    MatrixXd AAy;
    MatrixXd InvAAy;
       	
    MatrixXd AAz;
    MatrixXd InvAAz;  
    
    VectorXd dx;
    VectorXd dy;
    VectorXd dz;

	long double x1x, x1y, x1z, x2x, x2y, x2z, x3x, x3y, x3z, x4x, x4y, x4z;
	long double overlapBx = 0;
        long double overlapBy = 0;
        long double overlapBz = 0;

	for (int iperm = 0; iperm < nPerm; iperm++)   //sum over permutation
	{
	A2x = PM[iperm].transpose() * A2_0x * PM[iperm];
	B2x = PM[iperm].transpose() * B2_0x * PM[iperm];
        s2x = PM[iperm].transpose() * s2_0x;
        
        A2y = PM[iperm].transpose() * A2_0y * PM[iperm];
	B2y = PM[iperm].transpose() * B2_0y * PM[iperm];
        s2y = PM[iperm].transpose() * s2_0y;
        
        A2z = PM[iperm].transpose() * A2_0z * PM[iperm];
	B2z = PM[iperm].transpose() * B2_0z * PM[iperm];
        s2z = PM[iperm].transpose() * s2_0z;
        
        
	AAx = A1x + B1x + A2x + B2x;
	InvAAx = AAx.inverse();
        
        AAy = A1y + B1y + A2y + B2y;
	InvAAy= AAy.inverse();
        
        AAz = A1z + B1z + A2z + B2z;
	InvAAz = AAz.inverse();
        
	x1x = sqrt(1.0 / AAx.determinant());
        x1y = sqrt(1.0 / AAy.determinant());
        x1z = sqrt(1.0 / AAz.determinant());
        
	x2x = -0.5*s2x.transpose()*B2x*s2x;
	x3x = -0.5*s1x.transpose()*B1x*s1x;   
        
        x2y = -0.5*s2y.transpose()*B2y*s2y;
	x3y = -0.5*s1y.transpose()*B1y*s1y;  
        
        x2z = -0.5*s2z.transpose()*B2z*s2z;
	x3z = -0.5*s1z.transpose()*B1z*s1z;  
        
        dx = B1x*s1x + B2x*s2x;
        dy = B1y*s1y + B2y*s2y;
        dz = B1z*s1z + B2z*s2z;
        
        x4x = 0.5*dx.transpose()*InvAAx*dx;
        x4y = 0.5*dy.transpose()*InvAAy*dy;
        x4z = 0.5*dz.transpose()*InvAAz*dz;
        
        overlapBx = x1x*exp(x2x+x3x+x4x);
        overlapBy = x1y*exp(x2y+x3y+x4y);
        overlapBz = x1z*exp(x2z+x3z+x4z);


		overlap = overlap + parity[iperm] * stme[iperm*nPairs*nop] * overlapBx*overlapBy*overlapBz;
        }
	return overlap;
}




//==============================================================================================
double MatrixElement::energy(std::vector<MatrixXd> state1, std::vector<MatrixXd> state2)
{
    MatrixXd A1x = state1[0];
	MatrixXd B1x = state1[1];
	VectorXd s1x = state1[2].diagonal();

	MatrixXd A2_0x = state2[0];
	MatrixXd B2_0x = state2[1];
	VectorXd s2_0x = state2[2].diagonal();

	MatrixXd A2x;
	MatrixXd B2x;
        VectorXd s2x;
    
    
        MatrixXd A1y = state1[3];
	MatrixXd B1y = state1[4];
	VectorXd s1y = state1[5].diagonal();

	MatrixXd A2_0y = state2[3];
	MatrixXd B2_0y = state2[4];
	VectorXd s2_0y = state2[5].diagonal();

	MatrixXd A2y;
	MatrixXd B2y;
        VectorXd s2y;
    
    
        MatrixXd A1z = state1[6];
	MatrixXd B1z = state1[7];
	VectorXd s1z = state1[8].diagonal();

	MatrixXd A2_0z = state2[6];
	MatrixXd B2_0z = state2[7];
	VectorXd s2_0z = state2[8].diagonal();
    
        MatrixXd A2z;
	MatrixXd B2z;
        VectorXd s2z;

	MatrixXd AAx;
	MatrixXd InvAAx;
    
    	
        MatrixXd AAy;
 	MatrixXd InvAAy;
    
    	
        MatrixXd AAz;
	MatrixXd InvAAz;
  
  double PotEnergy = 0, KinEnergy = 0, MagneticEnergy = 0, HarmonicEnergy=0, PotEnergy3B = 0;


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


       long double KinEnergyBx, KinEnergyBy, KinEnergyBz, PotEnergyBx, PotEnergyBy, PotEnergyBz;
//===========magnetic===================
        MatrixXd Magneticx;
	MatrixXd Magneticy;
        MatrixXd Magneticz;
        double Mx, My, Mz;

//======================================
       
        VectorXd dx;
        VectorXd dy;
        VectorXd dz;
    
        
        VectorXd yyx;
        VectorXd yyy;
        VectorXd yyz;

	    long double x1x, x1y, x1z, x2x, x2y, x2z, x3x, x3y, x3z, x4x, x4y, x4z, x9x, x9y, x9z, x10x, x10y, x10z;
	    long double v3bx, v3by, v3bz;
        long double y1x, y1y, y1z, y2x, y2y, y2z, y3x, y3y, y3z;
        long double z1x, z1y, z1z, z2, z3, z4, z5;
        long double sx, sy, sz, rx, ry, rz;

                //===========3-body===================
                double PotEnergy3BP = 0;
                MatrixXd BBx(2,2);
                MatrixXd BBBx(2,2);
                MatrixXd BBy(2,2);
                MatrixXd BBBy(2,2);
                MatrixXd BBz(2,2);
                MatrixXd BBBz(2,2);
                VectorXd ex(2);
                VectorXd ey(2);
                VectorXd ez(2);
                double xx1x, xx1y, xx1z, xx2x, xx2y, xx2z, xx7x, xx7y, xx7z;
                //====================================


  
	for (int iperm = 0; iperm < nPerm; iperm++)   //sum over permutation
	{
        A2x = PM[iperm].transpose() * A2_0x * PM[iperm];
	B2x = PM[iperm].transpose() * B2_0x * PM[iperm];
        s2x = PM[iperm].transpose() * s2_0x;
        
        A2y = PM[iperm].transpose() * A2_0y * PM[iperm];
        B2y = PM[iperm].transpose() * B2_0y * PM[iperm];
        s2y = PM[iperm].transpose() * s2_0y;
        
        A2z = PM[iperm].transpose() * A2_0z * PM[iperm];
	B2z = PM[iperm].transpose() * B2_0z * PM[iperm];
        s2z = PM[iperm].transpose() * s2_0z;
        
        
	AAx = A1x + B1x + A2x + B2x;
	InvAAx = AAx.inverse();
        
        AAy = A1y + B1y + A2y + B2y;
	InvAAy= AAy.inverse();
        
        AAz = A1z + B1z + A2z + B2z;
	InvAAz = AAz.inverse();
        
	x1x = sqrt(1.0 / AAx.determinant());
        x1y = sqrt(1.0 / AAy.determinant());
        x1z = sqrt(1.0 / AAz.determinant());
        
	x2x = -0.5*s2x.transpose()*B2x*s2x;
	x3x = -0.5*s1x.transpose()*B1x*s1x;   
        
        x2y = -0.5*s2y.transpose()*B2y*s2y;
	x3y = -0.5*s1y.transpose()*B1y*s1y;  
        
        x2z = -0.5*s2z.transpose()*B2z*s2z;
	x3z = -0.5*s1z.transpose()*B1z*s1z;  

              dx = B1x*s1x + B2x*s2x;
              dy = B1y*s1y + B2y*s2y;
              dz = B1z*s1z + B2z*s2z;
              
              x4x = 0.5*dx.transpose()*InvAAx*dx;
              x4y = 0.5*dy.transpose()*InvAAy*dy;
              x4z = 0.5*dz.transpose()*InvAAz*dz;

                      y3x = x1x*exp(x2x+x3x+x4x);
                      y3y = x1y*exp(x2y+x3y+x4y);
                      y3z = x1z*exp(x2z+x3z+x4z);
        
        
//2B potential enegry=======================================
	  for (int ipt = 0; ipt < npt; ipt++)
      {
	    for (int iop = 0; iop < nop; iop++)
        {
	      for (int ipair = 0; ipair < nPairs; ipair++)   //sum over pairs
		  {
		     sx = TBP[ipair].transpose() * InvAAx * TBP[ipair];
                     sy = TBP[ipair].transpose() * InvAAy * TBP[ipair];
                     sz = TBP[ipair].transpose() * InvAAz * TBP[ipair];

                     x9x = sqrt(1.0 / (2.0*apot(iop, ipt) *sx + 1));
                     x9y = sqrt(1.0 / (2.0*apot(iop, ipt) *sy + 1));
                     x9z = sqrt(1.0 / (2.0*apot(iop, ipt) *sz + 1));
                             

              rx = TBP[ipair].transpose() * InvAAx * dx;
              ry = TBP[ipair].transpose() * InvAAy * dy;
              rz = TBP[ipair].transpose() * InvAAz * dz;

              x10x = -apot(iop, ipt)*rx*rx *(1.0/(2.0*apot(iop, ipt)*sx + 1));
              x10y = -apot(iop, ipt)*ry*ry *(1.0/(2.0*apot(iop, ipt)*sy + 1));
              x10z = -apot(iop, ipt)*rz*rz *(1.0/(2.0*apot(iop, ipt)*sz + 1));
                                                 
              PotEnergyBx = y3x*x9x*exp(x10x);
              PotEnergyBy = y3y*x9y*exp(x10y);
              PotEnergyBz = y3z*x9z*exp(x10z);

		      
		      PotEnergy = PotEnergy + parity[iperm] *stme[(iperm*nPairs+ipair)*nop+iop]*vpot(iop, ipt)*PotEnergyBx*PotEnergyBy*PotEnergyBz;
		  }
	    }
	  }
//========================================================
//==========kinetic energy================================

                      y1x = ((A1x + B1x) * InvAAx * (A2x + B2x)).trace();
                      y1y = ((A1y + B1y) * InvAAy * (A2y + B2y)).trace();
                      y1z = ((A1z + B1z) * InvAAz * (A2z + B2z)).trace();

                      yyx = (A2x+B2x)*InvAAx*(B1x*s1x)-(A1x+B1x)*InvAAx*(B2x*s2x);
                      yyy = (A2y+B2y)*InvAAy*(B1y*s1y)-(A1y+B1y)*InvAAy*(B2y*s2y);
                      yyz = (A2z+B2z)*InvAAz*(B1z*s1z)-(A1z+B1z)*InvAAz*(B2z*s2z);

                      y2x = yyx.transpose()*yyx;
                      y2y = yyy.transpose()*yyy;
                      y2z = yyz.transpose()*yyz;

                      KinEnergyBx = y1x-y2x;
                      KinEnergyBy = y1y-y2y;
                      KinEnergyBz = y1z-y2z;
                      
                      KinEnergy = KinEnergy + parity[iperm] *stme[iperm*nPairs*nop] *(KinEnergyBx + KinEnergyBy + KinEnergyBz)*y3x*y3y*y3z;

//===================================================================
//====================Magnetic energy================================

Magneticx=(InvAAx*dx).asDiagonal();
Mx=(Magneticx*Magneticx).trace()+InvAAx.trace();

Magneticy=(InvAAy*dy).asDiagonal();
My=(Magneticy*Magneticy).trace()+InvAAy.trace();

Magneticz=(InvAAz*dz).asDiagonal();
Mz=(Magneticz*Magneticz).trace()+InvAAz.trace();

MagneticEnergy = MagneticEnergy + parity[iperm] *stme[iperm*nPairs*nop]*y3x*y3y*y3z*My;
HarmonicEnergy = HarmonicEnergy + parity[iperm] *stme[iperm*nPairs*nop]*y3x*y3y*y3z*(Mx+My+Mz);
//===================================================================

if(npar>2)
{
		  PotEnergy3BP = 0;
		  for (int i = 0; i < npar; i++)
          {
		  for (int j = i + 1; j < npar; j++)
          {
		  for (int k = j + 1; k < npar; k++)
          {
		    for (int cyc = 0; cyc < 3; cyc++)
            {
		      if (cyc == 0) {i1 = i; j1 = j; k1 = k;}
		      if (cyc == 1) {i1 = j; j1 = k; k1 = i;}
		      if (cyc == 2) {i1 = k; j1 = i; k1 = j;}
		      Cik = v.SingleParticle(i1, k1);
		      Cjk = v.SingleParticle(j1, k1);

		      Bx(0,0) = Cik.transpose()*InvAAx*Cik;
		      Bx(0,1) = Cik.transpose()*InvAAx*Cjk;
		      Bx(1,0) = Cjk.transpose()*InvAAx*Cik;
		      Bx(1,1) = Cjk.transpose()*InvAAx*Cjk;

		      By(0,0) = Cik.transpose()*InvAAy*Cik;
		      By(0,1) = Cik.transpose()*InvAAy*Cjk;
		      By(1,0) = Cjk.transpose()*InvAAy*Cik;
		      By(1,1) = Cjk.transpose()*InvAAy*Cjk;

		      Bz(0,0) = Cik.transpose()*InvAAz*Cik;
		      Bz(0,1) = Cik.transpose()*InvAAz*Cjk;
		      Bz(1,0) = Cjk.transpose()*InvAAz*Cik;
		      Bz(1,1) = Cjk.transpose()*InvAAz*Cjk;

		      BBx = I+2.0*apot3b*Bx;
                      BBy = I+2.0*apot3b*By;
                      BBz = I+2.0*apot3b*Bz;

		      BBBx = Bx.inverse()*(I-BBx.inverse());
                      BBBy = By.inverse()*(I-BBy.inverse());
                      BBBz = Bz.inverse()*(I-BBz.inverse());

		      xx1x = sqrt(1.0 / BBx.determinant());
                      xx1y = sqrt(1.0 / BBy.determinant());
                      xx1z = sqrt(1.0 / BBz.determinant());
                             

              dx = B1x*s1x + B2x*s2x;
              dy = B1y*s1y + B2y*s2y;
              dz = B1z*s1z + B2z*s2z;

              xx2x = 0.5*dx.transpose()*InvAAx*dx;
              xx2y = 0.5*dy.transpose()*InvAAy*dy;
              xx2z = 0.5*dz.transpose()*InvAAz*dz;
  
	      ex(0)=dx.transpose()*InvAAx*Cik;
	      ex(1)=dx.transpose()*InvAAx*Cjk;

	      ey(0)=dy.transpose()*InvAAy*Cik;
	      ey(1)=dy.transpose()*InvAAy*Cjk;

	      ez(0)=dz.transpose()*InvAAz*Cik;
	      ez(1)=dz.transpose()*InvAAz*Cjk;

	      xx7x=-0.5*ex.transpose()*BBBx*ex;
              xx7y=-0.5*ey.transpose()*BBBy*ey;
              xx7z=-0.5*ez.transpose()*BBBz*ez;

	      v3bx = x1x*exp(x2x+x3x)*xx1x*exp(xx2x+xx7x);
              v3by = x1y*exp(x2y+x3y)*xx1y*exp(xx2y+xx7y);
              v3bz = x1z*exp(x2z+x3z)*xx1z*exp(xx2z+xx7z);


	      PotEnergy3BP = PotEnergy3BP + v3bx*v3by*v3bz;
		    }
		  }
       }
    }

                  PotEnergy3B = PotEnergy3B + parity[iperm] * stme[iperm*nPairs*nop] * PotEnergy3BP;

}  //  end if(npar>2) 

               
	}   // end sum over permutation
 
        KinEnergy=0.5*h2m*KinEnergy;
        PotEnergy3B = vpot3b * PotEnergy3B;

      MagneticEnergy = 0.5 * h2m * eB * eB * MagneticEnergy; 
      HarmonicEnergy = 0.5 * h2m * momega * momega * HarmonicEnergy;

//cout<<endl<<KinEnergy<<endl;
//cout<<endl<<PotEnergy<<endl;
//cout<<endl<<PotEnergy3B<<endl;
//cout<<endl<<MagneticEnergy<<endl;
//cout<<endl<<HarmonicEnergy<<endl<<endl<<endl;
        
  return PotEnergy+KinEnergy+PotEnergy3B+MagneticEnergy+HarmonicEnergy;
}

//=============================================================================








vector<VectorXd> SumOverCells(int dmax, int npar)
{
	int Nconf = pow(2 * dmax + 1, npar); //number of configuration
	vector<VectorXd> d(Nconf);
	VectorXd b(npar);
	for (int i = 0; i < npar; i++)  b(i) = -dmax; 

	for (int i = 0; i < Nconf; i++)
	{
		d[i] = b;
		b(0)++;
		for (int k = 0; k < npar - 1; k++)
		{
			if (b(k) > dmax)
			{
				b(k) = -dmax;
				b(k + 1)++;
			}
		}
	}
	return d;
}

//=============================================================================

vector<VectorXd> TwoBodyPairs(int npar)
{
	vector<VectorXd> Cij(npar*(npar - 1) / 2);
	CoordinatsTransformation v(npar);
	int ipair = 0;
	for (int i = 0; i < npar; i++)	{
		for (int j = i + 1; j < npar; j++){
			Cij[ipair] = v.SingleParticle(i, j);
			ipair++;
		}
	}
	return Cij;
}
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
        int keypr = 1;

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
















/*
//harmonic=====================================

//=============================================================================
//double MatrixElement::overlap(std::vector<MatrixXd> state1, std::vector<MatrixXd> state2)
{
	MatrixXd A1 = state1[0];
	MatrixXd B1 = state1[1];
	MatrixXd A2 = state2[0];
	MatrixXd B2 = state2[1];

	return pow((A1 + B1 + A2 + B2).determinant(),-1.5);
}

//=============================================================================
//double MatrixElement::energy(std::vector<MatrixXd> state1, std::vector<MatrixXd> state2)
{
  MatrixXd A1 = state1[0];
  MatrixXd B1 = state1[1];
  MatrixXd A2 = state2[0];
  MatrixXd B2 = state2[1];

  
  MatrixXd AA = A1 + B1 + A2 + B2;
  MatrixXd invAA = AA.inverse(); 
  double overlap = pow(AA.determinant(), -1.5); 
  MatrixXd TT = (A1 + B1) * invAA * (A2 + B2);
  double PotEnergy = 0, KinEnergy = 0;
  Vector2d c(1, -1);
  KinEnergy = 0.5 * h2m * TT.trace();
  PotEnergy = (0.25/h2m) * c.transpose() * invAA * c;

  return overlap*(PotEnergy+KinEnergy);
}

// end harmoinic====================================

*/


MatrixElement::~MatrixElement()
{
}

