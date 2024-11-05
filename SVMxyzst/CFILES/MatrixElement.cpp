#include "MatrixElement.h"
#include "Input.h"
#include "Operators.h"
#include "CoordinatsTransformation.h"
#include "Permutation.h"
#include "BasisState.h"
#include <cmath>
#include <Eigen> 
#include <vector>
using namespace Eigen;
using namespace std;

vector<VectorXd> SumOverCells(int dmax, int npar);
vector<VectorXd> TwoBodyPairs(int npar);
//=============================================================================
MatrixElement::MatrixElement(Input &input)
{
    cout << "\t Initialize MatrixElement\n";
    npar = input.npar;
    nts_states = input.nts_states;
    h2m = input.h2m;
    eB=input.eB;
    eBspin=input.eBspin;
    //momega=input.momega;
    momega=0.0;
//=======================
    KinEnergy_cof      = 0.5 * h2m;
    MagneticEnergy_cof = 0.25 * pow(npar, -1)* 0.5 * h2m * eB * eB; 
    harmonic_cof       = 0.25 * pow(npar, -1)* 0.5 * h2m * momega * momega;
    MagneticSpin_cof   = -0.5* h2m * eBspin;
//=======================
    isospin   = input.isospin;
    spin      = input.spin;
//=======================
    dmax   = input.dmax;
    nop    = input.nop;
    npt    = input.npt;
    ibf    = input.bosefermi;  // ibf=2 for bosons.  ibf=1 for fermions
    apot3b = input.apot3b;
    vpot3b = input.vpot3b;
    nPairs = npar*(npar - 1) / 2;
    nPerm  = factorial(npar);
    nConf  = pow(2 * dmax + 1, npar); //number of cells configurations

    Permutation P(npar); 
    PM                = P.perm_matrix;
    PV                = P.permutation;
    if(ibf==1) parity = P.parity; //fermions
    if (ibf == 2) for (int iperm = 0; iperm < nPerm; iperm++) parity.push_back(1); //bosons
	
    mass = MatrixXd::Zero(npar,npar);
    for (int imas = 0; imas < npar; imas++) mass(imas,imas)=input.mass[imas];

    TBP = TwoBodyPairs(npar);

    PrepareSpinIsospinME(input);
    PreparePotential(input);
    //PrepareMagneticSpinME(input);
    //PrepareMagneticMePair(input);
                     
}
//=============================================================================
double MatrixElement::overlap(BasisState &state1, BasisState &state2){
    MatrixXd A1x  = state1.Ax;
    MatrixXd A1y  = state1.Ay;
    MatrixXd A1z  = state1.Az;
    int      its1 = state1.ts;

    MatrixXd A2_0x = state2.Ax;
    MatrixXd A2_0y = state2.Ay;
    MatrixXd A2_0z = state2.Az;
    int      its2  = state2.ts;

    MatrixXd A2x;
    MatrixXd A2y;
    MatrixXd A2z;

    double detx, dety, detz;
    long double xx;

    double overlap=0;

    int ijts = its1*nts_states+its2;
    for (int iperm = 0; iperm < nPerm; iperm++){   //sum over permutation
	
	A2x = PM[iperm].transpose() * A2_0x * PM[iperm];
	A2y = PM[iperm].transpose() * A2_0y * PM[iperm];
	A2z = PM[iperm].transpose() * A2_0z * PM[iperm];
	
	detx=(A1x+A2x).determinant();
	dety=(A1y+A2y).determinant();
	detz=(A1z+A2z).determinant();
      
	xx=1.0/sqrt(detx*dety*detz);
//	     overlap = overlap + parity[iperm]* stme[iperm*nPairs*nop] *x;
	overlap = overlap + parity[iperm]* stmeop[ijts].perm[iperm].op2b[0].me(0) *xx;
    }
    return overlap;
}

//==============================================================================================
double MatrixElement::energy(BasisState &state1, BasisState &state2){    
    MatrixXd A1x  = state1.Ax;
    MatrixXd A1y  = state1.Ay;
    MatrixXd A1z  = state1.Az;
    int      its1 = state1.ts;

    //cout<<endl<<A1x<<endl<<endl<<A1x.determinant()<<endl<<endl;

    MatrixXd A2_0x = state2.Ax;
    MatrixXd A2_0y = state2.Ay;
    MatrixXd A2_0z = state2.Az;
    int      its2  = state2.ts;

    MatrixXd A2x, A2y, A2z;
    MatrixXd InvAAx, InvAAy, InvAAz;
    MatrixXd InvAAmx, InvAAmy, InvAAmz;
    MatrixXd TTx, TTy, TTz;

    MatrixXd qInvAAx = MatrixXd::Zero(npar, npar);
    MatrixXd qInvAAy = MatrixXd::Zero(npar, npar);
    MatrixXd qInvAAz = MatrixXd::Zero(npar, npar);

    double detx, dety, detz;
    double sx,sy,sz, xme; 
    long double x, y, z;
    double mi, mj;
  
    double PotEnergy = 0, KinEnergy = 0, MagneticEnergy=0;
    double MagneticEnergyT=0,MagneticSpin=0, harmonic=0, PotEnergy3B = 0;
    double PotEnergy3BP = 0;

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
    int ijts = its1*nts_states+its2;
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
      
      	x=1/sqrt(detx*dety*detz);

      	// Kinetic Energy 
      	TTx = A1x * InvAAx * A2x;
      	TTy = A1y * InvAAy * A2y;
      	TTz = A1z * InvAAz * A2z;
//     	KinEnergy = KinEnergy + parity[iperm] * stme[iperm*nPairs*nop] *(TTx.trace()+TTy.trace()+TTz.trace())*x;
      	KinEnergy = KinEnergy + parity[iperm] * stmeop[ijts].perm[iperm].op2b[0].me(0)
	    *(TTx.trace()+TTy.trace()+TTz.trace())*x;
	//==================================================================================
/*      	// Magnetic Energy (single-particle version)
      	for (int ipar = 0 ; ipar < npar ; ipar++){
	    MagneticEnergyT = MagneticEnergyT + parity[iperm] * 
		magnetic_charge_me[iperm*npar+ipar] * (InvAAx(ipar,ipar)+InvAAy(ipar,ipar)) * x;
      	}
	//==================================================================================
	
	//=============MagneticEnergy=======================================================  
	int ipair = 0;			
	for (int i = 0; i < npar; i++) {
	    for (int j = i+1; j < npar; j++) {
		sx =    magnetic_me_pair[iperm*nPairs+ipair](0) * InvAAx(i,i) 
		    +   magnetic_me_pair[iperm*nPairs+ipair](1)*InvAAx(j,j)
		    - 2*magnetic_me_pair[iperm*nPairs+ipair](2)*InvAAx(i,j);
		sy =    magnetic_me_pair[iperm*nPairs+ipair](0) * InvAAy(i,i) 
		    +   magnetic_me_pair[iperm*nPairs+ipair](1)*InvAAy(j,j)
		    - 2*magnetic_me_pair[iperm*nPairs+ipair](2)*InvAAy(i,j);
		sz =    magnetic_me_pair[iperm*nPairs+ipair](0) * InvAAz(i,i) 
		    +   magnetic_me_pair[iperm*nPairs+ipair](1)*InvAAz(j,j)
		    - 2*magnetic_me_pair[iperm*nPairs+ipair](2)*InvAAz(i,j);
		ipair++;
		
		MagneticEnergy = MagneticEnergy + parity[iperm]*(sx+sy)*x;
		harmonic = harmonic + parity[iperm] * sz * x;
		
		//MagneticEnergy = MagneticEnergy + parity[iperm]*magnetic_me_pair[iperm*nPairs]*(sx+sy)*x;
		//harmonic = harmonic + parity[iperm] * magnetic_me_pair[iperm*nPairs+ipair] *x* sz;
	    }
	    
	} 
	
	//====================MagneticSpinEnergy============================================
		
	for (int ipar = 0; ipar < npar; ipar++){   //sum over single particles
	    MagneticSpin = MagneticSpin + parity[iperm] * magnetic_spin_me[iperm*npar+ipar]*x;    
	}
		
		//=========================2-body===================================================
        for (int ipair = 0; ipair < nPairs; ipair++){   //sum over pairs
	    sx = TBP[ipair].transpose() * InvAAx * TBP[ipair];
	    sy = TBP[ipair].transpose() * InvAAy * TBP[ipair];
	    sz = TBP[ipair].transpose() * InvAAz * TBP[ipair];
	    for (int iop = 0; iop < nop; iop++){
//		xme = parity[iperm] * stme[(iperm*nPairs+ipair)*nop+iop];
		xme = parity[iperm] * stmeop[ijts].perm[iperm].op2b[iop].me(ipair);
		for (int ipt = 0; ipt < npt; ipt++){
		    y=1/sqrt( (2.0*apot(iop, ipt)*sx + 1)
			     *(2.0*apot(iop, ipt)*sy + 1)
			     *(2.0*apot(iop, ipt)*sz + 1));
		    PotEnergy = PotEnergy + xme * vpot(iop, ipt)*x*y;
		}
	    }
	}
		//===============================================================
*/
	if (npar>2){
	    PotEnergy3BP = 0;
	    for (int i = 0; i < npar; i++){
		for (int j = i + 1; j < npar; j++){
		    for (int k = j + 1; k < npar; k++){
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

                        	z=1/sqrt(  (2.0*apot3b*Bx + I).determinant()
					 * (2.0*apot3b*By + I).determinant()
					  *(2.0*apot3b*Bz + I).determinant());                          
				PotEnergy3BP = PotEnergy3BP + z*x;
			    }
		    }
		}
	    }
//	  		PotEnergy3B = PotEnergy3B + parity[iperm] * stme[iperm*nPairs*nop] *PotEnergy3BP;
	    PotEnergy3B = PotEnergy3B + parity[iperm] * stmeop[ijts].perm[iperm].op2b[0].me(0)
		*PotEnergy3BP;
	}
	  	//====================================
    }    

    PotEnergy3B     = vpot3b             * PotEnergy3B;

    KinEnergy       = KinEnergy_cof      * KinEnergy;
    //MagneticEnergy  = MagneticEnergy_cof * MagneticEnergy;
    //MagneticEnergyT = MagneticEnergy_cof * MagneticEnergyT;
    //harmonic        = harmonic_cof       * harmonic;
    //MagneticSpin    = MagneticSpin_cof   * MagneticSpin;

    return PotEnergy + KinEnergy + PotEnergy3B; //+ MagneticEnergyT + harmonic + MagneticSpin
}
//=============================================================================

//=============================================================================
vector<VectorXd> TwoBodyPairs(int npar){
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
/*
vector<VectorXd> MatrixElement::ChargeVec()
{
	vector<VectorXd> charge_vec(nPerm);
	for (int iperm = 0; iperm < nPerm; iperm++)
          {
              for (int ipar = 0; ipar < npar; ipar++)
                {
                      charge_vec[iperm](ipar) = One_magneticme(ipar, PV[iperm]); 
                }                          
	  }
	return charge_vec;
}
*/
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
void MatrixElement::PrepareSpinIsospinME(Input &input){
    Operators operators(input);

//  Define the ST matrix element structure
    nop1b=1;
    nop2b=nop;
    int nst = input.ts_states.size();
    stmeop.resize(nst*nst);
    for (int ij = 0; ij < nst*nst; ij++){
	stmeop[ij].perm.resize(nPerm);
	for (int iperm = 0; iperm < nPerm; iperm++){
	    stmeop[ij].perm[iperm].op1b.resize(nop1b);
	    stmeop[ij].perm[iperm].op2b.resize(nop2b);
	    for (int iop1b = 0; iop1b < nop1b; iop1b++){
		stmeop[ij].perm[iperm].op1b[iop1b].me.resize(npar);
	    }
	    for (int iop2b = 0; iop2b < nop2b; iop2b++){
		stmeop[ij].perm[iperm].op2b[iop2b].me.resize(nPairs);
	    }
	}
    }
// Calculate the 2-body matrix elements
    for (int istl = 0; istl < nst; istl++){
    for (int istr = 0; istr < nst; istr++){
	SpinIsospinState st_l = input.ts_states[istl];
	SpinIsospinState st_r = input.ts_states[istr];
	int ij = istl*nst+istr;
	for (int iperm = 0; iperm < nPerm; iperm++){
	    stmeop[ij].perm[iperm].op2b = operators.OST_2bme(st_l,st_r,PV[iperm]);
	}
    }
    }
}
//=============================================================================
//   PrepareMagneticSpinME
void MatrixElement::PrepareMagneticSpinME(Input &input){
    magnetic_spin_me.resize(npar*nPerm);
    magnetic_charge_me.resize(npar*nPerm);
    for (int ipar = 0; ipar < npar; ipar++){
	for (int iperm = 0; iperm < nPerm; iperm++){
            magnetic_charge_me[iperm*npar+ipar] = One_magnetic_me(ipar, PV[iperm], 1);
            magnetic_spin_me[iperm*npar+ipar] = One_magnetic_me(ipar, PV[iperm], 2);
	}
    }
}
//=======================================
double MatrixElement::One_magnetic_me(int ipar, VectorXi Perm, int i_spin_coupling){
    double x = 0;
    for (int icmp = 0; icmp < spin.ncmp; icmp++){
	for (int jcmp = 0; jcmp < spin.ncmp; jcmp++){
	    x = x + spin.coef[icmp] * spin.coef[jcmp] 
		* SpinOp(spin.cmp[icmp].sz,spin.cmp[jcmp].sz, npar, ipar, Perm, i_spin_coupling);
	}
    }

    double y = 0;
    for (int icmp = 0; icmp < isospin.ncmp; icmp++){
	for (int jcmp = 0; jcmp < isospin.ncmp; jcmp++){
	    y = y + isospin.coef[icmp] * isospin.coef[jcmp]
		* IsospinOp(isospin.cmp[icmp].tz, isospin.cmp[jcmp].tz, 
			    npar, ipar, Perm, i_spin_coupling);
	}
    }
    return x*y;
}
//=============================================================================
double MatrixElement::SpinOp(std::vector<int> sz1, std::vector<int> sz2, int npar, int ipar, 
			     VectorXi Perm, int i_spin_coupling){
    for (int ip = 0; ip < npar; ip++) if (sz1[ip] != sz2[Perm(ip)]) return 0.0;
    if(i_spin_coupling==1){
	return 1.0;
    }
    if(i_spin_coupling==2){
	//	1 for spin up,  2 for spin down
	if (sz2[Perm(ipar)] == 1)  return 0.5;
	if (sz2[Perm(ipar)] == 2)  return -0.5;
    }
}
//=============================================================================
double MatrixElement::IsospinOp(std::vector<int> tz1, std::vector<int> tz2, int npar, int ipar, 
		 VectorXi Perm, int i_spin_coupling)
{
	for (int ip = 0; ip < npar; ip++) if (tz1[ip] != tz2[Perm(ip)]) return 0.0;
	if(i_spin_coupling==1)
 	{
    	//	1 for proton q=1,   2 for netron q=0
		if (tz2[Perm(ipar)] == 1)  return 1.0;  //q1
		if (tz2[Perm(ipar)] == 2)  return 0.0;  //q2
 	}
	if(i_spin_coupling==2)
 	{
    	//	1 for proton g_p=5.586,   2 for netron g_n=-3.826
		if (tz2[Perm(ipar)] == 1)  return 5.586;
		if (tz2[Perm(ipar)] == 2)  return -3.826;
	}
}
//==============================EndPrepareMagneticME==============================================

//=============================================================================
int MatrixElement::factorial(int npar)
{
	int N;
	if (npar <= 1) return 1;
	N = npar * factorial(npar - 1);
	return N;
}
//=============================================================================
int MatrixElement::sign(double x)
{
        int N;
	if (x >= 0) N = 1;
	if (x < 0) N = -1;
	return N;
}
//=============================================================================
void MatrixElement::PrepareMagneticMePair(Input &input){
    magnetic_me_pair.resize(nPairs*nPerm);
    int ipair = -1;
    for (int ip = 0; ip < npar; ip++){
	for (int jp = ip + 1; jp < npar; jp++){
	    ipair++;
	    for (int iperm = 0; iperm < nPerm; iperm++){   
		magnetic_me_pair[iperm*nPairs+ipair]=VectorXd::Zero(3); 
		magnetic_me_pair[iperm*nPairs+ipair](0) = One_magnetic_me_pair(ip, jp, PV[iperm], 0);  
		magnetic_me_pair[iperm*nPairs+ipair](1) = One_magnetic_me_pair(ip, jp, PV[iperm], 1); 
		magnetic_me_pair[iperm*nPairs+ipair](2) = One_magnetic_me_pair(ip, jp, PV[iperm], 2);              	      	}
	}
    }
}
//=============================================================================
double MatrixElement::One_magnetic_me_pair(int i, int j, VectorXi Perm, int index)
{
   	double x = 0;
	for (int icmp = 0; icmp < spin.ncmp; icmp++){
		for (int jcmp = 0; jcmp < spin.ncmp; jcmp++){
	    	x = x + spin.coef[icmp] * spin.coef[jcmp] 
	        	* SpinOp_pair(spin.cmp[icmp].sz,spin.cmp[jcmp].sz, npar, i, j, Perm);
	  	}
	}

	double  y = 0;
	for (int icmp = 0; icmp < isospin.ncmp; icmp++){
		for (int jcmp = 0; jcmp < isospin.ncmp; jcmp++){
	    	y = y + isospin.coef[icmp] * isospin.coef[jcmp]
                  * IsospinOp_pair(isospin.cmp[icmp].tz, isospin.cmp[jcmp].tz, npar, i, j, Perm, index);
	  	}
	}
	return x*y;
}
//=============================================================================
double MatrixElement::SpinOp_pair(std::vector<int> sz1, std::vector<int> sz2, int npar, int i, int j, 
				  VectorXi Perm)
{
	for (int ipar = 0; ipar < npar; ipar++)
	{
		if ((ipar != i) && (ipar != j))
		{
			if (sz1[ipar] != sz2[Perm(ipar)]) return 0;
		}
	}

    double x = 0;
	if ((sz1[i] == sz2[Perm(i)]) && (sz1[j] == sz2[Perm(j)])) x = 1.0;
	return x;
}
//=============================================================================
double MatrixElement::IsospinOp_pair(std::vector<int> tz1, std::vector<int> tz2, int npar, 
				     int i, int j, VectorXi Perm, int index){
    for (int ipar = 0; ipar < npar; ipar++){
	       if ((ipar != i) && (ipar != j))
		{
			if (tz1[ipar] != tz2[Perm(ipar)]) return 0;
		}
	}
    double x = 0;
    double qn = 0.00;
    double qp = 1;
        //1 for proton qp=1,   2 for netron qn=0.001
	if ((tz1[i] == tz2[Perm(i)]) && (tz1[j] == tz2[Perm(j)]))
           {
              if(index= 0)
                {
                   if ((tz2[Perm(i)] == 1) &&(tz2[Perm(j)] == 1))  x = qp*qp;
                   if ((tz2[Perm(i)] == 1) &&(tz2[Perm(j)] == 2))  x = qp*qp;
                   if ((tz2[Perm(i)] == 2) &&(tz2[Perm(j)] == 1))  x = qn*qn;
                   if ((tz2[Perm(i)] == 2) &&(tz2[Perm(j)] == 2))  x = qn*qn;
                }
              if(index= 1)
                {
                   if ((tz2[Perm(i)] == 1) &&(tz2[Perm(j)] == 1))  x = qp*qp;
                   if ((tz2[Perm(i)] == 1) &&(tz2[Perm(j)] == 2))  x = qn*qn;
                   if ((tz2[Perm(i)] == 2) &&(tz2[Perm(j)] == 1))  x = qp*qp;
                   if ((tz2[Perm(i)] == 2) &&(tz2[Perm(j)] == 2))  x = qn*qn;
                }
              if(index= 2)
                {
                   if ((tz2[Perm(i)] == 1) &&(tz2[Perm(j)] == 1))  x = qp*qp;
                   if ((tz2[Perm(i)] == 1) &&(tz2[Perm(j)] == 2))  x = qp*qn;
                   if ((tz2[Perm(i)] == 2) &&(tz2[Perm(j)] == 1))  x = qp*qn;
                   if ((tz2[Perm(i)] == 2) &&(tz2[Perm(j)] == 2))  x = qn*qn;
                }
           }

	return x;
}


/*

                   if ((tz2[Perm(i)] == 1) &&(tz2[Perm(j)] == 1))  x = 1;
                   if ((tz2[Perm(i)] == 1) &&(tz2[Perm(j)] == 2))  x = 1;
                   if ((tz2[Perm(i)] == 2) &&(tz2[Perm(j)] == 1))  x = q;
                   if ((tz2[Perm(i)] == 2) &&(tz2[Perm(j)] == 2))  x = q;
                }
              if(index= 1)
                {
                   if ((tz2[Perm(i)] == 1) &&(tz2[Perm(j)] == 1))  x = 1;
                   if ((tz2[Perm(i)] == 1) &&(tz2[Perm(j)] == 2))  x = q;
                   if ((tz2[Perm(i)] == 2) &&(tz2[Perm(j)] == 1))  x = 1;
                   if ((tz2[Perm(i)] == 2) &&(tz2[Perm(j)] == 2))  x = q;
                }
              if(index= 2)
                {
                   if ((tz2[Perm(i)] == 1) &&(tz2[Perm(j)] == 1))  x = 1;
                   if ((tz2[Perm(i)] == 1) &&(tz2[Perm(j)] == 2))  x = q;
                   if ((tz2[Perm(i)] == 2) &&(tz2[Perm(j)] == 1))  x = q;
                   if ((tz2[Perm(i)] == 2) &&(tz2[Perm(j)] == 2))  x = q*q;
                }

*/
//=============================================================================
MatrixElement::~MatrixElement(){}
//=============================================================================
