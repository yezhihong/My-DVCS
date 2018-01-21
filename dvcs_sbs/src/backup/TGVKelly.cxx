/*--------------------------------------- TGVKelly.cxx --------------------------------------*
 |                                                                                      |
 | Author : H. Moutarde (CEA-Saclay, IRFU/SPhN).                                        |
 | v1.4, July, 25th 2008.                                                               |
 |                                                                                      |
 | Computation of the differential cross sections of the Bethe-Heitler and Deeply       |
 | Virtual Compton Scattering processes (and also of the interference of the two        |
 | processes).                                                                          |
 |                                                                                      |
 | The formulas are derived from a Mathematica package written by P. Guichon            |
 | (CEA-Saclay, IRFU/SPhN) and M. Vanderhaegen (W&M College), see arXiv:0808.????       |
 |                                                                                      | 
 | This library is not standalone, since some Root librairies are necessary.            |
 |                                                                                      |
 | All expressions are written in a class named TGVKelly.                                    | 
 |                                                                                      |
 *--------------------------------------------------------------------------------------*/ 
 
 
 
/*------------------------------- Necessary libraries  ---------------------------------*/
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "TGVKelly.h"
using namespace std;



/*--------------------------------------- Constructors ---------------------------------*/

TGVKelly::TGVKelly() // Default 
{
	if ( NoPrint == kFALSE )
	{
		cout << "TGVKelly : Call of TGVKelly default constructor" << endl;
	} // end if NoPrint
	SetBeamEnergy(5.77);
} // end default TGVKelly



TGVKelly::TGVKelly( Double_t EBeam, Bool_t Valid = kFALSE, Bool_t NoPr = kTRUE )
{
	Validation = Valid;
	NoPrint = NoPr;
	InitBeamEnergy = kFALSE;
	InitKinematics = kFALSE;
	InitMvcs = kFALSE;
	InitExactBHCrossSections = kFALSE;
	InitExactVCSAndInterfCrossSections = kFALSE;
	InitLeadingBHCrossSections = kFALSE;
	InitLeadingVCSAndInterfCrossSections = kFALSE;
	M = 0.93827231;
		
		
	if ( NoPrint == kFALSE )
	{
		cout << "TGVKelly : Call of TGVKelly constructor" << endl;
	} // end if NoPrint
	
	
	if ( Validation == kTRUE && NoPrint == kFALSE )
	{
		cout << "    TGVKelly : Package validation mode : all variables printed in working directory ! " << endl;
		cout << "        TGVKelly : Kinematics printed in CheckKinematics.dat" << endl;
		cout << "        TGVKelly : Form Factors printed in CheckFormFactors.dat" << endl;
		cout << "        TGVKelly : Mvcs printed in CheckMvcs.dat" << endl;
		cout << "        TGVKelly : Exact Urs printed in CheckUrExact.dat" << endl;
		cout << "        TGVKelly : Exact Jem printed in CheckJemExact.dat" << endl;
		cout << "        TGVKelly : Exact SigmaBHs printed in CheckSigmaBHExact.dat" << endl;
		cout << "        TGVKelly : Exact SigmaVCSs printed in CheckSigmaVCSExact.dat" << endl;	
		cout << "        TGVKelly : Exact SigmaIs printed in CheckSigmaIExact.dat" << endl;
		cout << "        TGVKelly : Leading Urs printed in CheckUrLeading.dat" << endl;
		cout << "        TGVKelly : Leading Jem printed in CheckJemLeading.dat" << endl;
		cout << "        TGVKelly : Leading SigmaBHs printed in CheckSigmaBHLeading.dat" << endl;
		cout << "        TGVKelly : Leading SigmaVCSs printed in CheckSigmaVCSLeading.dat" << endl;	
		cout << "        TGVKelly : Leading SigmaIs printed in CheckSigmaILeading.dat" << endl;	
	} // end if Validation & NoPrint
		
	SetBeamEnergy(EBeam);
	
	if ( NoPrint == kFALSE )
	{
		cout << "    TGVKelly : Initialization status : " << endl;
		cout << "        TGVKelly : Kinematics initialization : " << InitKinematics << endl;
		cout << "        TGVKelly : Helicity amplitude initialization : " << InitMvcs << endl;
		cout << "        TGVKelly : Exact Bethe Heitler Cross Section Initialization : " << InitExactBHCrossSections << endl;	
		cout << "        TGVKelly : Exact Virtual Compton Scattering and Interference Cross Sections Initialization : " << InitExactVCSAndInterfCrossSections << endl;
		cout << "        TGVKelly : Leading Bethe Heitler Cross Section Initialization : " << InitLeadingBHCrossSections << endl;	
		cout << "        TGVKelly : Leading Virtual Compton Scattering and Interference Cross Sections Initialization : " << InitLeadingVCSAndInterfCrossSections << endl;
	} // end if NoPrint
} // end TGVKelly



/*-------------------------------------- Destructor ------------------------------------*/

TGVKelly::~TGVKelly()
{
	if ( NoPrint == kFALSE )
	{
		cout << "TGVKelly : Call of TGVKelly desctructor" << endl;
	} // end if NoPrint
} // end ~TGVKelly



/*---------------------------------- Methods' definition -------------------------------*
 |   - SetBeamEnergy(EBeam)                                                             |
 |	 - MakeKinematics()                                                                 |
 |	 - MakeHelicityAmplitudes()                                                         | 
 |   - MakeExactBHCrossSections()                                                       | 
 |   - MakeExactVCSAndInterfCrossSections()                                             |
 |   - MakeLeadingBHCrossSections()                                                     |
 |   - MakeLeadingVCSAndInterfCrossSections()                                           |
 |	 - DdirectDcrossed(phi)                                                             |
 |	 - SqrAmplBH(Q2,xB,t,phi,BeamHeli,TargetPolar,ReH,ImH,ReE,ImE,ReHT,ImHT,ReET,ImET)  |
 |	 - SqrAmplVCS(Q2,xB,t,phi,BeamHeli,TargetPolar,ReH,ImH,ReE,ImE,ReHT,ImHT,ReET,ImET) |
 |	 - SqrAmplInterf(Q2,xB,t,phi,BeamHeli,TargetPolar,BeamCharge,ReH,ImH,ReE,ImE,ReHT,  |
 |           ImHT,ReET,ImET)                                                            |
 *--------------------------------------------------------------------------------------*/
 
 
 
/*--------------------------- Function SetBeamEnergy(EBeam) ----------------------------*
 | Sets the beam energy in the laboratory frame.                                        |
 *--------------------------------------------------------------------------------------*/
 
void TGVKelly::SetBeamEnergy( Double_t EBeam )
{
	if ( InitBeamEnergy == kFALSE)
	{
		ELab = EBeam;
	
				
		// Flag
		
		InitBeamEnergy = kTRUE;
		if ( NoPrint == kFALSE )
		{
			cout << "TGVKelly : Beam Energy Initialization : " << InitBeamEnergy << endl;
			cout << "    TGVKelly : ELab = " << EBeam << " GeV" << endl;		
		} // end if NoPrint
		
	} // end if InitBeamEnergy
	
} // end SetBeamEnergy



/*--------------- Function MakeKinematics(Q2Input,xBInput,tInput) ----------------------*
 | Sets the global variables xB, Q2 and t and computes the 4-vector components of the   |
 | particules involved in the center of mass frame :                                    |
 |   - virtual photon's 4-momentum qCM,                                                 |
 |   - real photon's 4-momentum qpCM,                                                   |
 |   - incoming proton's 4-momentum pCM,                                                |
 |   - outgoing proton's 4-momentum ppCM,                                               |
 |   - thetag angle between the trajectories of real and virtual photons.               |
 | It also evaluates the parameter Omega and the invariant phase space PhaseSpace.      |
 *--------------------------------------------------------------------------------------*/

void TGVKelly::MakeKinematics( Double_t Q2Input, Double_t xBInput, Double_t tInput )
{
	if ( InitBeamEnergy == kFALSE )
	{
 		cout << "TGVKelly : Computation of Kinematics" << endl;
		cout << "    TGVKelly : Beam energy initialization : " << InitKinematics << endl;
		cout << "    TGVKelly : Incomplete initialization !" << endl;
		exit(-1);		
	} // end if InitBeamEnergy
	
	if (InitKinematics == kFALSE)
	{
		Double_t y; // Variable used in the computation of Omega;
		Double_t e = 0.30282211985966434;
		Double_t xBMin, xBMax; // Boundaries on the xB physical region
		
		
		// Q2, xB, t as global variables and derived quantities
		
		Q2 = Q2Input;
		xB = xBInput;
		t = tInput;
		
		Q = TMath::Sqrt(Q2);
		s = TMath::Power(M,2) - TMath::Power(Q,2) + TMath::Power(Q,2)/xB;
		
		
		// Test : realistic kinematic configuration ?
		
		xBMin = 2.*ELab*Q2/(M*(4*TMath::Power(ELab,2)-Q2)); 
		xBMax = Q2/(Q2-TMath::Power(M,2));
// The value of xBMin comes from the requirement of omega to be real, and the value of xBMax expresses the fact that s >= 0.		


		if( !( xBMin < xB && xB < xBMax ) ) 
		{
	    	cout << "TGVKelly : Unrealistic kinematic configuration : xB isn't in the physical region !" << endl;
    		cout << "    TGVKelly : xB = " << xB << endl;
    		cout << "    TGVKelly : xBMin = " << xBMin << endl;
    		cout << "    TGVKelly : xBMax = " << xBMax << endl;
    		exit(-1);
		  }

	
		// Overall initialization
	
		qCM.SetPxPyPzE(0.,0.,0.,0.);
		qpCM.SetPxPyPzE(0.,0.,0.,0.);
		pCM.SetPxPyPzE(0.,0.,0.,0.);
		ppCM.SetPxPyPzE(0.,0.,0.,0.);
	
	
		// Timelike coordinate

		qCM.SetE((TMath::Power(Q,2)*(1 - 2*xB))/(2.*TMath::Sqrt(TMath::Power(M,2) + (TMath::Power(Q,2)*(1 - xB))/xB)*xB));
		qpCM.SetE(-(TMath::Power(Q,2)*(-1 + xB))/(2.*TMath::Sqrt(TMath::Power(M,2) + (TMath::Power(Q,2)*(1 - xB))/xB)*xB));
		pCM.SetE((TMath::Power(Q,2) + 2*TMath::Power(M,2)*xB)/
   (2.*TMath::Sqrt(TMath::Power(M,2) + (TMath::Power(Q,2)*(1 - xB))/xB)*xB));
		ppCM.SetE((2*TMath::Power(M,2) + (TMath::Power(Q,2)*(1 - xB))/xB)/
   (2.*TMath::Sqrt(TMath::Power(M,2) + (TMath::Power(Q,2)*(1 - xB))/xB)));
	
	
		// Spacelike coordinates
	
		qpCM.SetPx(TMath::Sqrt((-(TMath::Power(M,2)*TMath::Power(t,2)*TMath::Power(xB,2)) + 
      TMath::Power(Q,2)*t*xB*(t*(-1 + xB) - 2*TMath::Power(M,2)*xB) + 
      TMath::Power(Q,4)*(t*(-1 + xB) - TMath::Power(M,2)*TMath::Power(xB,2)))/
    (TMath::Power(Q,4) + 4*TMath::Power(M,2)*TMath::Power(Q,2)*TMath::Power(xB,2))));
		qCM.SetPz(TMath::Sqrt((TMath::Power(Q,4) + 4*TMath::Power(M,2)*TMath::Power(Q,2)*TMath::Power(xB,2))/
     (xB*(TMath::Power(Q,2)*(1 - xB) + TMath::Power(M,2)*xB)))/2.);
		qpCM.SetPz((TMath::Power(Q,4)*(-1 + xB) - 2*TMath::Power(M,2)*t*TMath::Power(xB,2) + 
     2*TMath::Power(Q,2)*xB*(t*(-1 + xB) - TMath::Power(M,2)*xB))/
   (2.*xB*(TMath::Power(Q,2)*(-1 + xB) - TMath::Power(M,2)*xB)*
     TMath::Sqrt((TMath::Power(Q,4) + 4*TMath::Power(M,2)*TMath::Power(Q,2)*TMath::Power(xB,2))/
       (xB*(TMath::Power(Q,2)*(1 - xB) + TMath::Power(M,2)*xB)))));
		pCM.SetPz(-qCM.Pz());
		ppCM.SetPz(-qpCM.Pz());
	
		qpPerp=qpCM.Px();
	
		thetag=TMath::ACos(qpCM.Pz()/qpCM.E());
		
		
		// Omega
		
		y = (-TMath::Power(Q,2) + 4*ELab*M*xB)/
   TMath::Sqrt(TMath::Power(Q,4) + 4*TMath::Power(M,2)*TMath::Power(Q,2)*TMath::Power(xB,2));
		if ( y < 1. )
		{
			cout << "TGVKelly : ELab " << ELab << " GeV is too small. TMath::CosH[ome]= " << y << endl;
			exit(-1);
		} // end if y
		Omega = TMath::Log(y + TMath::Sqrt(-1 + y)*TMath::Sqrt(1 + y));
		
		
		// PhaseSpace
		
		PhaseSpace = (TMath::Power(e,6)*TMath::Power(Q,2))/
   (4096.*TMath::Power(TMath::Pi(),5)*TMath::Sqrt(4*TMath::Power(M,2)*TMath::Power(Q,2) + TMath::Power(Q,4)/TMath::Power(xB,2))*
     TMath::Power(xB,2)*TMath::Power(TMath::Power(Q,2)/(4.*xB) + 
       (TMath::Sqrt(4*TMath::Power(M,2)*TMath::Power(Q,2) + TMath::Power(Q,4)/TMath::Power(xB,2))*TMath::CosH(Omega))/4.,
      2));
		
		
		// Flag
		
		InitKinematics = kTRUE;
    
    
    	// Validation mode only
    
    	if ( Validation == kTRUE )
    	{	
    		ofstream outfile;    		
    		
    		outfile.open("CheckKinematics.dat");
			outfile << Q << endl;
			outfile <<  qpPerp << endl;
			outfile <<  s << endl;
			outfile << t << endl;
			outfile <<  thetag << endl;
			outfile <<  xB << endl;
			outfile <<  pCM.E() << endl;
			outfile <<  ppCM.E() << endl;
			outfile <<  qCM.E() << endl;
			outfile <<  qCM.Pz() << endl;
			outfile <<  qpCM.E() << endl;
			outfile <<  qpCM.Px() << endl;
			outfile <<  qpCM.Pz() << endl;
			outfile.close();
		} // end if Validation
		
	} // end if InitKinematics	
			
} // end MakeKinematics



/*----- Function MakeVCSHelicityAmplitudes(ReH,ImH,ReE,ImE,ReHT,ImHT,ReET,ImET) --------*
 | Computes the hadronic helicity amplitudes RMvcs and IMvcs at leading twist.          |
 | The formula depend on the real and imaginary part of the Compton Form Factor :       |
 | ReH, ImH, ReE, ImE, ReHT, ImHT, ReET, ImET.                                          |
 *--------------------------------------------------------------------------------------*/

void TGVKelly::MakeVCSHelicityAmplitudes( Double_t ReH, Double_t ImH, Double_t ReE, Double_t ImE, Double_t ReHT, Double_t ImHT, Double_t ReET, Double_t ImET )
{
//	if ( InitMvcs == kFALSE )
//	{
	 
        RMvcs[0][0]=0;
        RMvcs[0][1]=0;
        RMvcs[0][2]=-((-4*ReH*(-1 + xB) - ReE*TMath::Power(xB,2))/(TMath::Sqrt(1 - xB)*(-2 + xB)));
        RMvcs[1][0]=0;
        RMvcs[1][1]=0;
        RMvcs[1][2]=-((qpPerp*ReET*xB)/(TMath::Sqrt(1 - xB)*(-2*M + M*xB)));
        RMvcs[2][0]=0;
        RMvcs[2][1]=0;
        RMvcs[2][2]=(qpPerp*ReE)/(M*TMath::Sqrt(1 - xB));
        RMvcs[3][0]=0;
        RMvcs[3][1]=0;
        RMvcs[3][2]=-((-4*ReHT*(-1 + xB) - ReET*TMath::Power(xB,2))/(TMath::Sqrt(1 - xB)*(-2 + xB)));

        IMvcs[0][0]=0;
        IMvcs[0][1]=0;
        IMvcs[0][2]=-((-4*ImH*(-1 + xB) - ImE*TMath::Power(xB,2))/(TMath::Sqrt(1 - xB)*(-2 + xB)));
        IMvcs[1][0]=0;
        IMvcs[1][1]=0;
        IMvcs[1][2]=-((ImET*qpPerp*xB)/(TMath::Sqrt(1 - xB)*(-2*M + M*xB)));
        IMvcs[2][0]=0;
        IMvcs[2][1]=0;
        IMvcs[2][2]=(ImE*qpPerp)/(M*TMath::Sqrt(1 - xB));
        IMvcs[3][0]=0;
        IMvcs[3][1]=0;
        IMvcs[3][2]=-((-4*ImHT*(-1 + xB) - ImET*TMath::Power(xB,2))/(TMath::Sqrt(1 - xB)*(-2 + xB)));
    	
    	
    	// Flag
    	
    	InitMvcs = kTRUE;
    
    
    	// Validation mode only
    
    	if ( Validation == kTRUE )
    	{	
    		ofstream outfile;    		
    		
    		outfile.open("CheckMvcs.dat");
    		for(Int_t i=0;i<12;i++)
			{
				outfile << IMvcs[div(i,3).quot][div(i,3).rem] << endl;
			} // end for i
			for(Int_t i=12;i<24;i++)
			{
				outfile << RMvcs[div(i-12,3).quot][div(i-12,3).rem] << endl;
			} // end for i
			outfile.close();
		} // end if Validation
		
//	} // end if InitMvcs
	
} // end MakeHelicityAmplitudes


Double_t TGVKelly::KellyE(Double_t q2)
  /* JJ Kelly PRC 70, 068202 (2004)
   */
{
  Int_t ia=1, ib=3;
  Double_t a[1]={-0.24};
  Double_t b[3]={10.98, 12.82, 21.97};
  Double_t Mp = 0.938271998;
  Double_t tau = -q2/(4.*pow(Mp,2));
  Double_t GEKn = 1.0;
  Double_t GEKd = 1.0;
  for (int ja=0; ja<ia; ja++){
    GEKn+=a[ja]*pow(tau,ja+1);
  }
  for (int jb=0; jb<ib; jb++){
    GEKd+=b[jb]*pow(tau,jb+1);
  }
  return GEKn/GEKd;
}



Double_t TGVKelly::KellyM(Double_t q2)
  /* JJ Kelly PRC 70, 068202 (2004)
     Magnetic Form Factor fit
  */
{
  Int_t ia=1, ib=3;
  Double_t a[1]={0.12};
  Double_t b[3]={10.97, 18.86, 6.55};
  Double_t Mp = 0.938271998;
  Double_t tau = -q2/(4.*pow(Mp,2));
  Double_t GMKn = 1.0;
  Double_t GMKd = 1.0;
  for (int ja=0; ja<ia; ja++)
    {
      GMKn+=a[ja]*pow(tau,ja+1);
    }
  for (int jb=0; jb<ib; jb++)
    {
      GMKd+=b[jb]*pow(tau,jb+1);
    }
  return 2.79285*GMKn/GMKd;
}


	
/*----------------------- Function MakeExactBHCrossSections() --------------------------*
 | Computes all the stuff to evaluate the cross section assuming the hadronic helicity  |
 | amplitudes are given, i.e. initializes :                                             |
 |   - the helicity amplitudes Jem,                                                     |
 |   - the SigmaBHPol's.                                                                |
 | No approximations apart those in the computation of the hadronic helicity amplitudes |
 | in terms of the Compton Form Factors.                                                |
 *--------------------------------------------------------------------------------------*/
 
void TGVKelly::MakeExactBHCrossSections()
{
	if ( InitExactBHCrossSections == kFALSE )
	{
		Double_t F1; // Dirac form factor
		Double_t F2; // Pauli form factor
		Double_t Ge, Gm; // Sachs' parametrization	
	
/*		F1=(4.*TMath::Power(M,2) - 2.79285*t)/
   (TMath::Power(1. - 1.4084507042253522*t,2)*(4.*TMath::Power(M,2) - 1.*t)); 
		F2=(7.1714*TMath::Power(M,2))/(TMath::Power(1 - 1.4084507042253522*t,2)*(4*TMath::Power(M,2) - t)); 
		Gm=2.79285/TMath::Power(1 - 1.4084507042253522*t,2);
		Ge=TMath::Power(1 - 1.4084507042253522*t,-2); */
	
		Gm=KellyM( t );
		Ge=KellyE( t ); 
		F1=((-1.)*t/4./TMath::Power(0.938271998,2.)*Gm+Ge)/(1.-t/4./TMath::Power(0.938271998,2.)); 
  		F2=( Gm-Ge)/(1.-t/4./TMath::Power(0.938271998,2.));	
	

/*----------------- Helicity amplitudes of the interference process --------------------*/

        Jem[0][0]=(-2*TMath::Sqrt(2.)*qCM.Pz()*qpCM.Px()*F2*
      (TMath::Sqrt(pCM.E() - M)*TMath::Sqrt(ppCM.E() - M) - 
        TMath::Sqrt(pCM.E() + M)*TMath::Sqrt(ppCM.E() + M))*TMath::Sqrt(s)*TMath::Cos(thetag/2.))/
    (M*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t)) - 
   (2*TMath::Sqrt(2.)*(qCM.Pz() + qpCM.E())*Gm*
      (TMath::Sqrt(ppCM.E() - M)*TMath::Sqrt(pCM.E() + M) + 
        TMath::Sqrt(pCM.E() - M)*TMath::Sqrt(ppCM.E() + M))*TMath::Sqrt(s)*TMath::Sin(thetag/2.))/
    TMath::Sqrt(TMath::Power(Q,4) - 4*s*t);
        Jem[0][1]=-((F2*(TMath::Sqrt(pCM.E() - M)*TMath::Sqrt(ppCM.E() - M) - 
          TMath::Sqrt(pCM.E() + M)*TMath::Sqrt(ppCM.E() + M))*
        (2*TMath::Power(M,2) + TMath::Power(Q,2) + 2*s)*TMath::Sqrt(-(s*t))*TMath::Cos(thetag/2.))/
      (M*TMath::Sqrt(s)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t))) - 
   (8*Gm*(TMath::Sqrt(ppCM.E() - M)*TMath::Sqrt(pCM.E() + M) + 
        TMath::Sqrt(pCM.E() - M)*TMath::Sqrt(ppCM.E() + M))*s*TMath::Sqrt(-t)*
      ((qCM.Pz() - qpCM.Pz())*TMath::Cos(thetag/2.) - 
        qpCM.Px()*TMath::Sin(thetag/2.)))/(TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Jem[0][2]=(2*TMath::Sqrt(2.)*qCM.Pz()*qpCM.Px()*F2*
      (TMath::Sqrt(pCM.E() - M)*TMath::Sqrt(ppCM.E() - M) - 
        TMath::Sqrt(pCM.E() + M)*TMath::Sqrt(ppCM.E() + M))*TMath::Sqrt(s)*TMath::Cos(thetag/2.))/
    (M*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t)) + 
   (2*TMath::Sqrt(2.)*(qCM.Pz() + qpCM.E())*Gm*
      (TMath::Sqrt(ppCM.E() - M)*TMath::Sqrt(pCM.E() + M) + 
        TMath::Sqrt(pCM.E() - M)*TMath::Sqrt(ppCM.E() + M))*TMath::Sqrt(s)*TMath::Sin(thetag/2.))/
    TMath::Sqrt(TMath::Power(Q,4) - 4*s*t);
        Jem[1][0]=TMath::Sqrt(2.)*Gm*(-(TMath::Sqrt(ppCM.E() - M)*TMath::Sqrt(pCM.E() + M)) + 
     TMath::Sqrt(pCM.E() - M)*TMath::Sqrt(ppCM.E() + M))*TMath::Cos(thetag/2.);
        Jem[1][1]=0;
        Jem[1][2]=TMath::Sqrt(2.)*Gm*(-(TMath::Sqrt(ppCM.E() - M)*TMath::Sqrt(pCM.E() + M)) + 
     TMath::Sqrt(pCM.E() - M)*TMath::Sqrt(ppCM.E() + M))*TMath::Cos(thetag/2.);
        Jem[2][0]=(2*TMath::Sqrt(2.)*(qCM.Pz() - qpCM.E())*Gm*
      (TMath::Sqrt(ppCM.E() - M)*TMath::Sqrt(pCM.E() + M) - 
        TMath::Sqrt(pCM.E() - M)*TMath::Sqrt(ppCM.E() + M))*TMath::Sqrt(s)*TMath::Cos(thetag/2.))/
    TMath::Sqrt(TMath::Power(Q,4) - 4*s*t) - (2*TMath::Sqrt(2.)*qCM.Pz()*qpCM.Px()*F2*
      (TMath::Sqrt(pCM.E() - M)*TMath::Sqrt(ppCM.E() - M) + 
        TMath::Sqrt(pCM.E() + M)*TMath::Sqrt(ppCM.E() + M))*TMath::Sqrt(s)*TMath::Sin(thetag/2.))/
    (M*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Jem[2][1]=-((F2*(TMath::Sqrt(pCM.E() - M)*TMath::Sqrt(ppCM.E() - M) + 
          TMath::Sqrt(pCM.E() + M)*TMath::Sqrt(ppCM.E() + M))*
        (2*TMath::Power(M,2) + TMath::Power(Q,2) + 2*s)*TMath::Sqrt(-(s*t))*TMath::Sin(thetag/2.))/
      (M*TMath::Sqrt(s)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t))) - 
   (8*Gm*(TMath::Sqrt(ppCM.E() - M)*TMath::Sqrt(pCM.E() + M) - 
        TMath::Sqrt(pCM.E() - M)*TMath::Sqrt(ppCM.E() + M))*s*TMath::Sqrt(-t)*
      (qpCM.Px()*TMath::Cos(thetag/2.) + 
        (qCM.Pz() - qpCM.Pz())*TMath::Sin(thetag/2.)))/
    (TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Jem[2][2]=(-2*TMath::Sqrt(2.)*(qCM.Pz() - qpCM.E())*Gm*
      (TMath::Sqrt(ppCM.E() - M)*TMath::Sqrt(pCM.E() + M) - 
        TMath::Sqrt(pCM.E() - M)*TMath::Sqrt(ppCM.E() + M))*TMath::Sqrt(s)*TMath::Cos(thetag/2.))/
    TMath::Sqrt(TMath::Power(Q,4) - 4*s*t) + (2*TMath::Sqrt(2.)*qCM.Pz()*qpCM.Px()*F2*
      (TMath::Sqrt(pCM.E() - M)*TMath::Sqrt(ppCM.E() - M) + 
        TMath::Sqrt(pCM.E() + M)*TMath::Sqrt(ppCM.E() + M))*TMath::Sqrt(s)*TMath::Sin(thetag/2.))/
    (M*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Jem[3][0]=TMath::Sqrt(2.)*Gm*(TMath::Sqrt(ppCM.E() - M)*TMath::Sqrt(pCM.E() + M) + 
     TMath::Sqrt(pCM.E() - M)*TMath::Sqrt(ppCM.E() + M))*TMath::Sin(thetag/2.);
        Jem[3][1]=0;
        Jem[3][2]=TMath::Sqrt(2.)*Gm*(TMath::Sqrt(ppCM.E() - M)*TMath::Sqrt(pCM.E() + M) + 
     TMath::Sqrt(pCM.E() - M)*TMath::Sqrt(ppCM.E() + M))*TMath::Sin(thetag/2.);


/*--------------------------------BH cross sections ------------------------------------*/

        SigmaBHPol0[0]=(-32*TMath::Power(Ge,2)*TMath::Power(M,2)*(3*TMath::Power(M,8)*t + 
        3*TMath::Power(M,6)*(TMath::Power(Q,4) + TMath::Power(Q,2)*t - 4*s*t) + 
        TMath::Power(M,2)*(3*TMath::Power(Q,4)*TMath::Power(TMath::Power(Q,2) + s,2) - 
           (TMath::Power(Q,2) + s)*(TMath::Power(Q,4) + 3*TMath::Power(Q,2)*s + 12*TMath::Power(s,2))*t + 
           2*(TMath::Power(Q,4) + TMath::Power(Q,2)*s - 3*TMath::Power(s,2))*TMath::Power(t,2)) + 
        TMath::Power(TMath::Power(Q,2) + s,2)*t*(3*s*(s + t) + TMath::Power(Q,2)*(3*s + t)) + 
        TMath::Power(M,4)*(4*TMath::Power(Q,6) + TMath::Power(Q,2)*t*(3*s + t) - 
           TMath::Power(Q,4)*(6*s + t) + 3*s*t*(6*s + t))) + 
     4*TMath::Power(Gm,2)*t*(6*TMath::Power(M,8)*t + 
        t*(3*TMath::Power(TMath::Power(Q,2) + s,2)*
            (TMath::Power(Q,4) + 2*TMath::Power(Q,2)*s + 2*TMath::Power(s,2)) + 
           2*s*(TMath::Power(Q,2) + s)*(2*TMath::Power(Q,2) + 3*s)*t + 
           (3*TMath::Power(Q,4) + 4*TMath::Power(Q,2)*s + 3*TMath::Power(s,2))*TMath::Power(t,2)) - 
        2*TMath::Power(M,6)*(3*TMath::Power(Q,4) - 11*TMath::Power(Q,2)*t + 6*t*(2*s + t)) + 
        TMath::Power(M,4)*(-8*TMath::Power(Q,6) - 26*TMath::Power(Q,2)*t*(s + t) + 
           TMath::Power(Q,4)*(12*s + 25*t) + 
           3*t*(12*TMath::Power(s,2) + 10*s*t + TMath::Power(t,2))) - 
        2*TMath::Power(M,2)*(3*TMath::Power(Q,8) + TMath::Power(Q,6)*(6*s - 5*t) + 
           3*s*t*TMath::Power(2*s + t,2) + 
           TMath::Power(Q,2)*t*(7*TMath::Power(s,2) + 2*s*t - 3*TMath::Power(t,2)) + 
           TMath::Power(Q,4)*(3*TMath::Power(s,2) - 5*s*t + 7*TMath::Power(t,2)))))/
   ((TMath::Power(M,4) + 2*TMath::Power(M,2)*(TMath::Power(Q,2) - s) + TMath::Power(TMath::Power(Q,2) + s,2))*
     (4*TMath::Power(M,2) - t));
        SigmaBHPol0[1]=(-32*TMath::Power(Ge,2)*TMath::Power(M,2)*(TMath::Power(M,8)*t + 
        TMath::Power(M,6)*(TMath::Power(Q,4) + TMath::Power(Q,2)*t - 4*s*t) + 
        TMath::Power(M,2)*(TMath::Power(Q,4)*TMath::Power(TMath::Power(Q,2) + s,2) + 
           (TMath::Power(Q,2) - s)*(TMath::Power(Q,2) + s)*(5*TMath::Power(Q,2) + 4*s)*t - 
           2*(TMath::Power(Q,4) + TMath::Power(Q,2)*s + TMath::Power(s,2))*TMath::Power(t,2)) + 
        TMath::Power(TMath::Power(Q,2) + s,2)*t*(TMath::Power(Q,2)*(s - t) + s*(s + t)) + 
        TMath::Power(M,4)*(-4*TMath::Power(Q,6) + TMath::Power(Q,2)*(s - t)*t + s*t*(6*s + t) + 
           TMath::Power(Q,4)*(-2*s + 5*t))) + 
     4*TMath::Power(Gm,2)*t*(2*TMath::Power(M,8)*t + 
        t*(TMath::Power(TMath::Power(Q,2) + s,2)*
            (TMath::Power(Q,4) + 2*TMath::Power(Q,2)*s + 2*TMath::Power(s,2)) + 
           2*s*(-2*TMath::Power(Q,4) - TMath::Power(Q,2)*s + TMath::Power(s,2))*t + 
           (TMath::Power(Q,4) - 4*TMath::Power(Q,2)*s + TMath::Power(s,2))*TMath::Power(t,2)) - 
        2*TMath::Power(M,6)*(TMath::Power(Q,4) - 9*TMath::Power(Q,2)*t + 2*t*(2*s + t)) + 
        TMath::Power(M,4)*(8*TMath::Power(Q,6) - 2*TMath::Power(Q,2)*t*(15*s + 7*t) + 
           TMath::Power(Q,4)*(4*s + 19*t) + t*(12*TMath::Power(s,2) + 10*s*t + TMath::Power(t,2)))\
         - 2*TMath::Power(M,2)*(TMath::Power(Q,8) + TMath::Power(Q,6)*(2*s + t) + 
           s*t*TMath::Power(2*s + t,2) - 
           TMath::Power(Q,2)*t*(3*TMath::Power(s,2) + 10*s*t + TMath::Power(t,2)) + 
           TMath::Power(Q,4)*(TMath::Power(s,2) - 7*s*t + 5*TMath::Power(t,2)))))/
   ((TMath::Power(M,4) + 2*TMath::Power(M,2)*(TMath::Power(Q,2) - s) + TMath::Power(TMath::Power(Q,2) + s,2))*
     (4*TMath::Power(M,2) - t));
        SigmaBHPol0[2]=(-64*TMath::Power(Ge,2)*TMath::Power(M,2)*Q*qpPerp*(TMath::Power(M,2) - TMath::Power(Q,2) - s)*
      (2*TMath::Power(M,2)*TMath::Power(Q,2) - (TMath::Power(M,2) + TMath::Power(Q,2) + s)*t) + 
     16*TMath::Power(Gm,2)*Q*qpPerp*t*(TMath::Power(M,4)*(-2*TMath::Power(Q,2) + 3*t) + 
        t*(TMath::Power(Q,2)*(s - t) + s*(s + t)) + 
        TMath::Power(M,2)*(2*TMath::Power(Q,4) - t*(4*s + t) + TMath::Power(Q,2)*(2*s + 5*t))))/
   (TMath::Sqrt(TMath::Power(M,4) + 2*TMath::Power(M,2)*(TMath::Power(Q,2) - s) + 
       TMath::Power(TMath::Power(Q,2) + s,2))*(4*TMath::Power(M,2) - t));
        SigmaBHPol0[3]=(-64*TMath::Power(Ge,2)*TMath::Power(M,4)*TMath::Power(Q,2)*
      (TMath::Power(M,4)*t + s*t*(TMath::Power(Q,2) + s + t) + 
        TMath::Power(M,2)*(TMath::Power(Q,4) + TMath::Power(Q,2)*t - 2*s*t)) - 
     8*TMath::Power(Gm,2)*TMath::Power(Q,2)*(2*TMath::Power(M,2) - t)*t*
      (TMath::Power(M,4)*t + s*t*(TMath::Power(Q,2) + s + t) + 
        TMath::Power(M,2)*(TMath::Power(Q,4) + TMath::Power(Q,2)*t - 2*s*t)))/
   ((TMath::Power(M,4) + 2*TMath::Power(M,2)*(TMath::Power(Q,2) - s) + TMath::Power(TMath::Power(Q,2) + s,2))*
     (4*TMath::Power(M,2) - t));
        SigmaBHPolX[0]=(64*Ge*Gm*M*qpPerp*(-(TMath::Power(Q,2)*(TMath::Power(Q,2) + s)) + 
        TMath::Power(M,2)*(TMath::Power(Q,2) - t) + s*t)*
      (2*TMath::Power(M,2)*TMath::Power(Q,2) - (TMath::Power(M,2) + TMath::Power(Q,2) + s)*t) + 
     32*TMath::Power(Gm,2)*M*qpPerp*t*(-TMath::Power(Q,6) - 3*TMath::Power(Q,4)*s - 
        2*TMath::Power(M,4)*(TMath::Power(Q,2) - t) + s*t*(2*s + t) - 
        TMath::Power(Q,2)*(2*TMath::Power(s,2) + TMath::Power(t,2)) - 
        TMath::Power(M,2)*(TMath::Power(Q,4) - 4*TMath::Power(Q,2)*(s + t) + t*(4*s + t))))/
   (TMath::Sqrt(TMath::Power(M,4) + 2*TMath::Power(M,2)*(TMath::Power(Q,2) - s) + 
       TMath::Power(TMath::Power(Q,2) + s,2))*(4*TMath::Power(M,2) - t));
        SigmaBHPolX[1]=(32*Ge*Gm*M*Q*(-2*TMath::Power(M,2)*TMath::Power(Q,2) + (TMath::Power(M,2) + TMath::Power(Q,2) + s)*t)*
      (4*TMath::Power(M,2)*TMath::Power(Q,4) + 
        (2*TMath::Power(M,2) + TMath::Power(Q,2) - 2*s)*(TMath::Power(M,2) - TMath::Power(Q,2) - s)*t + 
        (TMath::Power(M,2) + TMath::Power(Q,2) + 3*s)*TMath::Power(t,2)) + 
     64*TMath::Power(Gm,2)*M*Q*(TMath::Power(M,2) - s - t)*t*
      (TMath::Power(M,4)*t + s*t*(TMath::Power(Q,2) + s + t) + 
        TMath::Power(M,2)*(TMath::Power(Q,4) + TMath::Power(Q,2)*t - 2*s*t)))/
   ((TMath::Power(M,4) + 2*TMath::Power(M,2)*(TMath::Power(Q,2) - s) + TMath::Power(TMath::Power(Q,2) + s,2))*
     (4*TMath::Power(M,2) - t));
        SigmaBHPolY=32*Ge*Gm*M*Q*(TMath::Power(Q,2) - t)*t;
        SigmaBHPolZ[0]=(-128*Ge*Gm*TMath::Power(M,2)*(-(TMath::Power(Q,2)*(TMath::Power(Q,2) + s)) + 
        TMath::Power(M,2)*(TMath::Power(Q,2) - t) + s*t)*
      (TMath::Power(M,4)*t + s*t*(TMath::Power(Q,2) + s + t) + 
        TMath::Power(M,2)*(TMath::Power(Q,4) + TMath::Power(Q,2)*t - 2*s*t)) + 
     16*TMath::Power(Gm,2)*t*(2*TMath::Power(M,2)*TMath::Power(Q,2) - 
        (TMath::Power(M,2) + TMath::Power(Q,2) + s)*t)*
      (TMath::Power(Q,6) + 3*TMath::Power(Q,4)*s + 2*TMath::Power(M,4)*(TMath::Power(Q,2) - t) - 
        s*t*(2*s + t) + TMath::Power(Q,2)*(2*TMath::Power(s,2) + TMath::Power(t,2)) + 
        TMath::Power(M,2)*(TMath::Power(Q,4) - 4*TMath::Power(Q,2)*(s + t) + t*(4*s + t))))/
   ((TMath::Power(M,4) + 2*TMath::Power(M,2)*(TMath::Power(Q,2) - s) + TMath::Power(TMath::Power(Q,2) + s,2))*
     (4*TMath::Power(M,2) - t));
        SigmaBHPolZ[1]=(64*Ge*Gm*TMath::Power(M,2)*Q*qpPerp*(-4*TMath::Power(M,2)*TMath::Power(Q,4) - 
        (2*TMath::Power(M,2) + TMath::Power(Q,2) - 2*s)*(TMath::Power(M,2) - TMath::Power(Q,2) - s)*t - 
        (TMath::Power(M,2) + TMath::Power(Q,2) + 3*s)*TMath::Power(t,2)) + 
     32*TMath::Power(Gm,2)*Q*qpPerp*t*(-TMath::Power(M,2) + s + t)*
      ((TMath::Power(Q,2) + s)*t + TMath::Power(M,2)*(-2*TMath::Power(Q,2) + t)))/
   (TMath::Sqrt(TMath::Power(M,4) + 2*TMath::Power(M,2)*(TMath::Power(Q,2) - s) + 
       TMath::Power(TMath::Power(Q,2) + s,2))*(4*TMath::Power(M,2) - t));
    
    
   	 	// Flag
   	 	
   	 	InitExactBHCrossSections = kTRUE;
   	 	
   	 	
   	 	// Validation mode only
    
   	 	if ( Validation == kTRUE )
   	 	{	
    		ofstream outfile;    		
    		
    		outfile.open("CheckFormFactors.dat");
			outfile << F1 << endl;
			outfile <<  F2 << endl;
			outfile <<  Ge << endl;
			outfile << Gm << endl;
			outfile.close();
		
		
    		outfile.open("CheckJemExact.dat");
    		for(Int_t i=0;i<12;i++)
			{
				outfile << Jem[div(i,4).rem][div(i,4).quot] << endl;	
			} // end for i		
			outfile.close();
		
		
    		outfile.open("CheckSigmaBHExact.dat");
    		for(Int_t i=0;i<9;i++)
			{
				if (0 <= i && i < 4)
				{
					outfile << SigmaBHPol0[i] << endl;
				} // end if i
				if (4 <= i && i < 6)
				{
					outfile << SigmaBHPolX[i-4] << endl;
				} // end if i
				if (i == 6)
				{
					outfile << SigmaBHPolY << endl;
				} // end if i
				if (7<= i && i < 9)
				{
					outfile << SigmaBHPolZ[i-7] << endl;
				} // end if i
			} // end for i
			outfile.close();		
				
		} // end if Validation  
		
	} // end if InitExactBHCrossSections
	 
} // end MakeExactBHCrossSections


	
/*------------------ Function MakeExactVCSAndInterfCrossSections() ---------------------*
 | Computes all the stuff to evaluate the cross section assuming the hadronic helicity  |
 | amplitudes are given, i.e. initializes :                                             |
 |   - the helicity amplitudes Mvcs,                                                    |
 |   - the expansion coefficients Ur,                                                   |
 |   - the SigmaVCSPol's and SigmaIPol's.                                               |
 | No approximations apart those in the computation of the hadronic helicity amplitudes |
 | in terms of the Compton Form Factors.                                                |
 *--------------------------------------------------------------------------------------*/
 
void TGVKelly::MakeExactVCSAndInterfCrossSections()
{
//	if ( InitExactVCSAndInterfCrossSections == kFALSE )
//	{

/*-------------- Harmonic expansion coefficients of the VCS cross section --------------*/

        Ur[0]=-((TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t)*
      (TMath::Power(Q,4) + 4*qCM.E()*(4*qCM.E() - 3*qCM.Pz())*
         TMath::Power(qpPerp,2) + (qCM.E() + qCM.Pz())*
         (qCM.E() + qpCM.E())*t + 
        TMath::Power(Q,2)*((qCM.E() - qpCM.E())*
            (qCM.E() + qCM.Pz() + 2*qpCM.E()) + 
           16*TMath::Power(qpPerp,2) + t)))/
   (64.*qCM.Pz()*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[1]=((TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t)*
     (TMath::Power(Q,4) - 4*qCM.E()*qCM.Pz()*TMath::Power(qpPerp,2) + 
       (qCM.E() - qCM.Pz())*(qCM.E() + qpCM.E())*t + 
       TMath::Power(Q,2)*((qCM.E() - qpCM.E())*
           (qCM.E() - qCM.Pz() + 2*qpCM.E()) + t)))/
   (64.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[2]=((TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t)*
     (TMath::Power(Q,4) + 4*qCM.E()*qCM.Pz()*TMath::Power(qpPerp,2) + 
       (qCM.E() + qCM.Pz())*(qCM.E() + qpCM.E())*t + 
       TMath::Power(Q,2)*((qCM.E() - qpCM.E())*
           (qCM.E() + qCM.Pz() + 2*qpCM.E()) + t)))/
   (64.*qCM.Pz()*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[3]=(qpPerp*(TMath::Power(Q,4) + TMath::Power(Q,2)*
        (12*TMath::Power(qCM.E(),2) - 16*qCM.E()*qCM.Pz() - 
          22*qCM.E()*qpCM.E() + 26*qCM.Pz()*qpCM.E() - 7*t) + 
       4*qCM.E()*(qCM.E() - 4*qCM.Pz())*t)*
     (TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t))
    /(128.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[4]=(3*qpPerp*(2*(qCM.E() + qCM.Pz())*qpCM.E() + TMath::Power(Q,2) + t)*
     (TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t))
    /(128.*qCM.Pz()*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[5]=(qpPerp*(5*TMath::Power(Q,4) + TMath::Power(Q,2)*
        (2*(6*TMath::Power(qCM.E(),2) - 7*qCM.E()*qpCM.E() + 
             qCM.Pz()*qpCM.E()) - 3*t) + 4*TMath::Power(qCM.E(),2)*t)*
     (TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t))
    /(128.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[6]=-((TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t)*
      (TMath::Power(Q,4) + TMath::Power(Q,2)*((qCM.E() - qpCM.E())*
            (qCM.E() - qCM.Pz() + 2*qpCM.E()) + 
           2*TMath::Power(qpPerp,2) + t) + 
        (qCM.E() - qCM.Pz())*
         (2*qCM.E()*TMath::Power(qpPerp,2) + (qCM.E() + qpCM.E())*t)))/
   (32.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[7]=-(qpPerp*(qCM.E()*(3*qCM.E() - qCM.Pz() - 2*qpCM.E())*
         TMath::Power(Q,2) + 2*TMath::Power(Q,4) + 
        qCM.E()*(qCM.E() - qCM.Pz())*t)*
      (TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t)
      )/(16.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*
     TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[8]=(qpPerp*TMath::Sqrt(s)*((4*TMath::Power(qCM.E(),2) - 
          3*qCM.E()*(qCM.Pz() - qpCM.E()) - 
          4*qCM.Pz()*qpCM.E())*t + 
       TMath::Power(Q,2)*(-3*(qCM.E() + qCM.Pz() - 2*qpCM.E())*
           (qCM.E() - qpCM.E()) + 7*t)))/
   (4.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[9]=-(qpPerp*TMath::Sqrt(s)*(qCM.E()*(qCM.Pz() + qpCM.E())*t + 
        TMath::Power(Q,2)*(-((qCM.E() - qCM.Pz() - 2*qpCM.E())*
              (qCM.E() - qpCM.E())) + t)))/
   (4.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[10]=(qpPerp*TMath::Sqrt(s)*(qCM.E()*(-qCM.Pz() + qpCM.E())*t + 
       TMath::Power(Q,2)*(-((qCM.E() + qCM.Pz() - 2*qpCM.E())*
             (qCM.E() - qpCM.E())) + t)))/
   (4.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[11]=(TMath::Sqrt(s)*(3*TMath::Power(Q,6) + TMath::Power(Q,4)*
        (3*TMath::Power(qCM.E(),2) - 7*qCM.E()*qCM.Pz() + 
          35*qCM.E()*qpCM.E() + 3*qCM.Pz()*qpCM.E() - 
          22*TMath::Power(qpCM.E(),2) + 12*TMath::Power(qpPerp,2) + 3*t) + 
       8*qCM.E()*(4*qCM.E() - qCM.Pz())*
        (qCM.E() + qpCM.E())*
        (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
       TMath::Power(Q,2)*(8*(4*qCM.E() - qCM.Pz())*qpCM.E()*
           (TMath::Power(qCM.E(),2) - TMath::Power(qpCM.E(),2)) + 
          2*(22*TMath::Power(qCM.E(),2) + 7*qCM.Pz()*qpCM.E() + 
             qCM.E()*(-15*qCM.Pz() + 16*qpCM.E()))*TMath::Power(qpPerp,2)
            + (3*TMath::Power(qCM.E(),2) - 
             7*qCM.E()*(qCM.Pz() - 5*qpCM.E()) + 
             qpCM.E()*(-11*qCM.Pz() + 16*qpCM.E()))*t)))/
   (16.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[12]=(-3*TMath::Sqrt(s)*(TMath::Power(Q,4) + 2*(2*TMath::Power(qCM.E(),2) - 
          qCM.E()*qCM.Pz() + qCM.Pz()*qpCM.E())*
        TMath::Power(qpPerp,2) + (qCM.E() - qCM.Pz())*
        (qCM.E() + qpCM.E())*t + 
       TMath::Power(Q,2)*((qCM.E() - qpCM.E())*
           (qCM.E() - qCM.Pz() + 2*qpCM.E()) + 
          4*TMath::Power(qpPerp,2) + t)))/
   (16.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[13]=-(TMath::Sqrt(s)*(TMath::Power(Q,6) + TMath::Power(Q,4)*
         (TMath::Power(qCM.E(),2) + 3*qCM.E()*qCM.Pz() + 
           qCM.E()*qpCM.E() + qCM.Pz()*qpCM.E() - 
           2*TMath::Power(qpCM.E(),2) + 4*TMath::Power(qpPerp,2) + t) + 
        8*qCM.E()*qCM.Pz()*(qCM.E() + qpCM.E())*
         (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
        TMath::Power(Q,2)*(8*qCM.Pz()*qpCM.E()*
            (TMath::Power(qCM.E(),2) - TMath::Power(qpCM.E(),2)) + 
           (4*TMath::Power(qCM.E(),2) + 22*qCM.E()*qCM.Pz() - 
              6*qCM.Pz()*qpCM.E())*TMath::Power(qpPerp,2) + 
           (TMath::Power(qCM.E(),2) + 7*qCM.Pz()*qpCM.E() + 
              qCM.E()*(3*qCM.Pz() + qpCM.E()))*t)))/
   (16.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[14]=((qCM.E() - qpCM.E())*qpPerp*TMath::Sqrt(s)*
     ((-qCM.E() + qCM.Pz() + 2*qpCM.E())*TMath::Power(Q,2) + 
       (-qCM.E() + qCM.Pz())*t))/
   (4.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[15]=-(TMath::Sqrt(s)*(TMath::Power(Q,6) + TMath::Power(Q,4)*
         (TMath::Power(qCM.E(),2) + 3*qCM.E()*qpCM.E() - 
           qCM.Pz()*qpCM.E() - 3*TMath::Power(qpCM.E(),2) + 
           4*TMath::Power(qpPerp,2) + t) + 
        2*qCM.E()*(qCM.E() - qCM.Pz())*
         (qCM.E() + qpCM.E())*
         (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
        TMath::Power(Q,2)*(2*(qCM.E() - qCM.Pz())*qpCM.E()*
            (TMath::Power(qCM.E(),2) - TMath::Power(qpCM.E(),2)) + 
           2*qCM.E()*(3*qCM.E() - 2*qCM.Pz() + qpCM.E())*
            TMath::Power(qpPerp,2) + (TMath::Power(qCM.E(),2) + 
              3*qCM.E()*qpCM.E() - qCM.Pz()*qpCM.E() + 
              TMath::Power(qpCM.E(),2))*t)))/
   (2.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[16]=(qpPerp*((qCM.E() - qCM.Pz() - 2*qpCM.E())*TMath::Power(Q,2) + 
       qCM.E()*t))/(8.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2));
        Ur[17]=-((3*qCM.E() - 7*qCM.Pz() + 10*qpCM.E())*TMath::Power(Q,4) + 
      8*qCM.E()*(4*qCM.E() - qCM.Pz())*
       (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
      TMath::Power(Q,2)*(8*(4*qCM.E() - qCM.Pz())*
          (qCM.E() - qpCM.E())*qpCM.E() + 
         (32*qCM.E() - 6*qCM.Pz())*TMath::Power(qpPerp,2) + 
         (3*qCM.E() - 7*qCM.Pz() + 16*qpCM.E())*t))/
   (32.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2));
        Ur[18]=(3*((-qCM.E() + qCM.Pz() + 2*qpCM.E())*TMath::Power(Q,2) + 
       2*qCM.Pz()*TMath::Power(qpPerp,2) + (-qCM.E() + qCM.Pz())*t))/
   (32.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2));
        Ur[19]=((qCM.E() + 3*qCM.Pz() - 2*qpCM.E())*TMath::Power(Q,4) + 
     TMath::Power(Q,2)*(-2*qCM.Pz()*
         (4*qpCM.E()*(-qCM.E() + qpCM.E()) + TMath::Power(qpPerp,2)) + 
        (qCM.E() + 3*qCM.Pz())*t) + 
     8*qCM.E()*qCM.Pz()*(qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t)
     )/(32.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2));
        Ur[20]=(qpPerp*((qCM.E() + qCM.Pz() - 2*qpCM.E())*TMath::Power(Q,2) + 
       (qCM.E() - 3*qCM.Pz())*t))/(8.*TMath::Sqrt(2.)*qCM.Pz());
        Ur[21]=(qpPerp*((-qCM.E() + qCM.Pz() + 2*qpCM.E())*TMath::Power(Q,2) + 
       (-qCM.E() + qCM.Pz())*t))/(8.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2));
        Ur[22]=(qpPerp*(3*(qCM.E() + qCM.Pz() - 2*qpCM.E())*TMath::Power(Q,2) + 
       (3*qCM.E() - qCM.Pz())*t))/(8.*TMath::Sqrt(2.)*qCM.Pz());
        Ur[23]=((qCM.E() - qpCM.E())*TMath::Power(Q,4) + 
     2*qCM.E()*(qCM.E() - qCM.Pz())*
      (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
     TMath::Power(Q,2)*(2*(qCM.E() - qCM.Pz())*(qCM.E() - qpCM.E()\
            )*qpCM.E() + 2*(qCM.E() + qCM.Pz())*
         TMath::Power(qpPerp,2) + (qCM.E() + qpCM.E())*t))/
   (4.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2));
        Ur[24]=(qpPerp*(TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 
       4*s*t)*(3*TMath::Power(Q,4) + 2*TMath::Power(qCM.E(),2)*t - 
       TMath::Power(Q,2)*(6*qCM.E()*(-qCM.E() + qpCM.E()) + t)))/
   (32.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*
     TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[25]=-(qpPerp*(2*qCM.E()*qpCM.E() + TMath::Power(Q,2) + t)*
      (TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t)
      )/(32.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[26]=((3*(qCM.E() - qpCM.E())*TMath::Power(Q,2) + 4*qCM.E()*TMath::Power(qpPerp,2) + 
       3*(qCM.E() + qpCM.E())*t)*
     (TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t))
    /(16.*TMath::Sqrt(2.)*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[27]=-(((qCM.E() - qpCM.E())*TMath::Power(Q,2) + 4*qCM.E()*TMath::Power(qpPerp,2) + 
        (qCM.E() + qpCM.E())*t)*
      (TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t)
      )/(16.*TMath::Sqrt(2.)*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[28]=-(qCM.E()*TMath::Power(qpPerp,2)*(TMath::Power(Q,4) + 
        2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t))/
   (4.*TMath::Sqrt(2.)*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[29]=-(TMath::Sqrt(s)*((qCM.E() + qpCM.E())*TMath::Power(Q,4) + 
        4*qCM.E()*(qCM.E() + qpCM.E())*
         (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
        TMath::Power(Q,2)*(4*TMath::Power(qCM.E(),2)*qpCM.E() - 
           4*TMath::Power(qpCM.E(),3) + 
           2*(5*qCM.E() - qpCM.E())*TMath::Power(qpPerp,2) + qCM.E()*t + 
           3*qpCM.E()*t)))/(8.*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[30]=-(TMath::Sqrt(s)*((qCM.E() - qpCM.E())*TMath::Power(Q,2) + 
        2*(qCM.E() - qpCM.E())*TMath::Power(qpPerp,2) + 
        (qCM.E() + qpCM.E())*t))/
   (8.*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[31]=-(qpPerp*TMath::Sqrt(s)*(qCM.E()*(2*qCM.E() - qpCM.E())*t + 
        TMath::Power(Q,2)*(TMath::Power(qCM.E(),2) - 3*qCM.E()*qpCM.E() + 
           2*TMath::Power(qpCM.E(),2) + t)))/
   (2.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[32]=-(qpPerp*TMath::Sqrt(s)*(qCM.E()*qpCM.E()*t + 
        TMath::Power(Q,2)*(-TMath::Power(qCM.E(),2) + 3*qCM.E()*qpCM.E() - 
           2*TMath::Power(qpCM.E(),2) + t)))/
   (2.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[33]=-(qpPerp*TMath::Sqrt(s)*(qCM.E()*(qCM.E() + qpCM.E())*t + 
        TMath::Power(Q,2)*(-TMath::Power(qCM.E(),2) + 3*qCM.E()*qpCM.E() - 
           2*TMath::Power(qpCM.E(),2) + 2*t)))/
   (2.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[34]=((qCM.E() - 2*qpCM.E())*TMath::Power(Q,5) + qCM.E()*TMath::Power(Q,3)*t)/
   (16.*qCM.Pz()*TMath::Power(Q,5));
        Ur[35]=(qpPerp*(TMath::Power(Q,2) + 2*t))/(4.*TMath::Power(Q,2));
        Ur[36]=-qpPerp/4.;
        Ur[37]=(qCM.E()*TMath::Power(Q,4) + 4*TMath::Power(qCM.E(),2)*
      (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
     TMath::Power(Q,2)*(4*qCM.E()*(qCM.E() - qpCM.E())*qpCM.E() + 
        4*qCM.E()*TMath::Power(qpPerp,2) + (qCM.E() + 2*qpCM.E())*t))/
   (8.*qCM.Pz()*TMath::Power(Q,2));
        Ur[38]=(qpPerp*(-TMath::Power(Q,2) + t))/(4.*TMath::Power(Q,2));
        Ur[39]=-((TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t)*
      (TMath::Power(Q,4) + 4*qCM.E()*(4*qCM.E() + 3*qCM.Pz())*
         TMath::Power(qpPerp,2) + (qCM.E() - qCM.Pz())*
         (qCM.E() + qpCM.E())*t + 
        TMath::Power(Q,2)*((qCM.E() - qpCM.E())*
            (qCM.E() - qCM.Pz() + 2*qpCM.E()) + 
           16*TMath::Power(qpPerp,2) + t)))/
   (64.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[40]=-(qpPerp*(TMath::Power(Q,4) + TMath::Power(Q,2)*
         (12*TMath::Power(qCM.E(),2) + 16*qCM.E()*qCM.Pz() - 
           22*qCM.E()*qpCM.E() - 26*qCM.Pz()*qpCM.E() - 7*t)\
         + 4*qCM.E()*(qCM.E() + 4*qCM.Pz())*t)*
      (TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t)
      )/(128.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*
     TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[41]=(-3*TMath::Power(Q,2)*qpPerp*(2*(qCM.E() - qCM.Pz())*qpCM.E() + 
       TMath::Power(Q,2) + t)*(TMath::Power(Q,4) + 
       2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t))/
   (128.*qCM.Pz()*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[42]=-(qpPerp*(TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 
        4*s*t)*(5*TMath::Power(Q,4) + 4*TMath::Power(qCM.E(),2)*t - 
        TMath::Power(Q,2)*(2*(-6*TMath::Power(qCM.E(),2) + 7*qCM.E()*qpCM.E() + 
              qCM.Pz()*qpCM.E()) + 3*t)))/
   (128.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[43]=-((TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t)*
      (TMath::Power(Q,4) + TMath::Power(Q,2)*((qCM.E() - qpCM.E())*
            (qCM.E() + qCM.Pz() + 2*qpCM.E()) + 
           2*TMath::Power(qpPerp,2) + t) + 
        (qCM.E() + qCM.Pz())*
         (2*qCM.E()*TMath::Power(qpPerp,2) + (qCM.E() + qpCM.E())*t)))/
   (32.*qCM.Pz()*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[44]=(qpPerp*(qCM.E()*(3*qCM.E() + qCM.Pz() - 2*qpCM.E())*
        TMath::Power(Q,2) + 2*TMath::Power(Q,4) + qCM.E()*(qCM.E() + qCM.Pz())*t
       )*(TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 
       4*s*t))/(16.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*
     TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[45]=-(qpPerp*TMath::Sqrt(s)*((4*TMath::Power(qCM.E(),2) + 4*qCM.Pz()*qpCM.E() + 
           3*qCM.E()*(qCM.Pz() + qpCM.E()))*t + 
        TMath::Power(Q,2)*(-3*(qCM.E() - qCM.Pz() - 2*qpCM.E())*
            (qCM.E() - qpCM.E()) + 7*t)))/
   (4.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[46]=(TMath::Sqrt(s)*(3*TMath::Power(Q,6) + TMath::Power(Q,4)*
        (3*TMath::Power(qCM.E(),2) + 7*qCM.E()*qCM.Pz() + 
          35*qCM.E()*qpCM.E() - 3*qCM.Pz()*qpCM.E() - 
          22*TMath::Power(qpCM.E(),2) + 12*TMath::Power(qpPerp,2) + 3*t) + 
       8*qCM.E()*(4*qCM.E() + qCM.Pz())*
        (qCM.E() + qpCM.E())*
        (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
       TMath::Power(Q,2)*(8*(4*qCM.E() + qCM.Pz())*qpCM.E()*
           (TMath::Power(qCM.E(),2) - TMath::Power(qpCM.E(),2)) + 
          2*(22*TMath::Power(qCM.E(),2) - 7*qCM.Pz()*qpCM.E() + 
             qCM.E()*(15*qCM.Pz() + 16*qpCM.E()))*TMath::Power(qpPerp,2)\
           + (3*TMath::Power(qCM.E(),2) + 
             7*qCM.E()*(qCM.Pz() + 5*qpCM.E()) + 
             qpCM.E()*(11*qCM.Pz() + 16*qpCM.E()))*t)))/
   (16.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,4)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[47]=(-3*TMath::Power(Q,2)*TMath::Sqrt(s)*(TMath::Power(Q,4) + 
       2*(2*TMath::Power(qCM.E(),2) + qCM.E()*qCM.Pz() - 
          qCM.Pz()*qpCM.E())*TMath::Power(qpPerp,2) + 
       (qCM.E() + qCM.Pz())*(qCM.E() + qpCM.E())*t + 
       TMath::Power(Q,2)*((qCM.E() - qpCM.E())*
           (qCM.E() + qCM.Pz() + 2*qpCM.E()) + 
          4*TMath::Power(qpPerp,2) + t)))/
   (16.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[48]=-(TMath::Sqrt(s)*(TMath::Power(Q,6) + TMath::Power(Q,4)*
         (TMath::Power(qCM.E(),2) - 3*qCM.E()*qCM.Pz() + 
           qCM.E()*qpCM.E() - qCM.Pz()*qpCM.E() - 
           2*TMath::Power(qpCM.E(),2) + 4*TMath::Power(qpPerp,2) + t) - 
        8*qCM.E()*qCM.Pz()*(qCM.E() + qpCM.E())*
         (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
        TMath::Power(Q,2)*(8*qCM.Pz()*qpCM.E()*
            (-TMath::Power(qCM.E(),2) + TMath::Power(qpCM.E(),2)) + 
           (4*TMath::Power(qCM.E(),2) - 22*qCM.E()*qCM.Pz() + 
              6*qCM.Pz()*qpCM.E())*TMath::Power(qpPerp,2) + 
           (TMath::Power(qCM.E(),2) - 7*qCM.Pz()*qpCM.E() + 
              qCM.E()*(-3*qCM.Pz() + qpCM.E()))*t)))/
   (16.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,4)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[49]=((qCM.E() - qpCM.E())*qpPerp*TMath::Sqrt(s)*
     ((qCM.E() + qCM.Pz() - 2*qpCM.E())*TMath::Power(Q,2) + 
       (qCM.E() + qCM.Pz())*t))/
   (4.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[50]=-(TMath::Sqrt(s)*(TMath::Power(Q,6) + TMath::Power(Q,4)*
         (TMath::Power(qCM.E(),2) + 3*qCM.E()*qpCM.E() + 
           qCM.Pz()*qpCM.E() - 3*TMath::Power(qpCM.E(),2) + 
           4*TMath::Power(qpPerp,2) + t) + 
        2*qCM.E()*(qCM.E() + qCM.Pz())*
         (qCM.E() + qpCM.E())*
         (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
        TMath::Power(Q,2)*(2*(qCM.E() + qCM.Pz())*qpCM.E()*
            (TMath::Power(qCM.E(),2) - TMath::Power(qpCM.E(),2)) + 
           2*qCM.E()*(3*qCM.E() + 2*qCM.Pz() + qpCM.E())*
            TMath::Power(qpPerp,2) + (TMath::Power(qCM.E(),2) + 
              3*qCM.E()*qpCM.E() + 
              qpCM.E()*(qCM.Pz() + qpCM.E()))*t)))/
   (2.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,4)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[51]=(qpPerp*((qCM.E() + qCM.Pz() - 2*qpCM.E())*TMath::Power(Q,2) + 
       qCM.E()*t))/(8.*TMath::Sqrt(2.)*qCM.Pz());
        Ur[52]=((3*qCM.E() + 7*qCM.Pz() + 10*qpCM.E())*TMath::Power(Q,4) + 
     8*qCM.E()*(4*qCM.E() + qCM.Pz())*
      (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
     TMath::Power(Q,2)*(8*(4*qCM.E() + qCM.Pz())*
         (qCM.E() - qpCM.E())*qpCM.E() + 
        (32*qCM.E() + 6*qCM.Pz())*TMath::Power(qpPerp,2) + 
        (3*qCM.E() + 7*qCM.Pz() + 16*qpCM.E())*t))/
   (32.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,4));
        Ur[53]=(3*TMath::Power(Q,2)*((qCM.E() + qCM.Pz() - 2*qpCM.E())*TMath::Power(Q,2) + 
       2*qCM.Pz()*TMath::Power(qpPerp,2) + (qCM.E() + qCM.Pz())*t))/
   (32.*TMath::Sqrt(2.)*qCM.Pz());
        Ur[54]=((-qCM.E() + 3*qCM.Pz() + 2*qpCM.E())*TMath::Power(Q,4) - 
     TMath::Power(Q,2)*(2*qCM.Pz()*
         (4*qpCM.E()*(-qCM.E() + qpCM.E()) + TMath::Power(qpPerp,2)) + 
        (qCM.E() - 3*qCM.Pz())*t) + 
     8*qCM.E()*qCM.Pz()*(qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t)
     )/(32.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,4));
        Ur[55]=(qpPerp*((qCM.E() - qCM.Pz() - 2*qpCM.E())*TMath::Power(Q,2) + 
       (qCM.E() + 3*qCM.Pz())*t))/(8.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2));
        Ur[56]=-(qpPerp*((qCM.E() + qCM.Pz() - 2*qpCM.E())*TMath::Power(Q,2) + 
        (qCM.E() + qCM.Pz())*t))/(8.*TMath::Sqrt(2.)*qCM.Pz());
        Ur[57]=(qpPerp*(3*(qCM.E() - qCM.Pz() - 2*qpCM.E())*TMath::Power(Q,2) + 
       (3*qCM.E() + qCM.Pz())*t))/(8.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2));
        Ur[58]=-((qCM.E() - qpCM.E())*TMath::Power(Q,4) + 
      2*qCM.E()*(qCM.E() + qCM.Pz())*
       (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
      TMath::Power(Q,2)*(2*(qCM.E() + qCM.Pz())*(qCM.E() - qpCM.E()\
             )*qpCM.E() + 
         2*(qCM.E() - qCM.Pz())*TMath::Power(qpPerp,2) + 
         (qCM.E() + qpCM.E())*t))/(4.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,4));
        Ur[59]=-(qpPerp*(TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 
        4*s*t)*(11*TMath::Power(Q,4) + 
        4*qCM.E()*(qCM.E() - 4*qCM.Pz())*t + 
        TMath::Power(Q,2)*(2*(6*TMath::Power(qCM.E(),2) + 3*qCM.Pz()*qpCM.E() - 
              qCM.E()*(8*qCM.Pz() + qpCM.E())) + 3*t)))/
   (128.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[60]=-(qpPerp*(TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 
        4*s*t)*(7*TMath::Power(Q,4) + 4*TMath::Power(qCM.E(),2)*t - 
        TMath::Power(Q,2)*(2*(-6*TMath::Power(qCM.E(),2) + 5*qCM.E()*qpCM.E() + 
              qCM.Pz()*qpCM.E()) + t)))/
   (128.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[61]=-((TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t)*
      (TMath::Power(Q,4) + 2*qCM.E()*(3*qCM.E() - qCM.Pz())*
         TMath::Power(qpPerp,2) + (qCM.E() + qCM.Pz())*
         (qCM.E() + qpCM.E())*t + 
        TMath::Power(Q,2)*((qCM.E() - qpCM.E())*
            (qCM.E() + qCM.Pz() + 2*qpCM.E()) + 
           6*TMath::Power(qpPerp,2) + t)))/
   (32.*qCM.Pz()*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[62]=((TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t)*
     (TMath::Power(Q,4) - 2*qCM.E()*(qCM.E() - 3*qCM.Pz())*
        TMath::Power(qpPerp,2) + (qCM.E() + qCM.Pz())*
        (qCM.E() + qpCM.E())*t + 
       TMath::Power(Q,2)*((qCM.E() - qpCM.E())*
           (qCM.E() + qCM.Pz() + 2*qpCM.E()) - 
          2*TMath::Power(qpPerp,2) + t)))/
   (32.*qCM.Pz()*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[63]=(qpPerp*(TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 
       4*s*t)*(TMath::Power(Q,4) + qCM.E()*(qCM.E() - qCM.Pz())*t - 
       TMath::Power(Q,2)*(-3*TMath::Power(qCM.E(),2) - 2*qCM.Pz()*qpCM.E() + 
          qCM.E()*(qCM.Pz() + 4*qpCM.E()) + t)))/
   (16.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[64]=-(TMath::Sqrt(s)*(13*TMath::Power(Q,6) + TMath::Power(Q,4)*
         (13*TMath::Power(qCM.E(),2) + 3*qCM.E()*qCM.Pz() + 
           45*qCM.E()*qpCM.E() - 7*qCM.Pz()*qpCM.E() - 
           42*TMath::Power(qpCM.E(),2) + 52*TMath::Power(qpPerp,2) + 13*t) + 
        8*qCM.E()*(4*qCM.E() - qCM.Pz())*
         (qCM.E() + qpCM.E())*
         (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
        TMath::Power(Q,2)*(8*(4*qCM.E() - qCM.Pz())*qpCM.E()*
            (TMath::Power(qCM.E(),2) - TMath::Power(qpCM.E(),2)) + 
           2*(42*TMath::Power(qCM.E(),2) - 3*qCM.Pz()*qpCM.E() + 
              qCM.E()*(-5*qCM.Pz() + 16*qpCM.E()))*TMath::Power(qpPerp,2)
             + (13*TMath::Power(qCM.E(),2) + 
              3*qCM.E()*(qCM.Pz() + 15*qpCM.E()) + 
              qpCM.E()*(-qCM.Pz() + 16*qpCM.E()))*t)))/
   (16.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[65]=-(TMath::Sqrt(s)*(TMath::Power(Q,6) + TMath::Power(Q,4)*
         (TMath::Power(qCM.E(),2) - qCM.E()*qCM.Pz() + 
           qCM.E()*qpCM.E() - 3*qCM.Pz()*qpCM.E() - 
           2*TMath::Power(qpCM.E(),2) + 4*TMath::Power(qpPerp,2) + t) - 
        8*qCM.E()*qCM.Pz()*(qCM.E() + qpCM.E())*
         (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
        TMath::Power(Q,2)*(8*qCM.Pz()*qpCM.E()*
            (-TMath::Power(qCM.E(),2) + TMath::Power(qpCM.E(),2)) + 
           2*(2*TMath::Power(qCM.E(),2) - 9*qCM.E()*qCM.Pz() + 
              qCM.Pz()*qpCM.E())*TMath::Power(qpPerp,2) + 
           (TMath::Power(qCM.E(),2) - 5*qCM.Pz()*qpCM.E() + 
              qCM.E()*(-qCM.Pz() + qpCM.E()))*t)))/
   (16.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[66]=(qpPerp*TMath::Sqrt(s)*((3*TMath::Power(qCM.E(),2) - 3*qCM.Pz()*qpCM.E() + 
          qCM.E()*(-qCM.Pz() + qpCM.E()))*t + 
       TMath::Power(Q,2)*(-((qCM.E() + qCM.Pz() - 2*qpCM.E())*
             (qCM.E() - qpCM.E())) + 4*t)))/
   (4.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[67]=(qpPerp*TMath::Sqrt(s)*((TMath::Power(qCM.E(),2) - 
          3*qCM.E()*(qCM.Pz() - qpCM.E()) - 
          qCM.Pz()*qpCM.E())*t + 
       TMath::Power(Q,2)*(-3*(qCM.E() + qCM.Pz() - 2*qpCM.E())*
           (qCM.E() - qpCM.E()) + 4*t)))/
   (4.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[68]=(TMath::Sqrt(s)*(-((qCM.E()*(qCM.Pz() - 2*qpCM.E()) + 
            TMath::Power(qpCM.E(),2))*TMath::Power(Q,4)) + 
       2*qCM.E()*(qCM.E() - qCM.Pz())*
        (qCM.E() + qpCM.E())*
        (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
       TMath::Power(Q,2)*(2*(qCM.E() - qCM.Pz())*qpCM.E()*
           (TMath::Power(qCM.E(),2) - TMath::Power(qpCM.E(),2)) + 
          2*(TMath::Power(qCM.E(),2) + qCM.Pz()*qpCM.E() + 
             qCM.E()*(-3*qCM.Pz() + qpCM.E()))*TMath::Power(qpPerp,2) + 
          (-(qCM.E()*qCM.Pz()) + 2*qCM.E()*qpCM.E() - 
             2*qCM.Pz()*qpCM.E() + TMath::Power(qpCM.E(),2))*t)))/
   (2.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[69]=(qpPerp*(3*(qCM.E() + qCM.Pz() - 2*qpCM.E())*TMath::Power(Q,2) + 
       (3*qCM.E() - 4*qCM.Pz())*t))/(8.*TMath::Sqrt(2.)*qCM.Pz());
        Ur[70]=((13*qCM.E() + 3*qCM.Pz() - 10*qpCM.E())*TMath::Power(Q,4) + 
     8*qCM.E()*(4*qCM.E() - qCM.Pz())*
      (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
     TMath::Power(Q,2)*(8*(4*qCM.E() - qCM.Pz())*
         (qCM.E() - qpCM.E())*qpCM.E() + 
        2*(16*qCM.E() + 7*qCM.Pz())*TMath::Power(qpPerp,2) + 
        (13*qCM.E() + 3*qCM.Pz() + 16*qpCM.E())*t))/
   (32.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2));
        Ur[71]=((qCM.E() - qCM.Pz() - 2*qpCM.E())*TMath::Power(Q,4) + 
     TMath::Power(Q,2)*(2*qCM.Pz()*
         (4*qpCM.E()*(-qCM.E() + qpCM.E()) + 3*TMath::Power(qpPerp,2)) + 
        (qCM.E() - qCM.Pz())*t) - 
     8*qCM.E()*qCM.Pz()*(qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t)
     )/(32.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2));
        Ur[72]=((qCM.Pz() - qpCM.E())*TMath::Power(Q,4) + 
     TMath::Power(Q,2)*(-2*(qCM.E() - qCM.Pz())*(qCM.E() - qpCM.E()\
            )*qpCM.E() - 2*qCM.E()*TMath::Power(qpPerp,2) + 
        (qCM.Pz() - qpCM.E())*t) - 
     2*qCM.E()*(qCM.E() - qCM.Pz())*
      (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t))/
   (4.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2));
        Ur[73]=(qpCM.E()*qpPerp*(TMath::Power(Q,4) + 
       2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t))/
   (16.*TMath::Sqrt(2.)*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[74]=((TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t)*
     (3*TMath::Power(Q,4) + qCM.E()*
        (8*qCM.E()*TMath::Power(qpPerp,2) + 3*(qCM.E() + qpCM.E())*t) + 
       TMath::Power(Q,2)*(8*TMath::Power(qpPerp,2) + 
          3*(TMath::Power(qCM.E(),2) + qCM.E()*qpCM.E() - 
             2*TMath::Power(qpCM.E(),2) + t))))/
   (16.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*
     TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[75]=-((TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t)*
      (TMath::Power(Q,4) + qCM.E()*(qCM.E() + qpCM.E())*t + 
        TMath::Power(Q,2)*(TMath::Power(qCM.E(),2) + qCM.E()*qpCM.E() - 
           2*TMath::Power(qpCM.E(),2) + t)))/
   (16.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*
     TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[76]=-(qpPerp*((qCM.E() - qpCM.E())*TMath::Power(Q,2) + qCM.E()*t)*
      (TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t)
      )/(8.*TMath::Sqrt(2.)*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[77]=(qCM.Pz()*TMath::Power(qpPerp,2)*(TMath::Power(Q,4) + 
       2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t))/
   (4.*TMath::Sqrt(2.)*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[78]=-(TMath::Sqrt(s)*(TMath::Power(Q,4) + TMath::Power(Q,2)*
         (TMath::Power(qCM.E(),2) + qCM.E()*qpCM.E() - 
           2*TMath::Power(qpCM.E(),2) + 4*TMath::Power(qpPerp,2) + t) + 
        qCM.E()*(4*qCM.E()*TMath::Power(qpPerp,2) + 
           (qCM.E() + qpCM.E())*t)))/
   (8.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[79]=-(qpPerp*TMath::Sqrt(s)*((qCM.E() - qpCM.E())*TMath::Power(Q,2) + 
        (qCM.E() - 2*qpCM.E())*t))/
   (2.*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[80]=(qpPerp*TMath::Sqrt(s)*((qCM.E() - qpCM.E())*TMath::Power(Q,2) + qCM.E()*t))/
   (2.*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[81]=(TMath::Sqrt(s)*(TMath::Power(Q,6) + TMath::Power(Q,4)*
        (TMath::Power(qCM.E(),2) + 5*qCM.E()*qpCM.E() - 
          4*TMath::Power(qpCM.E(),2) + 4*TMath::Power(qpPerp,2) + t) + 
       4*TMath::Power(qCM.E(),2)*(qCM.E() + qpCM.E())*
        (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
       TMath::Power(Q,2)*(4*TMath::Power(qCM.E(),3)*qpCM.E() - 
          4*qCM.E()*TMath::Power(qpCM.E(),3) + 
          4*qCM.E()*(2*qCM.E() + qpCM.E())*TMath::Power(qpPerp,2) + 
          TMath::Power(qCM.E(),2)*t + 5*qCM.E()*qpCM.E()*t + 
          2*TMath::Power(qpCM.E(),2)*t)))/
   (4.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[82]=(qpPerp*TMath::Sqrt(s)*((qCM.E() - qpCM.E())*TMath::Power(Q,2) + 
       (qCM.E() + qpCM.E())*t))/
   (2.*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[83]=-(TMath::Power(Q,4) + TMath::Power(Q,2)*(4*(qCM.E() - qpCM.E())*qpCM.E() - 
         2*TMath::Power(qpPerp,2) + t) + 
      4*qCM.E()*(qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t))/
   (16.*TMath::Power(Q,2));
        Ur[84]=(TMath::Power(Q,2) + 2*TMath::Power(qpPerp,2) + t)/(16.*TMath::Power(Q,2));
        Ur[85]=(qpPerp*((qCM.E() - 2*qpCM.E())*TMath::Power(Q,2) + qCM.E()*t))/
   (4.*qCM.Pz()*TMath::Power(Q,2));
        Ur[86]=-(qpPerp*(TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 
        4*s*t)*(11*TMath::Power(Q,4) + 
        4*qCM.E()*(qCM.E() + 4*qCM.Pz())*t + 
        TMath::Power(Q,2)*(12*TMath::Power(qCM.E(),2) + 16*qCM.E()*qCM.Pz() - 
           2*qCM.E()*qpCM.E() - 6*qCM.Pz()*qpCM.E() + 3*t)))/
   (128.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[87]=-(qpPerp*(TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 
        4*s*t)*(7*TMath::Power(Q,4) + 4*TMath::Power(qCM.E(),2)*t - 
        TMath::Power(Q,2)*(-2*(6*TMath::Power(qCM.E(),2) - 5*qCM.E()*qpCM.E() + 
              qCM.Pz()*qpCM.E()) + t)))/
   (128.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[88]=((TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t)*
     (TMath::Power(Q,4) + 2*qCM.E()*(3*qCM.E() + qCM.Pz())*
        TMath::Power(qpPerp,2) + (qCM.E() - qCM.Pz())*
        (qCM.E() + qpCM.E())*t + 
       TMath::Power(Q,2)*((qCM.E() - qpCM.E())*
           (qCM.E() - qCM.Pz() + 2*qpCM.E()) + 
          6*TMath::Power(qpPerp,2) + t)))/
   (32.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[89]=-((TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 4*s*t)*
      (TMath::Power(Q,4) - 2*qCM.E()*(qCM.E() + 3*qCM.Pz())*
         TMath::Power(qpPerp,2) + (qCM.E() - qCM.Pz())*
         (qCM.E() + qpCM.E())*t + 
        TMath::Power(Q,2)*((qCM.E() - qpCM.E())*
            (qCM.E() - qCM.Pz() + 2*qpCM.E()) - 
           2*TMath::Power(qpPerp,2) + t)))/
   (32.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[90]=(qpPerp*(TMath::Power(Q,4) + 2*(qCM.E() - qpCM.E())*TMath::Power(Q,2)*TMath::Sqrt(s) - 
       4*s*t)*(TMath::Power(Q,4) + qCM.E()*(qCM.E() + qCM.Pz())*t - 
       TMath::Power(Q,2)*(-3*TMath::Power(qCM.E(),2) - qCM.E()*qCM.Pz() + 
          4*qCM.E()*qpCM.E() + 2*qCM.Pz()*qpCM.E() + t)))/
   (16.*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(s)*TMath::Sqrt(-t)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[91]=(TMath::Sqrt(s)*(13*TMath::Power(Q,6) + TMath::Power(Q,4)*
        (13*TMath::Power(qCM.E(),2) - 3*qCM.E()*qCM.Pz() + 
          45*qCM.E()*qpCM.E() + 7*qCM.Pz()*qpCM.E() - 
          42*TMath::Power(qpCM.E(),2) + 52*TMath::Power(qpPerp,2) + 13*t) + 
       8*qCM.E()*(4*qCM.E() + qCM.Pz())*
        (qCM.E() + qpCM.E())*
        (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
       TMath::Power(Q,2)*(8*(4*qCM.E() + qCM.Pz())*qpCM.E()*
           (TMath::Power(qCM.E(),2) - TMath::Power(qpCM.E(),2)) + 
          2*(42*TMath::Power(qCM.E(),2) + 3*qCM.Pz()*qpCM.E() + 
             qCM.E()*(5*qCM.Pz() + 16*qpCM.E()))*TMath::Power(qpPerp,2)\
           + (13*TMath::Power(qCM.E(),2) - 
             3*qCM.E()*(qCM.Pz() - 15*qpCM.E()) + 
             qpCM.E()*(qCM.Pz() + 16*qpCM.E()))*t)))/
   (16.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,4)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[92]=(TMath::Sqrt(s)*(TMath::Power(Q,6) + TMath::Power(Q,4)*
        (TMath::Power(qCM.E(),2) + qCM.E()*qCM.Pz() + 
          qCM.E()*qpCM.E() + 3*qCM.Pz()*qpCM.E() - 
          2*TMath::Power(qpCM.E(),2) + 4*TMath::Power(qpPerp,2) + t) + 
       8*qCM.E()*qCM.Pz()*(qCM.E() + qpCM.E())*
        (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
       TMath::Power(Q,2)*(8*qCM.Pz()*qpCM.E()*
           (TMath::Power(qCM.E(),2) - TMath::Power(qpCM.E(),2)) + 
          2*(2*TMath::Power(qCM.E(),2) + 9*qCM.E()*qCM.Pz() - 
             qCM.Pz()*qpCM.E())*TMath::Power(qpPerp,2) + 
          (TMath::Power(qCM.E(),2) + 5*qCM.Pz()*qpCM.E() + 
             qCM.E()*(qCM.Pz() + qpCM.E()))*t)))/
   (16.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,4)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[93]=(qpPerp*TMath::Sqrt(s)*((3*TMath::Power(qCM.E(),2) + 3*qCM.Pz()*qpCM.E() + 
          qCM.E()*(qCM.Pz() + qpCM.E()))*t + 
       TMath::Power(Q,2)*(-((qCM.E() - qCM.Pz() - 2*qpCM.E())*
             (qCM.E() - qpCM.E())) + 4*t)))/
   (4.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[94]=(qpPerp*TMath::Sqrt(s)*((TMath::Power(qCM.E(),2) + qCM.Pz()*qpCM.E() + 
          3*qCM.E()*(qCM.Pz() + qpCM.E()))*t + 
       TMath::Power(Q,2)*(-3*(qCM.E() - qCM.Pz() - 2*qpCM.E())*
           (qCM.E() - qpCM.E()) + 4*t)))/
   (4.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[95]=-(TMath::Sqrt(s)*((-TMath::Power(qpCM.E(),2) + 
           qCM.E()*(qCM.Pz() + 2*qpCM.E()))*TMath::Power(Q,4) + 
        2*qCM.E()*(qCM.E() + qCM.Pz())*
         (qCM.E() + qpCM.E())*
         (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
        TMath::Power(Q,2)*(2*(qCM.E() + qCM.Pz())*qpCM.E()*
            (TMath::Power(qCM.E(),2) - TMath::Power(qpCM.E(),2)) + 
           2*(TMath::Power(qCM.E(),2) - qCM.Pz()*qpCM.E() + 
              qCM.E()*(3*qCM.Pz() + qpCM.E()))*TMath::Power(qpPerp,2) + 
           (qpCM.E()*(2*qCM.Pz() + qpCM.E()) + 
              qCM.E()*(qCM.Pz() + 2*qpCM.E()))*t)))/
   (2.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,4)*TMath::Sqrt(TMath::Power(Q,4) - 4*s*t));
        Ur[96]=(qpPerp*((-3*qCM.E() + 3*qCM.Pz() + 6*qpCM.E())*TMath::Power(Q,2) - 
       (3*qCM.E() + 4*qCM.Pz())*t))/
   (8.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,2));
        Ur[97]=((13*qCM.E() - 3*qCM.Pz() - 10*qpCM.E())*TMath::Power(Q,4) + 
     8*qCM.E()*(4*qCM.E() + qCM.Pz())*
      (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
     TMath::Power(Q,2)*(8*(4*qCM.E() + qCM.Pz())*
         (qCM.E() - qpCM.E())*qpCM.E() + 
        2*(16*qCM.E() - 7*qCM.Pz())*TMath::Power(qpPerp,2) + 
        (13*qCM.E() - 3*qCM.Pz() + 16*qpCM.E())*t))/
   (32.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,4));
        Ur[98]=((qCM.E() + qCM.Pz() - 2*qpCM.E())*TMath::Power(Q,4) + 
     TMath::Power(Q,2)*(-2*qCM.Pz()*
         (4*qpCM.E()*(-qCM.E() + qpCM.E()) + 3*TMath::Power(qpPerp,2)) + 
        (qCM.E() + qCM.Pz())*t) + 
     8*qCM.E()*qCM.Pz()*(qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t)
     )/(32.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,4));
        Ur[99]=-((qCM.Pz() + qpCM.E())*TMath::Power(Q,4) + 
      2*qCM.E()*(qCM.E() + qCM.Pz())*
       (qCM.E()*TMath::Power(qpPerp,2) + qpCM.E()*t) + 
      TMath::Power(Q,2)*(2*(qCM.E() + qCM.Pz())*(qCM.E() - qpCM.E()\
             )*qpCM.E() + 2*qCM.E()*TMath::Power(qpPerp,2) + 
         (qCM.Pz() + qpCM.E())*t))/
   (4.*TMath::Sqrt(2.)*qCM.Pz()*TMath::Power(Q,4));
	
	
/*-------------------------------- VCS cross sections ----------------------------------*/

        SigmaVCSPol0[0]=(3*TMath::Power(IMvcs[0][0],2) - 2*TMath::Power(IMvcs[0][1],2) + 
     3*TMath::Power(IMvcs[0][2],2) + 3*TMath::Power(IMvcs[1][0],2) - 
     2*TMath::Power(IMvcs[1][1],2) + 3*TMath::Power(IMvcs[1][2],2) + 
     3*TMath::Power(IMvcs[2][0],2) - 2*TMath::Power(IMvcs[2][1],2) + 
     3*TMath::Power(IMvcs[2][2],2) + 3*TMath::Power(IMvcs[3][0],2) - 
     2*TMath::Power(IMvcs[3][1],2) + 3*TMath::Power(IMvcs[3][2],2) + 
     3*TMath::Power(RMvcs[0][0],2) - 2*TMath::Power(RMvcs[0][1],2) + 
     3*TMath::Power(RMvcs[0][2],2) + 3*TMath::Power(RMvcs[1][0],2) - 
     2*TMath::Power(RMvcs[1][1],2) + 3*TMath::Power(RMvcs[1][2],2) + 
     3*TMath::Power(RMvcs[2][0],2) - 2*TMath::Power(RMvcs[2][1],2) + 
     3*TMath::Power(RMvcs[2][2],2) + 3*TMath::Power(RMvcs[3][0],2) - 
     2*TMath::Power(RMvcs[3][1],2) + 3*TMath::Power(RMvcs[3][2],2))/4.;
        SigmaVCSPol0[1]=(TMath::Power(IMvcs[0][0],2) + 2*TMath::Power(IMvcs[0][1],2) + 
     TMath::Power(IMvcs[0][2],2) + TMath::Power(IMvcs[1][0],2) + 
     2*TMath::Power(IMvcs[1][1],2) + TMath::Power(IMvcs[1][2],2) + 
     TMath::Power(IMvcs[2][0],2) + 2*TMath::Power(IMvcs[2][1],2) + 
     TMath::Power(IMvcs[2][2],2) + TMath::Power(IMvcs[3][0],2) + 
     2*TMath::Power(IMvcs[3][1],2) + TMath::Power(IMvcs[3][2],2) + 
     TMath::Power(RMvcs[0][0],2) + 2*TMath::Power(RMvcs[0][1],2) + 
     TMath::Power(RMvcs[0][2],2) + TMath::Power(RMvcs[1][0],2) + 
     2*TMath::Power(RMvcs[1][1],2) + TMath::Power(RMvcs[1][2],2) + 
     TMath::Power(RMvcs[2][0],2) + 2*TMath::Power(RMvcs[2][1],2) + 
     TMath::Power(RMvcs[2][2],2) + TMath::Power(RMvcs[3][0],2) + 
     2*TMath::Power(RMvcs[3][1],2) + TMath::Power(RMvcs[3][2],2))/4.;
        SigmaVCSPol0[2]=(IMvcs[0][1]*(IMvcs[0][0] - IMvcs[0][2]) + 
     IMvcs[1][1]*(IMvcs[1][0] - IMvcs[1][2]) + 
     IMvcs[2][1]*(IMvcs[2][0] - IMvcs[2][2]) + 
     IMvcs[3][1]*(IMvcs[3][0] - IMvcs[3][2]) + 
     RMvcs[0][1]*(RMvcs[0][0] - RMvcs[0][2]) + 
     RMvcs[1][1]*(RMvcs[1][0] - RMvcs[1][2]) + 
     RMvcs[2][1]*(RMvcs[2][0] - RMvcs[2][2]) + 
     RMvcs[3][1]*(RMvcs[3][0] - RMvcs[3][2]))/TMath::Sqrt(2.);
        SigmaVCSPol0[3]=(IMvcs[0][0]*IMvcs[0][2] + IMvcs[1][0]*IMvcs[1][2] + 
     IMvcs[2][0]*IMvcs[2][2] + IMvcs[3][0]*IMvcs[3][2] + 
     RMvcs[0][0]*RMvcs[0][2] + RMvcs[1][0]*RMvcs[1][2] + 
     RMvcs[2][0]*RMvcs[2][2] + RMvcs[3][0]*RMvcs[3][2])/2.;
        SigmaVCSPol0[4]=TMath::Sqrt(2.)*((IMvcs[0][0] - IMvcs[0][2])*RMvcs[0][1] + 
     IMvcs[0][1]*(-RMvcs[0][0] + RMvcs[0][2]) + 
     (IMvcs[1][0] - IMvcs[1][2])*RMvcs[1][1] + 
     IMvcs[1][1]*(-RMvcs[1][0] + RMvcs[1][2]) + 
     (IMvcs[2][0] - IMvcs[2][2])*RMvcs[2][1] + 
     IMvcs[2][1]*(-RMvcs[2][0] + RMvcs[2][2]) + 
     (IMvcs[3][0] - IMvcs[3][2])*RMvcs[3][1] + 
     IMvcs[3][1]*(-RMvcs[3][0] + RMvcs[3][2]));
        SigmaVCSPolX[0]=2*(IMvcs[0][0]*IMvcs[1][0] - IMvcs[0][2]*IMvcs[1][2] - 
     IMvcs[2][0]*IMvcs[3][0] + IMvcs[2][2]*IMvcs[3][2] + 
     RMvcs[0][0]*RMvcs[1][0] - RMvcs[0][2]*RMvcs[1][2] - 
     RMvcs[2][0]*RMvcs[3][0] + RMvcs[2][2]*RMvcs[3][2]);
        SigmaVCSPolX[1]=TMath::Sqrt(2.)*((IMvcs[0][0] + IMvcs[0][2])*IMvcs[1][1] + 
     IMvcs[0][1]*(IMvcs[1][0] + IMvcs[1][2]) - 
     (IMvcs[2][0] + IMvcs[2][2])*IMvcs[3][1] - 
     IMvcs[2][1]*(IMvcs[3][0] + IMvcs[3][2]) + 
     (RMvcs[0][0] + RMvcs[0][2])*RMvcs[1][1] + 
     RMvcs[0][1]*(RMvcs[1][0] + RMvcs[1][2]) - 
     (RMvcs[2][0] + RMvcs[2][2])*RMvcs[3][1] - 
     RMvcs[2][1]*(RMvcs[3][0] + RMvcs[3][2]));
        SigmaVCSPolX[2]=((IMvcs[1][0] + IMvcs[1][2])*RMvcs[0][1] - 
     IMvcs[1][1]*(RMvcs[0][0] + RMvcs[0][2]) + 
     (IMvcs[0][0] + IMvcs[0][2])*RMvcs[1][1] - 
     IMvcs[0][1]*(RMvcs[1][0] + RMvcs[1][2]) - 
     (IMvcs[3][0] + IMvcs[3][2])*RMvcs[2][1] + 
     IMvcs[3][1]*(RMvcs[2][0] + RMvcs[2][2]) - 
     (IMvcs[2][0] + IMvcs[2][2])*RMvcs[3][1] + 
     IMvcs[2][1]*(RMvcs[3][0] + RMvcs[3][2]))/TMath::Sqrt(2.);
        SigmaVCSPolX[3]=(-(IMvcs[1][2]*RMvcs[0][0]) + IMvcs[1][0]*RMvcs[0][2] - 
     IMvcs[0][2]*RMvcs[1][0] + IMvcs[0][0]*RMvcs[1][2] + 
     IMvcs[3][2]*RMvcs[2][0] - IMvcs[3][0]*RMvcs[2][2] + 
     IMvcs[2][2]*RMvcs[3][0] - IMvcs[2][0]*RMvcs[3][2])/2.;
        SigmaVCSPolY[0]=TMath::Sqrt(2.)*((IMvcs[0][0] - IMvcs[0][2])*IMvcs[2][1] + 
     IMvcs[0][1]*(-IMvcs[2][0] + IMvcs[2][2]) + 
     (IMvcs[1][0] - IMvcs[1][2])*IMvcs[3][1] + 
     IMvcs[1][1]*(-IMvcs[3][0] + IMvcs[3][2]) + 
     (RMvcs[0][0] - RMvcs[0][2])*RMvcs[2][1] + 
     RMvcs[0][1]*(-RMvcs[2][0] + RMvcs[2][2]) + 
     (RMvcs[1][0] - RMvcs[1][2])*RMvcs[3][1] + 
     RMvcs[1][1]*(-RMvcs[3][0] + RMvcs[3][2]));
        SigmaVCSPolY[1]=(3*IMvcs[2][0]*RMvcs[0][0] - 2*IMvcs[2][1]*RMvcs[0][1] + 
     3*IMvcs[2][2]*RMvcs[0][2] + 3*IMvcs[3][0]*RMvcs[1][0] - 
     2*IMvcs[3][1]*RMvcs[1][1] + 3*IMvcs[3][2]*RMvcs[1][2] - 
     3*IMvcs[0][0]*RMvcs[2][0] + 2*IMvcs[0][1]*RMvcs[2][1] - 
     3*IMvcs[0][2]*RMvcs[2][2] - 3*IMvcs[1][0]*RMvcs[3][0] + 
     2*IMvcs[1][1]*RMvcs[3][1] - 3*IMvcs[1][2]*RMvcs[3][2])/2.;
        SigmaVCSPolY[2]=(IMvcs[2][0]*RMvcs[0][0] + 2*IMvcs[2][1]*RMvcs[0][1] + 
     IMvcs[2][2]*RMvcs[0][2] + IMvcs[3][0]*RMvcs[1][0] + 
     2*IMvcs[3][1]*RMvcs[1][1] + IMvcs[3][2]*RMvcs[1][2] - 
     IMvcs[0][0]*RMvcs[2][0] - 2*IMvcs[0][1]*RMvcs[2][1] - 
     IMvcs[0][2]*RMvcs[2][2] - IMvcs[1][0]*RMvcs[3][0] - 
     2*IMvcs[1][1]*RMvcs[3][1] - IMvcs[1][2]*RMvcs[3][2])/2.;
        SigmaVCSPolY[3]=((IMvcs[2][0] - IMvcs[2][2])*RMvcs[0][1] + 
     IMvcs[2][1]*(RMvcs[0][0] - RMvcs[0][2]) + 
     (IMvcs[3][0] - IMvcs[3][2])*RMvcs[1][1] + 
     IMvcs[3][1]*(RMvcs[1][0] - RMvcs[1][2]) + 
     (-IMvcs[0][0] + IMvcs[0][2])*RMvcs[2][1] + 
     IMvcs[0][1]*(-RMvcs[2][0] + RMvcs[2][2]) + 
     (-IMvcs[1][0] + IMvcs[1][2])*RMvcs[3][1] + 
     IMvcs[1][1]*(-RMvcs[3][0] + RMvcs[3][2]))/TMath::Sqrt(2.);
        SigmaVCSPolY[4]=(IMvcs[2][2]*RMvcs[0][0] + IMvcs[2][0]*RMvcs[0][2] + 
     IMvcs[3][2]*RMvcs[1][0] + IMvcs[3][0]*RMvcs[1][2] - 
     IMvcs[0][2]*RMvcs[2][0] - IMvcs[0][0]*RMvcs[2][2] - 
     IMvcs[1][2]*RMvcs[3][0] - IMvcs[1][0]*RMvcs[3][2])/2.;
        SigmaVCSPolZ[0]=2*(IMvcs[1][0]*IMvcs[2][0] - IMvcs[1][2]*IMvcs[2][2] + 
     IMvcs[0][0]*IMvcs[3][0] - IMvcs[0][2]*IMvcs[3][2] + 
     RMvcs[1][0]*RMvcs[2][0] - RMvcs[1][2]*RMvcs[2][2] + 
     RMvcs[0][0]*RMvcs[3][0] - RMvcs[0][2]*RMvcs[3][2]);
        SigmaVCSPolZ[1]=TMath::Sqrt(2.)*((IMvcs[1][0] + IMvcs[1][2])*IMvcs[2][1] + 
     IMvcs[1][1]*(IMvcs[2][0] + IMvcs[2][2]) + 
     (IMvcs[0][0] + IMvcs[0][2])*IMvcs[3][1] + 
     IMvcs[0][1]*(IMvcs[3][0] + IMvcs[3][2]) + 
     (RMvcs[1][0] + RMvcs[1][2])*RMvcs[2][1] + 
     RMvcs[1][1]*(RMvcs[2][0] + RMvcs[2][2]) + 
     (RMvcs[0][0] + RMvcs[0][2])*RMvcs[3][1] + 
     RMvcs[0][1]*(RMvcs[3][0] + RMvcs[3][2]));
        SigmaVCSPolZ[2]=((IMvcs[3][0] + IMvcs[3][2])*RMvcs[0][1] - 
     IMvcs[3][1]*(RMvcs[0][0] + RMvcs[0][2]) + 
     (IMvcs[2][0] + IMvcs[2][2])*RMvcs[1][1] - 
     IMvcs[2][1]*(RMvcs[1][0] + RMvcs[1][2]) + 
     (IMvcs[1][0] + IMvcs[1][2])*RMvcs[2][1] - 
     IMvcs[1][1]*(RMvcs[2][0] + RMvcs[2][2]) + 
     (IMvcs[0][0] + IMvcs[0][2])*RMvcs[3][1] - 
     IMvcs[0][1]*(RMvcs[3][0] + RMvcs[3][2]))/TMath::Sqrt(2.);
        SigmaVCSPolZ[3]=(-(IMvcs[3][2]*RMvcs[0][0]) + IMvcs[3][0]*RMvcs[0][2] - 
     IMvcs[2][2]*RMvcs[1][0] + IMvcs[2][0]*RMvcs[1][2] - 
     IMvcs[1][2]*RMvcs[2][0] + IMvcs[1][0]*RMvcs[2][2] - 
     IMvcs[0][2]*RMvcs[3][0] + IMvcs[0][0]*RMvcs[3][2])/2.;
		
	
/*--------------------------- Interference cross sections ------------------------------*/

		if ( InitExactBHCrossSections == kFALSE )
		{
			TGVKelly::MakeExactBHCrossSections();
			
			
			// Flag
			
			InitExactBHCrossSections = kTRUE;
		} // end if InitExactBHCrossSections
		

        SigmaIPol0[0]=Jem[0][1]*(-2*RMvcs[0][0]*Ur[0]*TMath::Power(Q,2) + 
      2*RMvcs[0][1]*Ur[24]*TMath::Power(Q,3) - 
      2*RMvcs[0][2]*Ur[39]*TMath::Power(Q,4)) + 
   Jem[2][1]*(-2*RMvcs[2][0]*Ur[0]*TMath::Power(Q,2) + 
      2*RMvcs[2][1]*Ur[24]*TMath::Power(Q,3) - 
      2*RMvcs[2][2]*Ur[39]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*RMvcs[0][0]*Ur[8]*TMath::Power(Q,2) + 
      2*RMvcs[0][1]*Ur[29]*TMath::Power(Q,3) - 
      2*RMvcs[0][2]*Ur[45]*TMath::Power(Q,4)) + 
   Jem[2][2]*(-2*RMvcs[2][0]*Ur[8]*TMath::Power(Q,2) + 
      2*RMvcs[2][1]*Ur[29]*TMath::Power(Q,3) - 
      2*RMvcs[2][2]*Ur[45]*TMath::Power(Q,4)) + 
   Jem[1][2]*(-2*RMvcs[1][0]*Ur[69]*TMath::Power(Q,2) + 
      2*RMvcs[1][1]*Ur[83]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[96]*TMath::Power(Q,4)) + 
   Jem[3][2]*(-2*RMvcs[3][0]*Ur[69]*TMath::Power(Q,2) + 
      2*RMvcs[3][1]*Ur[83]*TMath::Power(Q,3) - 
      2*RMvcs[3][2]*Ur[96]*TMath::Power(Q,4));
        SigmaIPol0[1]=Jem[0][1]*(-2*RMvcs[0][0]*Ur[2]*TMath::Power(Q,2) - 
      2*RMvcs[0][1]*Ur[24]*TMath::Power(Q,3) - 
      2*RMvcs[0][2]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[2][1]*(-2*RMvcs[2][0]*Ur[2]*TMath::Power(Q,2) - 
      2*RMvcs[2][1]*Ur[24]*TMath::Power(Q,3) - 
      2*RMvcs[2][2]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[1][2]*(-2*RMvcs[1][0]*Ur[51]*TMath::Power(Q,2) - 
      2*RMvcs[1][1]*Ur[83]*TMath::Power(Q,3) + 
      2*RMvcs[1][2]*Ur[16]*TMath::Power(Q,4)) + 
   Jem[3][2]*(-2*RMvcs[3][0]*Ur[51]*TMath::Power(Q,2) - 
      2*RMvcs[3][1]*Ur[83]*TMath::Power(Q,3) + 
      2*RMvcs[3][2]*Ur[16]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*RMvcs[0][0]*Ur[10]*TMath::Power(Q,2) - 
      2*RMvcs[0][1]*Ur[29]*TMath::Power(Q,3) - 
      2*RMvcs[0][2]*Ur[9]*TMath::Power(Q,4)) + 
   Jem[2][2]*(-2*RMvcs[2][0]*Ur[10]*TMath::Power(Q,2) - 
      2*RMvcs[2][1]*Ur[29]*TMath::Power(Q,3) - 
      2*RMvcs[2][2]*Ur[9]*TMath::Power(Q,4));
        SigmaIPol0[2]=Jem[0][1]*(-2*RMvcs[0][0]*Ur[3]*TMath::Power(Q,3) - 
      2*RMvcs[0][2]*Ur[40]*TMath::Power(Q,3) + 
      2*RMvcs[0][1]*Ur[26]*TMath::Power(Q,4)) + 
   Jem[2][1]*(-2*RMvcs[2][0]*Ur[3]*TMath::Power(Q,3) - 
      2*RMvcs[2][2]*Ur[40]*TMath::Power(Q,3) + 
      2*RMvcs[2][1]*Ur[26]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*RMvcs[0][0]*Ur[11]*TMath::Power(Q,3) + 
      2*RMvcs[0][1]*Ur[31]*TMath::Power(Q,4) - 
      2*RMvcs[0][2]*Ur[46]*TMath::Power(Q,5)) + 
   Jem[2][2]*(-2*RMvcs[2][0]*Ur[11]*TMath::Power(Q,3) + 
      2*RMvcs[2][1]*Ur[31]*TMath::Power(Q,4) - 
      2*RMvcs[2][2]*Ur[46]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*RMvcs[1][0]*Ur[70]*TMath::Power(Q,3) + 
      2*RMvcs[1][1]*Ur[85]*TMath::Power(Q,4) - 
      2*RMvcs[1][2]*Ur[97]*TMath::Power(Q,5)) + 
   Jem[3][2]*(-2*RMvcs[3][0]*Ur[70]*TMath::Power(Q,3) + 
      2*RMvcs[3][1]*Ur[85]*TMath::Power(Q,4) - 
      2*RMvcs[3][2]*Ur[97]*TMath::Power(Q,5));
        SigmaIPol0[3]=Jem[0][1]*(-2*RMvcs[0][2]*Ur[42]*TMath::Power(Q,3) - 
      2*RMvcs[0][0]*Ur[5]*TMath::Power(Q,3) + 
      2*RMvcs[0][1]*Ur[27]*TMath::Power(Q,4)) + 
   Jem[2][1]*(-2*RMvcs[2][2]*Ur[42]*TMath::Power(Q,3) - 
      2*RMvcs[2][0]*Ur[5]*TMath::Power(Q,3) + 
      2*RMvcs[2][1]*Ur[27]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*RMvcs[0][0]*Ur[13]*TMath::Power(Q,3) + 
      2*RMvcs[0][1]*Ur[32]*TMath::Power(Q,4) - 
      2*RMvcs[0][2]*Ur[48]*TMath::Power(Q,5)) + 
   Jem[2][2]*(-2*RMvcs[2][0]*Ur[13]*TMath::Power(Q,3) + 
      2*RMvcs[2][1]*Ur[32]*TMath::Power(Q,4) - 
      2*RMvcs[2][2]*Ur[48]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*RMvcs[1][0]*Ur[71]*TMath::Power(Q,3) - 
      2*RMvcs[1][1]*Ur[85]*TMath::Power(Q,4) - 
      2*RMvcs[1][2]*Ur[98]*TMath::Power(Q,5)) + 
   Jem[3][2]*(-2*RMvcs[3][0]*Ur[71]*TMath::Power(Q,3) - 
      2*RMvcs[3][1]*Ur[85]*TMath::Power(Q,4) - 
      2*RMvcs[3][2]*Ur[98]*TMath::Power(Q,5));
        SigmaIPol0[4]=Jem[0][1]*(-2*RMvcs[0][2]*Ur[2]*TMath::Power(Q,2) + 
      2*RMvcs[0][1]*Ur[25]*TMath::Power(Q,3) - 
      2*RMvcs[0][0]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[2][1]*(-2*RMvcs[2][2]*Ur[2]*TMath::Power(Q,2) + 
      2*RMvcs[2][1]*Ur[25]*TMath::Power(Q,3) - 
      2*RMvcs[2][0]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*RMvcs[0][2]*Ur[10]*TMath::Power(Q,2) - 
      2*RMvcs[0][0]*Ur[9]*TMath::Power(Q,4) + 
      2*RMvcs[0][1]*Ur[30]*TMath::Power(Q,5)) + 
   Jem[2][2]*(-2*RMvcs[2][2]*Ur[10]*TMath::Power(Q,2) - 
      2*RMvcs[2][0]*Ur[9]*TMath::Power(Q,4) + 
      2*RMvcs[2][1]*Ur[30]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*RMvcs[1][2]*Ur[51]*TMath::Power(Q,2) + 
      2*RMvcs[1][0]*Ur[16]*TMath::Power(Q,4) + 
      2*RMvcs[1][1]*Ur[84]*TMath::Power(Q,5)) + 
   Jem[3][2]*(-2*RMvcs[3][2]*Ur[51]*TMath::Power(Q,2) + 
      2*RMvcs[3][0]*Ur[16]*TMath::Power(Q,4) + 
      2*RMvcs[3][1]*Ur[84]*TMath::Power(Q,5));
        SigmaIPol0[5]=Jem[0][1]*(-2*RMvcs[0][2]*Ur[41]*Q - 
      2*RMvcs[0][0]*Ur[4]*TMath::Power(Q,3)) + 
   Jem[2][1]*(-2*RMvcs[2][2]*Ur[41]*Q - 
      2*RMvcs[2][0]*Ur[4]*TMath::Power(Q,3)) + 
   Jem[0][2]*(-2*RMvcs[0][2]*Ur[47]*Q - 
      2*RMvcs[0][0]*Ur[12]*TMath::Power(Q,5)) + 
   Jem[2][2]*(-2*RMvcs[2][2]*Ur[47]*Q - 
      2*RMvcs[2][0]*Ur[12]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*RMvcs[1][2]*Ur[53]*Q + 
      2*RMvcs[1][0]*Ur[18]*TMath::Power(Q,5)) + 
   Jem[3][2]*(-2*RMvcs[3][2]*Ur[53]*Q + 
      2*RMvcs[3][0]*Ur[18]*TMath::Power(Q,5));
        SigmaIPol0[6]=Jem[0][1]*(-2*IMvcs[0][1]*Ur[28]*TMath::Power(Q,2) + 
      2*IMvcs[0][2]*Ur[44]*TMath::Power(Q,3) + 
      2*IMvcs[0][0]*Ur[7]*TMath::Power(Q,3)) + 
   Jem[2][1]*(-2*IMvcs[2][1]*Ur[28]*TMath::Power(Q,2) + 
      2*IMvcs[2][2]*Ur[44]*TMath::Power(Q,3) + 
      2*IMvcs[2][0]*Ur[7]*TMath::Power(Q,3)) + 
   Jem[0][2]*(2*IMvcs[0][0]*Ur[15]*TMath::Power(Q,3) - 
      2*IMvcs[0][1]*Ur[33]*TMath::Power(Q,4) + 
      2*IMvcs[0][2]*Ur[50]*TMath::Power(Q,5)) + 
   Jem[2][2]*(2*IMvcs[2][0]*Ur[15]*TMath::Power(Q,3) - 
      2*IMvcs[2][1]*Ur[33]*TMath::Power(Q,4) + 
      2*IMvcs[2][2]*Ur[50]*TMath::Power(Q,5)) + 
   Jem[1][2]*(2*IMvcs[1][0]*Ur[72]*TMath::Power(Q,3) + 
      2*IMvcs[1][1]*Ur[85]*TMath::Power(Q,4) + 
      2*IMvcs[1][2]*Ur[99]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*IMvcs[3][0]*Ur[72]*TMath::Power(Q,3) + 
      2*IMvcs[3][1]*Ur[85]*TMath::Power(Q,4) + 
      2*IMvcs[3][2]*Ur[99]*TMath::Power(Q,5));
        SigmaIPol0[7]=Jem[0][1]*(2*IMvcs[0][2]*Ur[43]*TMath::Power(Q,2) - 
      4*IMvcs[0][1]*Ur[25]*TMath::Power(Q,3) + 
      2*IMvcs[0][0]*Ur[6]*TMath::Power(Q,4)) + 
   Jem[2][1]*(2*IMvcs[2][2]*Ur[43]*TMath::Power(Q,2) - 
      4*IMvcs[2][1]*Ur[25]*TMath::Power(Q,3) + 
      2*IMvcs[2][0]*Ur[6]*TMath::Power(Q,4)) + 
   Jem[0][2]*(2*IMvcs[0][2]*Ur[49]*TMath::Power(Q,2) + 
      2*IMvcs[0][0]*Ur[14]*TMath::Power(Q,4) - 
      4*IMvcs[0][1]*Ur[30]*TMath::Power(Q,5)) + 
   Jem[2][2]*(2*IMvcs[2][2]*Ur[49]*TMath::Power(Q,2) + 
      2*IMvcs[2][0]*Ur[14]*TMath::Power(Q,4) - 
      4*IMvcs[2][1]*Ur[30]*TMath::Power(Q,5)) + 
   Jem[1][2]*(2*IMvcs[1][2]*Ur[56]*TMath::Power(Q,2) - 
      2*IMvcs[1][0]*Ur[21]*TMath::Power(Q,4) - 
      4*IMvcs[1][1]*Ur[84]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*IMvcs[3][2]*Ur[56]*TMath::Power(Q,2) - 
      2*IMvcs[3][0]*Ur[21]*TMath::Power(Q,4) - 
      4*IMvcs[3][1]*Ur[84]*TMath::Power(Q,5));
        SigmaIPolX[0]=Jem[0][1]*(2*IMvcs[1][0]*Ur[59]*TMath::Power(Q,3) + 
      2*IMvcs[1][2]*Ur[86]*TMath::Power(Q,3) - 
      2*IMvcs[1][1]*Ur[74]*TMath::Power(Q,4)) + 
   Jem[2][1]*(-2*IMvcs[3][0]*Ur[59]*TMath::Power(Q,3) - 
      2*IMvcs[3][2]*Ur[86]*TMath::Power(Q,3) + 
      2*IMvcs[3][1]*Ur[74]*TMath::Power(Q,4)) + 
   Jem[1][2]*(2*IMvcs[0][0]*Ur[17]*TMath::Power(Q,3) - 
      2*IMvcs[0][1]*Ur[35]*TMath::Power(Q,4) + 
      2*IMvcs[0][2]*Ur[52]*TMath::Power(Q,5)) + 
   Jem[3][2]*(-2*IMvcs[2][0]*Ur[17]*TMath::Power(Q,3) + 
      2*IMvcs[2][1]*Ur[35]*TMath::Power(Q,4) - 
      2*IMvcs[2][2]*Ur[52]*TMath::Power(Q,5)) + 
   Jem[0][2]*(2*IMvcs[1][0]*Ur[64]*TMath::Power(Q,3) - 
      2*IMvcs[1][1]*Ur[79]*TMath::Power(Q,4) + 
      2*IMvcs[1][2]*Ur[91]*TMath::Power(Q,5)) + 
   Jem[2][2]*(-2*IMvcs[3][0]*Ur[64]*TMath::Power(Q,3) + 
      2*IMvcs[3][1]*Ur[79]*TMath::Power(Q,4) - 
      2*IMvcs[3][2]*Ur[91]*TMath::Power(Q,5));
        SigmaIPolX[1]=Jem[0][1]*(2*IMvcs[1][0]*Ur[60]*TMath::Power(Q,3) + 
      2*IMvcs[1][2]*Ur[87]*TMath::Power(Q,3) - 
      2*IMvcs[1][1]*Ur[75]*TMath::Power(Q,4)) + 
   Jem[2][1]*(-2*IMvcs[3][0]*Ur[60]*TMath::Power(Q,3) - 
      2*IMvcs[3][2]*Ur[87]*TMath::Power(Q,3) + 
      2*IMvcs[3][1]*Ur[75]*TMath::Power(Q,4)) + 
   Jem[1][2]*(2*IMvcs[0][0]*Ur[19]*TMath::Power(Q,3) - 
      2*IMvcs[0][1]*Ur[36]*TMath::Power(Q,4) + 
      2*IMvcs[0][2]*Ur[54]*TMath::Power(Q,5)) + 
   Jem[3][2]*(-2*IMvcs[2][0]*Ur[19]*TMath::Power(Q,3) + 
      2*IMvcs[2][1]*Ur[36]*TMath::Power(Q,4) - 
      2*IMvcs[2][2]*Ur[54]*TMath::Power(Q,5)) + 
   Jem[0][2]*(2*IMvcs[1][0]*Ur[65]*TMath::Power(Q,3) - 
      2*IMvcs[1][1]*Ur[80]*TMath::Power(Q,4) + 
      2*IMvcs[1][2]*Ur[92]*TMath::Power(Q,5)) + 
   Jem[2][2]*(-2*IMvcs[3][0]*Ur[65]*TMath::Power(Q,3) + 
      2*IMvcs[3][1]*Ur[80]*TMath::Power(Q,4) - 
      2*IMvcs[3][2]*Ur[92]*TMath::Power(Q,5));
        SigmaIPolX[2]=Jem[0][1]*(2*IMvcs[1][2]*Ur[2]*TMath::Power(Q,2) - 
      2*IMvcs[1][1]*Ur[73]*TMath::Power(Q,3) - 
      2*IMvcs[1][0]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[2][1]*(-2*IMvcs[3][2]*Ur[2]*TMath::Power(Q,2) + 
      2*IMvcs[3][1]*Ur[73]*TMath::Power(Q,3) + 
      2*IMvcs[3][0]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[1][2]*(2*IMvcs[0][2]*Ur[51]*TMath::Power(Q,2) + 
      2*IMvcs[0][0]*Ur[16]*TMath::Power(Q,4) - 
      2*IMvcs[0][1]*Ur[34]*TMath::Power(Q,5)) + 
   Jem[3][2]*(-2*IMvcs[2][2]*Ur[51]*TMath::Power(Q,2) - 
      2*IMvcs[2][0]*Ur[16]*TMath::Power(Q,4) + 
      2*IMvcs[2][1]*Ur[34]*TMath::Power(Q,5)) + 
   Jem[0][2]*(2*IMvcs[1][2]*Ur[10]*TMath::Power(Q,2) - 
      2*IMvcs[1][0]*Ur[9]*TMath::Power(Q,4) - 
      2*IMvcs[1][1]*Ur[78]*TMath::Power(Q,5)) + 
   Jem[2][2]*(-2*IMvcs[3][2]*Ur[10]*TMath::Power(Q,2) + 
      2*IMvcs[3][0]*Ur[9]*TMath::Power(Q,4) + 
      2*IMvcs[3][1]*Ur[78]*TMath::Power(Q,5));
        SigmaIPolX[3]=Jem[0][1]*(2*IMvcs[1][2]*Ur[41]*Q - 
      2*IMvcs[1][0]*Ur[4]*TMath::Power(Q,3)) + 
   Jem[2][1]*(-2*IMvcs[3][2]*Ur[41]*Q + 
      2*IMvcs[3][0]*Ur[4]*TMath::Power(Q,3)) + 
   Jem[0][2]*(2*IMvcs[1][2]*Ur[47]*Q - 
      2*IMvcs[1][0]*Ur[12]*TMath::Power(Q,5)) + 
   Jem[2][2]*(-2*IMvcs[3][2]*Ur[47]*Q + 
      2*IMvcs[3][0]*Ur[12]*TMath::Power(Q,5)) + 
   Jem[1][2]*(2*IMvcs[0][2]*Ur[53]*Q + 
      2*IMvcs[0][0]*Ur[18]*TMath::Power(Q,5)) + 
   Jem[3][2]*(-2*IMvcs[2][2]*Ur[53]*Q - 
      2*IMvcs[2][0]*Ur[18]*TMath::Power(Q,5));
        SigmaIPolX[4]=Jem[1][2]*(-2*RMvcs[0][0]*Ur[20]*TMath::Power(Q,2) + 
      2*RMvcs[0][1]*Ur[37]*TMath::Power(Q,3) - 
      2*RMvcs[0][2]*Ur[55]*TMath::Power(Q,4)) + 
   Jem[3][2]*(2*RMvcs[2][0]*Ur[20]*TMath::Power(Q,2) - 
      2*RMvcs[2][1]*Ur[37]*TMath::Power(Q,3) + 
      2*RMvcs[2][2]*Ur[55]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*RMvcs[1][0]*Ur[61]*TMath::Power(Q,2) + 
      2*RMvcs[1][1]*Ur[76]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[88]*TMath::Power(Q,4)) + 
   Jem[2][1]*(2*RMvcs[3][0]*Ur[61]*TMath::Power(Q,2) - 
      2*RMvcs[3][1]*Ur[76]*TMath::Power(Q,3) + 
      2*RMvcs[3][2]*Ur[88]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*RMvcs[1][0]*Ur[66]*TMath::Power(Q,2) + 
      2*RMvcs[1][1]*Ur[81]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[93]*TMath::Power(Q,4)) + 
   Jem[2][2]*(2*RMvcs[3][0]*Ur[66]*TMath::Power(Q,2) - 
      2*RMvcs[3][1]*Ur[81]*TMath::Power(Q,3) + 
      2*RMvcs[3][2]*Ur[93]*TMath::Power(Q,4));
        SigmaIPolX[5]=Jem[1][2]*(-2*RMvcs[0][0]*Ur[22]*TMath::Power(Q,2) - 
      2*RMvcs[0][1]*Ur[37]*TMath::Power(Q,3) - 
      2*RMvcs[0][2]*Ur[57]*TMath::Power(Q,4)) + 
   Jem[3][2]*(2*RMvcs[2][0]*Ur[22]*TMath::Power(Q,2) + 
      2*RMvcs[2][1]*Ur[37]*TMath::Power(Q,3) + 
      2*RMvcs[2][2]*Ur[57]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*RMvcs[1][0]*Ur[62]*TMath::Power(Q,2) - 
      2*RMvcs[1][1]*Ur[76]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[89]*TMath::Power(Q,4)) + 
   Jem[2][1]*(2*RMvcs[3][0]*Ur[62]*TMath::Power(Q,2) + 
      2*RMvcs[3][1]*Ur[76]*TMath::Power(Q,3) + 
      2*RMvcs[3][2]*Ur[89]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*RMvcs[1][0]*Ur[67]*TMath::Power(Q,2) - 
      2*RMvcs[1][1]*Ur[81]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[94]*TMath::Power(Q,4)) + 
   Jem[2][2]*(2*RMvcs[3][0]*Ur[67]*TMath::Power(Q,2) + 
      2*RMvcs[3][1]*Ur[81]*TMath::Power(Q,3) + 
      2*RMvcs[3][2]*Ur[94]*TMath::Power(Q,4));
        SigmaIPolX[6]=Jem[0][1]*(2*RMvcs[1][1]*Ur[77]*TMath::Power(Q,2) - 
      2*RMvcs[1][0]*Ur[63]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[90]*TMath::Power(Q,3)) + 
   Jem[2][1]*(-2*RMvcs[3][1]*Ur[77]*TMath::Power(Q,2) + 
      2*RMvcs[3][0]*Ur[63]*TMath::Power(Q,3) + 
      2*RMvcs[3][2]*Ur[90]*TMath::Power(Q,3)) + 
   Jem[1][2]*(-2*RMvcs[0][0]*Ur[23]*TMath::Power(Q,3) + 
      2*RMvcs[0][1]*Ur[38]*TMath::Power(Q,4) - 
      2*RMvcs[0][2]*Ur[58]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*RMvcs[2][0]*Ur[23]*TMath::Power(Q,3) - 
      2*RMvcs[2][1]*Ur[38]*TMath::Power(Q,4) + 
      2*RMvcs[2][2]*Ur[58]*TMath::Power(Q,5)) + 
   Jem[0][2]*(-2*RMvcs[1][0]*Ur[68]*TMath::Power(Q,3) + 
      2*RMvcs[1][1]*Ur[82]*TMath::Power(Q,4) - 
      2*RMvcs[1][2]*Ur[95]*TMath::Power(Q,5)) + 
   Jem[2][2]*(2*RMvcs[3][0]*Ur[68]*TMath::Power(Q,3) - 
      2*RMvcs[3][1]*Ur[82]*TMath::Power(Q,4) + 
      2*RMvcs[3][2]*Ur[95]*TMath::Power(Q,5));
        SigmaIPolX[7]=Jem[0][1]*(-2*RMvcs[1][2]*Ur[43]*TMath::Power(Q,2) + 
      4*RMvcs[1][1]*Ur[73]*TMath::Power(Q,3) + 
      2*RMvcs[1][0]*Ur[6]*TMath::Power(Q,4)) + 
   Jem[2][1]*(2*RMvcs[3][2]*Ur[43]*TMath::Power(Q,2) - 
      4*RMvcs[3][1]*Ur[73]*TMath::Power(Q,3) - 
      2*RMvcs[3][0]*Ur[6]*TMath::Power(Q,4)) + 
   Jem[1][2]*(-2*RMvcs[0][2]*Ur[56]*TMath::Power(Q,2) - 
      2*RMvcs[0][0]*Ur[21]*TMath::Power(Q,4) + 
      4*RMvcs[0][1]*Ur[34]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*RMvcs[2][2]*Ur[56]*TMath::Power(Q,2) + 
      2*RMvcs[2][0]*Ur[21]*TMath::Power(Q,4) - 
      4*RMvcs[2][1]*Ur[34]*TMath::Power(Q,5)) + 
   Jem[0][2]*(-2*RMvcs[1][2]*Ur[49]*TMath::Power(Q,2) + 
      2*RMvcs[1][0]*Ur[14]*TMath::Power(Q,4) + 
      4*RMvcs[1][1]*Ur[78]*TMath::Power(Q,5)) + 
   Jem[2][2]*(2*RMvcs[3][2]*Ur[49]*TMath::Power(Q,2) - 
      2*RMvcs[3][0]*Ur[14]*TMath::Power(Q,4) - 
      4*RMvcs[3][1]*Ur[78]*TMath::Power(Q,5));
        SigmaIPolY[0]=Jem[2][1]*(2*IMvcs[0][0]*Ur[0]*TMath::Power(Q,2) - 
      2*IMvcs[0][1]*Ur[24]*TMath::Power(Q,3) + 
      2*IMvcs[0][2]*Ur[39]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*IMvcs[2][0]*Ur[0]*TMath::Power(Q,2) + 
      2*IMvcs[2][1]*Ur[24]*TMath::Power(Q,3) - 
      2*IMvcs[2][2]*Ur[39]*TMath::Power(Q,4)) + 
   Jem[2][2]*(2*IMvcs[0][0]*Ur[8]*TMath::Power(Q,2) - 
      2*IMvcs[0][1]*Ur[29]*TMath::Power(Q,3) + 
      2*IMvcs[0][2]*Ur[45]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*IMvcs[2][0]*Ur[8]*TMath::Power(Q,2) + 
      2*IMvcs[2][1]*Ur[29]*TMath::Power(Q,3) - 
      2*IMvcs[2][2]*Ur[45]*TMath::Power(Q,4)) + 
   Jem[3][2]*(2*IMvcs[1][0]*Ur[69]*TMath::Power(Q,2) - 
      2*IMvcs[1][1]*Ur[83]*TMath::Power(Q,3) + 
      2*IMvcs[1][2]*Ur[96]*TMath::Power(Q,4)) + 
   Jem[1][2]*(-2*IMvcs[3][0]*Ur[69]*TMath::Power(Q,2) + 
      2*IMvcs[3][1]*Ur[83]*TMath::Power(Q,3) - 
      2*IMvcs[3][2]*Ur[96]*TMath::Power(Q,4));
        SigmaIPolY[1]=Jem[2][1]*(2*IMvcs[0][0]*Ur[2]*TMath::Power(Q,2) + 
      2*IMvcs[0][1]*Ur[24]*TMath::Power(Q,3) + 
      2*IMvcs[0][2]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*IMvcs[2][0]*Ur[2]*TMath::Power(Q,2) - 
      2*IMvcs[2][1]*Ur[24]*TMath::Power(Q,3) - 
      2*IMvcs[2][2]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[3][2]*(2*IMvcs[1][0]*Ur[51]*TMath::Power(Q,2) + 
      2*IMvcs[1][1]*Ur[83]*TMath::Power(Q,3) - 
      2*IMvcs[1][2]*Ur[16]*TMath::Power(Q,4)) + 
   Jem[1][2]*(-2*IMvcs[3][0]*Ur[51]*TMath::Power(Q,2) - 
      2*IMvcs[3][1]*Ur[83]*TMath::Power(Q,3) + 
      2*IMvcs[3][2]*Ur[16]*TMath::Power(Q,4)) + 
   Jem[2][2]*(2*IMvcs[0][0]*Ur[10]*TMath::Power(Q,2) + 
      2*IMvcs[0][1]*Ur[29]*TMath::Power(Q,3) + 
      2*IMvcs[0][2]*Ur[9]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*IMvcs[2][0]*Ur[10]*TMath::Power(Q,2) - 
      2*IMvcs[2][1]*Ur[29]*TMath::Power(Q,3) - 
      2*IMvcs[2][2]*Ur[9]*TMath::Power(Q,4));
        SigmaIPolY[2]=Jem[2][1]*(2*IMvcs[0][0]*Ur[3]*TMath::Power(Q,3) + 
      2*IMvcs[0][2]*Ur[40]*TMath::Power(Q,3) - 
      2*IMvcs[0][1]*Ur[26]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*IMvcs[2][0]*Ur[3]*TMath::Power(Q,3) - 
      2*IMvcs[2][2]*Ur[40]*TMath::Power(Q,3) + 
      2*IMvcs[2][1]*Ur[26]*TMath::Power(Q,4)) + 
   Jem[2][2]*(2*IMvcs[0][0]*Ur[11]*TMath::Power(Q,3) - 
      2*IMvcs[0][1]*Ur[31]*TMath::Power(Q,4) + 
      2*IMvcs[0][2]*Ur[46]*TMath::Power(Q,5)) + 
   Jem[0][2]*(-2*IMvcs[2][0]*Ur[11]*TMath::Power(Q,3) + 
      2*IMvcs[2][1]*Ur[31]*TMath::Power(Q,4) - 
      2*IMvcs[2][2]*Ur[46]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*IMvcs[1][0]*Ur[70]*TMath::Power(Q,3) - 
      2*IMvcs[1][1]*Ur[85]*TMath::Power(Q,4) + 
      2*IMvcs[1][2]*Ur[97]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*IMvcs[3][0]*Ur[70]*TMath::Power(Q,3) + 
      2*IMvcs[3][1]*Ur[85]*TMath::Power(Q,4) - 
      2*IMvcs[3][2]*Ur[97]*TMath::Power(Q,5));
        SigmaIPolY[3]=Jem[2][1]*(2*IMvcs[0][2]*Ur[42]*TMath::Power(Q,3) + 
      2*IMvcs[0][0]*Ur[5]*TMath::Power(Q,3) - 
      2*IMvcs[0][1]*Ur[27]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*IMvcs[2][2]*Ur[42]*TMath::Power(Q,3) - 
      2*IMvcs[2][0]*Ur[5]*TMath::Power(Q,3) + 
      2*IMvcs[2][1]*Ur[27]*TMath::Power(Q,4)) + 
   Jem[2][2]*(2*IMvcs[0][0]*Ur[13]*TMath::Power(Q,3) - 
      2*IMvcs[0][1]*Ur[32]*TMath::Power(Q,4) + 
      2*IMvcs[0][2]*Ur[48]*TMath::Power(Q,5)) + 
   Jem[0][2]*(-2*IMvcs[2][0]*Ur[13]*TMath::Power(Q,3) + 
      2*IMvcs[2][1]*Ur[32]*TMath::Power(Q,4) - 
      2*IMvcs[2][2]*Ur[48]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*IMvcs[1][0]*Ur[71]*TMath::Power(Q,3) + 
      2*IMvcs[1][1]*Ur[85]*TMath::Power(Q,4) + 
      2*IMvcs[1][2]*Ur[98]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*IMvcs[3][0]*Ur[71]*TMath::Power(Q,3) - 
      2*IMvcs[3][1]*Ur[85]*TMath::Power(Q,4) - 
      2*IMvcs[3][2]*Ur[98]*TMath::Power(Q,5));
        SigmaIPolY[4]=Jem[2][1]*(2*IMvcs[0][2]*Ur[2]*TMath::Power(Q,2) - 
      2*IMvcs[0][1]*Ur[25]*TMath::Power(Q,3) + 
      2*IMvcs[0][0]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*IMvcs[2][2]*Ur[2]*TMath::Power(Q,2) + 
      2*IMvcs[2][1]*Ur[25]*TMath::Power(Q,3) - 
      2*IMvcs[2][0]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[2][2]*(2*IMvcs[0][2]*Ur[10]*TMath::Power(Q,2) + 
      2*IMvcs[0][0]*Ur[9]*TMath::Power(Q,4) - 
      2*IMvcs[0][1]*Ur[30]*TMath::Power(Q,5)) + 
   Jem[0][2]*(-2*IMvcs[2][2]*Ur[10]*TMath::Power(Q,2) - 
      2*IMvcs[2][0]*Ur[9]*TMath::Power(Q,4) + 
      2*IMvcs[2][1]*Ur[30]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*IMvcs[1][2]*Ur[51]*TMath::Power(Q,2) - 
      2*IMvcs[1][0]*Ur[16]*TMath::Power(Q,4) - 
      2*IMvcs[1][1]*Ur[84]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*IMvcs[3][2]*Ur[51]*TMath::Power(Q,2) + 
      2*IMvcs[3][0]*Ur[16]*TMath::Power(Q,4) + 
      2*IMvcs[3][1]*Ur[84]*TMath::Power(Q,5));
        SigmaIPolY[5]=Jem[2][1]*(2*IMvcs[0][2]*Ur[41]*Q + 
      2*IMvcs[0][0]*Ur[4]*TMath::Power(Q,3)) + 
   Jem[0][1]*(-2*IMvcs[2][2]*Ur[41]*Q - 
      2*IMvcs[2][0]*Ur[4]*TMath::Power(Q,3)) + 
   Jem[2][2]*(2*IMvcs[0][2]*Ur[47]*Q + 
      2*IMvcs[0][0]*Ur[12]*TMath::Power(Q,5)) + 
   Jem[0][2]*(-2*IMvcs[2][2]*Ur[47]*Q - 
      2*IMvcs[2][0]*Ur[12]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*IMvcs[1][2]*Ur[53]*Q - 
      2*IMvcs[1][0]*Ur[18]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*IMvcs[3][2]*Ur[53]*Q + 
      2*IMvcs[3][0]*Ur[18]*TMath::Power(Q,5));
        SigmaIPolY[6]=Jem[2][1]*(-2*RMvcs[0][1]*Ur[28]*TMath::Power(Q,2) + 
      2*RMvcs[0][2]*Ur[44]*TMath::Power(Q,3) + 
      2*RMvcs[0][0]*Ur[7]*TMath::Power(Q,3)) + 
   Jem[0][1]*(2*RMvcs[2][1]*Ur[28]*TMath::Power(Q,2) - 
      2*RMvcs[2][2]*Ur[44]*TMath::Power(Q,3) - 
      2*RMvcs[2][0]*Ur[7]*TMath::Power(Q,3)) + 
   Jem[2][2]*(2*RMvcs[0][0]*Ur[15]*TMath::Power(Q,3) - 
      2*RMvcs[0][1]*Ur[33]*TMath::Power(Q,4) + 
      2*RMvcs[0][2]*Ur[50]*TMath::Power(Q,5)) + 
   Jem[0][2]*(-2*RMvcs[2][0]*Ur[15]*TMath::Power(Q,3) + 
      2*RMvcs[2][1]*Ur[33]*TMath::Power(Q,4) - 
      2*RMvcs[2][2]*Ur[50]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*RMvcs[1][0]*Ur[72]*TMath::Power(Q,3) + 
      2*RMvcs[1][1]*Ur[85]*TMath::Power(Q,4) + 
      2*RMvcs[1][2]*Ur[99]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*RMvcs[3][0]*Ur[72]*TMath::Power(Q,3) - 
      2*RMvcs[3][1]*Ur[85]*TMath::Power(Q,4) - 
      2*RMvcs[3][2]*Ur[99]*TMath::Power(Q,5));
        SigmaIPolY[7]=Jem[2][1]*(2*RMvcs[0][2]*Ur[43]*TMath::Power(Q,2) - 
      4*RMvcs[0][1]*Ur[25]*TMath::Power(Q,3) + 
      2*RMvcs[0][0]*Ur[6]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*RMvcs[2][2]*Ur[43]*TMath::Power(Q,2) + 
      4*RMvcs[2][1]*Ur[25]*TMath::Power(Q,3) - 
      2*RMvcs[2][0]*Ur[6]*TMath::Power(Q,4)) + 
   Jem[2][2]*(2*RMvcs[0][2]*Ur[49]*TMath::Power(Q,2) + 
      2*RMvcs[0][0]*Ur[14]*TMath::Power(Q,4) - 
      4*RMvcs[0][1]*Ur[30]*TMath::Power(Q,5)) + 
   Jem[0][2]*(-2*RMvcs[2][2]*Ur[49]*TMath::Power(Q,2) - 
      2*RMvcs[2][0]*Ur[14]*TMath::Power(Q,4) + 
      4*RMvcs[2][1]*Ur[30]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*RMvcs[1][2]*Ur[56]*TMath::Power(Q,2) - 
      2*RMvcs[1][0]*Ur[21]*TMath::Power(Q,4) - 
      4*RMvcs[1][1]*Ur[84]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*RMvcs[3][2]*Ur[56]*TMath::Power(Q,2) + 
      2*RMvcs[3][0]*Ur[21]*TMath::Power(Q,4) + 
      4*RMvcs[3][1]*Ur[84]*TMath::Power(Q,5));
        SigmaIPolZ[0]=Jem[2][1]*(2*IMvcs[1][0]*Ur[59]*TMath::Power(Q,3) + 
      2*IMvcs[1][2]*Ur[86]*TMath::Power(Q,3) - 
      2*IMvcs[1][1]*Ur[74]*TMath::Power(Q,4)) + 
   Jem[0][1]*(2*IMvcs[3][0]*Ur[59]*TMath::Power(Q,3) + 
      2*IMvcs[3][2]*Ur[86]*TMath::Power(Q,3) - 
      2*IMvcs[3][1]*Ur[74]*TMath::Power(Q,4)) + 
   Jem[3][2]*(2*IMvcs[0][0]*Ur[17]*TMath::Power(Q,3) - 
      2*IMvcs[0][1]*Ur[35]*TMath::Power(Q,4) + 
      2*IMvcs[0][2]*Ur[52]*TMath::Power(Q,5)) + 
   Jem[1][2]*(2*IMvcs[2][0]*Ur[17]*TMath::Power(Q,3) - 
      2*IMvcs[2][1]*Ur[35]*TMath::Power(Q,4) + 
      2*IMvcs[2][2]*Ur[52]*TMath::Power(Q,5)) + 
   Jem[2][2]*(2*IMvcs[1][0]*Ur[64]*TMath::Power(Q,3) - 
      2*IMvcs[1][1]*Ur[79]*TMath::Power(Q,4) + 
      2*IMvcs[1][2]*Ur[91]*TMath::Power(Q,5)) + 
   Jem[0][2]*(2*IMvcs[3][0]*Ur[64]*TMath::Power(Q,3) - 
      2*IMvcs[3][1]*Ur[79]*TMath::Power(Q,4) + 
      2*IMvcs[3][2]*Ur[91]*TMath::Power(Q,5));
        SigmaIPolZ[1]=Jem[2][1]*(2*IMvcs[1][0]*Ur[60]*TMath::Power(Q,3) + 
      2*IMvcs[1][2]*Ur[87]*TMath::Power(Q,3) - 
      2*IMvcs[1][1]*Ur[75]*TMath::Power(Q,4)) + 
   Jem[0][1]*(2*IMvcs[3][0]*Ur[60]*TMath::Power(Q,3) + 
      2*IMvcs[3][2]*Ur[87]*TMath::Power(Q,3) - 
      2*IMvcs[3][1]*Ur[75]*TMath::Power(Q,4)) + 
   Jem[3][2]*(2*IMvcs[0][0]*Ur[19]*TMath::Power(Q,3) - 
      2*IMvcs[0][1]*Ur[36]*TMath::Power(Q,4) + 
      2*IMvcs[0][2]*Ur[54]*TMath::Power(Q,5)) + 
   Jem[1][2]*(2*IMvcs[2][0]*Ur[19]*TMath::Power(Q,3) - 
      2*IMvcs[2][1]*Ur[36]*TMath::Power(Q,4) + 
      2*IMvcs[2][2]*Ur[54]*TMath::Power(Q,5)) + 
   Jem[2][2]*(2*IMvcs[1][0]*Ur[65]*TMath::Power(Q,3) - 
      2*IMvcs[1][1]*Ur[80]*TMath::Power(Q,4) + 
      2*IMvcs[1][2]*Ur[92]*TMath::Power(Q,5)) + 
   Jem[0][2]*(2*IMvcs[3][0]*Ur[65]*TMath::Power(Q,3) - 
      2*IMvcs[3][1]*Ur[80]*TMath::Power(Q,4) + 
      2*IMvcs[3][2]*Ur[92]*TMath::Power(Q,5));
        SigmaIPolZ[2]=Jem[2][1]*(2*IMvcs[1][2]*Ur[2]*TMath::Power(Q,2) - 
      2*IMvcs[1][1]*Ur[73]*TMath::Power(Q,3) - 
      2*IMvcs[1][0]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[0][1]*(2*IMvcs[3][2]*Ur[2]*TMath::Power(Q,2) - 
      2*IMvcs[3][1]*Ur[73]*TMath::Power(Q,3) - 
      2*IMvcs[3][0]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[3][2]*(2*IMvcs[0][2]*Ur[51]*TMath::Power(Q,2) + 
      2*IMvcs[0][0]*Ur[16]*TMath::Power(Q,4) - 
      2*IMvcs[0][1]*Ur[34]*TMath::Power(Q,5)) + 
   Jem[1][2]*(2*IMvcs[2][2]*Ur[51]*TMath::Power(Q,2) + 
      2*IMvcs[2][0]*Ur[16]*TMath::Power(Q,4) - 
      2*IMvcs[2][1]*Ur[34]*TMath::Power(Q,5)) + 
   Jem[2][2]*(2*IMvcs[1][2]*Ur[10]*TMath::Power(Q,2) - 
      2*IMvcs[1][0]*Ur[9]*TMath::Power(Q,4) - 
      2*IMvcs[1][1]*Ur[78]*TMath::Power(Q,5)) + 
   Jem[0][2]*(2*IMvcs[3][2]*Ur[10]*TMath::Power(Q,2) - 
      2*IMvcs[3][0]*Ur[9]*TMath::Power(Q,4) - 
      2*IMvcs[3][1]*Ur[78]*TMath::Power(Q,5));
        SigmaIPolZ[3]=Jem[2][1]*(2*IMvcs[1][2]*Ur[41]*Q - 
      2*IMvcs[1][0]*Ur[4]*TMath::Power(Q,3)) + 
   Jem[0][1]*(2*IMvcs[3][2]*Ur[41]*Q - 
      2*IMvcs[3][0]*Ur[4]*TMath::Power(Q,3)) + 
   Jem[2][2]*(2*IMvcs[1][2]*Ur[47]*Q - 
      2*IMvcs[1][0]*Ur[12]*TMath::Power(Q,5)) + 
   Jem[0][2]*(2*IMvcs[3][2]*Ur[47]*Q - 
      2*IMvcs[3][0]*Ur[12]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*IMvcs[0][2]*Ur[53]*Q + 
      2*IMvcs[0][0]*Ur[18]*TMath::Power(Q,5)) + 
   Jem[1][2]*(2*IMvcs[2][2]*Ur[53]*Q + 
      2*IMvcs[2][0]*Ur[18]*TMath::Power(Q,5));
        SigmaIPolZ[4]=Jem[3][2]*(-2*RMvcs[0][0]*Ur[20]*TMath::Power(Q,2) + 
      2*RMvcs[0][1]*Ur[37]*TMath::Power(Q,3) - 
      2*RMvcs[0][2]*Ur[55]*TMath::Power(Q,4)) + 
   Jem[1][2]*(-2*RMvcs[2][0]*Ur[20]*TMath::Power(Q,2) + 
      2*RMvcs[2][1]*Ur[37]*TMath::Power(Q,3) - 
      2*RMvcs[2][2]*Ur[55]*TMath::Power(Q,4)) + 
   Jem[2][1]*(-2*RMvcs[1][0]*Ur[61]*TMath::Power(Q,2) + 
      2*RMvcs[1][1]*Ur[76]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[88]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*RMvcs[3][0]*Ur[61]*TMath::Power(Q,2) + 
      2*RMvcs[3][1]*Ur[76]*TMath::Power(Q,3) - 
      2*RMvcs[3][2]*Ur[88]*TMath::Power(Q,4)) + 
   Jem[2][2]*(-2*RMvcs[1][0]*Ur[66]*TMath::Power(Q,2) + 
      2*RMvcs[1][1]*Ur[81]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[93]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*RMvcs[3][0]*Ur[66]*TMath::Power(Q,2) + 
      2*RMvcs[3][1]*Ur[81]*TMath::Power(Q,3) - 
      2*RMvcs[3][2]*Ur[93]*TMath::Power(Q,4));
        SigmaIPolZ[5]=Jem[3][2]*(-2*RMvcs[0][0]*Ur[22]*TMath::Power(Q,2) - 
      2*RMvcs[0][1]*Ur[37]*TMath::Power(Q,3) - 
      2*RMvcs[0][2]*Ur[57]*TMath::Power(Q,4)) + 
   Jem[1][2]*(-2*RMvcs[2][0]*Ur[22]*TMath::Power(Q,2) - 
      2*RMvcs[2][1]*Ur[37]*TMath::Power(Q,3) - 
      2*RMvcs[2][2]*Ur[57]*TMath::Power(Q,4)) + 
   Jem[2][1]*(-2*RMvcs[1][0]*Ur[62]*TMath::Power(Q,2) - 
      2*RMvcs[1][1]*Ur[76]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[89]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*RMvcs[3][0]*Ur[62]*TMath::Power(Q,2) - 
      2*RMvcs[3][1]*Ur[76]*TMath::Power(Q,3) - 
      2*RMvcs[3][2]*Ur[89]*TMath::Power(Q,4)) + 
   Jem[2][2]*(-2*RMvcs[1][0]*Ur[67]*TMath::Power(Q,2) - 
      2*RMvcs[1][1]*Ur[81]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[94]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*RMvcs[3][0]*Ur[67]*TMath::Power(Q,2) - 
      2*RMvcs[3][1]*Ur[81]*TMath::Power(Q,3) - 
      2*RMvcs[3][2]*Ur[94]*TMath::Power(Q,4));
        SigmaIPolZ[6]=Jem[2][1]*(2*RMvcs[1][1]*Ur[77]*TMath::Power(Q,2) - 
      2*RMvcs[1][0]*Ur[63]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[90]*TMath::Power(Q,3)) + 
   Jem[0][1]*(2*RMvcs[3][1]*Ur[77]*TMath::Power(Q,2) - 
      2*RMvcs[3][0]*Ur[63]*TMath::Power(Q,3) - 
      2*RMvcs[3][2]*Ur[90]*TMath::Power(Q,3)) + 
   Jem[3][2]*(-2*RMvcs[0][0]*Ur[23]*TMath::Power(Q,3) + 
      2*RMvcs[0][1]*Ur[38]*TMath::Power(Q,4) - 
      2*RMvcs[0][2]*Ur[58]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*RMvcs[2][0]*Ur[23]*TMath::Power(Q,3) + 
      2*RMvcs[2][1]*Ur[38]*TMath::Power(Q,4) - 
      2*RMvcs[2][2]*Ur[58]*TMath::Power(Q,5)) + 
   Jem[2][2]*(-2*RMvcs[1][0]*Ur[68]*TMath::Power(Q,3) + 
      2*RMvcs[1][1]*Ur[82]*TMath::Power(Q,4) - 
      2*RMvcs[1][2]*Ur[95]*TMath::Power(Q,5)) + 
   Jem[0][2]*(-2*RMvcs[3][0]*Ur[68]*TMath::Power(Q,3) + 
      2*RMvcs[3][1]*Ur[82]*TMath::Power(Q,4) - 
      2*RMvcs[3][2]*Ur[95]*TMath::Power(Q,5));
        SigmaIPolZ[7]=Jem[2][1]*(-2*RMvcs[1][2]*Ur[43]*TMath::Power(Q,2) + 
      4*RMvcs[1][1]*Ur[73]*TMath::Power(Q,3) + 
      2*RMvcs[1][0]*Ur[6]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*RMvcs[3][2]*Ur[43]*TMath::Power(Q,2) + 
      4*RMvcs[3][1]*Ur[73]*TMath::Power(Q,3) + 
      2*RMvcs[3][0]*Ur[6]*TMath::Power(Q,4)) + 
   Jem[3][2]*(-2*RMvcs[0][2]*Ur[56]*TMath::Power(Q,2) - 
      2*RMvcs[0][0]*Ur[21]*TMath::Power(Q,4) + 
      4*RMvcs[0][1]*Ur[34]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*RMvcs[2][2]*Ur[56]*TMath::Power(Q,2) - 
      2*RMvcs[2][0]*Ur[21]*TMath::Power(Q,4) + 
      4*RMvcs[2][1]*Ur[34]*TMath::Power(Q,5)) + 
   Jem[2][2]*(-2*RMvcs[1][2]*Ur[49]*TMath::Power(Q,2) + 
      2*RMvcs[1][0]*Ur[14]*TMath::Power(Q,4) + 
      4*RMvcs[1][1]*Ur[78]*TMath::Power(Q,5)) + 
   Jem[0][2]*(-2*RMvcs[3][2]*Ur[49]*TMath::Power(Q,2) + 
      2*RMvcs[3][0]*Ur[14]*TMath::Power(Q,4) + 
      4*RMvcs[3][1]*Ur[78]*TMath::Power(Q,5));
    
    
   	 	// Flag
   	 	
   	 	InitExactVCSAndInterfCrossSections = kTRUE;
   	 	
   	 	
   	 	// Validation mode only
    
   	 	if ( Validation == kTRUE )
   	 	{	
    		ofstream outfile;    		
    		
    		outfile.open("CheckUrExact.dat");
    		for(Int_t i=0;i<100;i++)
			{
				outfile << Ur[i] << endl;
			} // end for i		
			outfile.close();	
		
			
    		outfile.open("CheckSigmaVCSExact.dat");
			outfile << SigmaVCSPol0[0] << endl;
			outfile << SigmaVCSPolZ[0] << endl;
			outfile << SigmaVCSPol0[1] << endl;
			outfile << SigmaVCSPolX[0] << endl;
			outfile << SigmaVCSPolX[1] << endl;
			outfile << SigmaVCSPolZ[1] << endl;
			outfile << SigmaVCSPol0[2] << endl;
			outfile << SigmaVCSPol0[3] << endl;
			outfile << SigmaVCSPolY[0] << endl;
			outfile << SigmaVCSPolY[1] << endl;
			outfile << SigmaVCSPolY[2] << endl;
			outfile << SigmaVCSPolY[3] << endl;
			outfile << SigmaVCSPolY[4] << endl;
			outfile << SigmaVCSPol0[4] << endl;
			outfile << SigmaVCSPolX[2] << endl;
			outfile << SigmaVCSPolZ[2] << endl;
			outfile << SigmaVCSPolX[3] << endl;
			outfile << SigmaVCSPolZ[3] << endl;		
			outfile.close();	
		
			
    		outfile.open("CheckSigmaIExact.dat");
			for(Int_t i=0;i<32;i++)
			{
				if (i<8)
				{
					outfile << SigmaIPol0[i] << endl;
				} // end if i
				if ( 8 <= i && i < 16 )
				{
					outfile << SigmaIPolX[i-8] << endl;
				} // end if i
				if ( 16 <= i && i < 24 )
				{
					outfile << SigmaIPolY[i-16] << endl;
				} // end if i
				if ( 24 <= i && i < 32)
				{
					outfile << SigmaIPolZ[i-24] << endl;
				} // end if i
			} // end for i    	
		outfile.close();		
				
		} // end if Validation  
		
//	} // end if InitExactVCSAndInterfCrossSections
	 
} // end MakeExactVCSAndInterfCrossSections


	
/*------------------------ Function MakeLeadingBHCrossSections() -----------------------*
 | Computes all the stuff to evaluate the cross section assuming the hadronic helicity  |
 | amplitudes are given, i.e. initializes :                                             |
 |   - the helicity amplitudes Jem,                                                     |
 |   - the SigmaBHPol's.                                                                |
 | All quantities are evaluated at leading order in the 1/Q expansion (see Belitsky,    |
 | Mueller and Kirchner, Nucl.Phys. B629 (2002) 323-392, ArXiv:hep-ph/0112108v2)        |
 *--------------------------------------------------------------------------------------*/
 
void TGVKelly::MakeLeadingBHCrossSections()
{
	if ( InitLeadingBHCrossSections == kFALSE )
	{
		Double_t F1; // Dirac form factor
		Double_t F2; // Pauli form factor
		Double_t Ge, Gm; // Sachs' parametrization	
	
		F1=(4.*TMath::Power(M,2) - 2.79285*t)/
   (TMath::Power(1. - 1.4084507042253522*t,2)*(4.*TMath::Power(M,2) - 1.*t)); 
		F2=(7.1714*TMath::Power(M,2))/(TMath::Power(1 - 1.4084507042253522*t,2)*(4*TMath::Power(M,2) - t)); 
		Gm=2.79285/TMath::Power(1 - 1.4084507042253522*t,2);
		Ge=TMath::Power(1 - 1.4084507042253522*t,-2); 
	
	

/*----------------- Helicity amplitudes of the interference process --------------------*/

        Jem[0][2]=-((TMath::Sqrt(2.)*(F2 - Gm)*qpPerp*(-2 + xB))/(TMath::Sqrt(1 - xB)*xB));
        Jem[0][2]=-(((F2*TMath::Power(-2 + xB,2) + 4*Gm*(-1 + xB))*
       TMath::Sqrt(TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2)))/((-1 + xB)*xB));
        Jem[0][2]=(TMath::Sqrt(2.)*(F2 - Gm)*qpPerp*(-2 + xB))/(TMath::Sqrt(1 - xB)*xB);
        Jem[1][2]=(TMath::Sqrt(2.)*Gm*M*xB)/TMath::Sqrt(1 - xB);
        Jem[1][2]=0;
        Jem[1][2]=(TMath::Sqrt(2.)*Gm*M*xB)/TMath::Sqrt(1 - xB);
        Jem[2][2]=-((TMath::Sqrt(2.)*(F2*TMath::Power(qpPerp,2) + Gm*TMath::Power(M,2)*TMath::Power(xB,2)))/
     (M*TMath::Sqrt(1 - xB)*xB));
        Jem[2][2]=-((F2*qpPerp*(-2 + xB)*TMath::Sqrt(TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2)))/
     (M*(-1 + xB)*xB));
        Jem[2][2]=(TMath::Sqrt(2.)*(F2*TMath::Power(qpPerp,2) + Gm*TMath::Power(M,2)*TMath::Power(xB,2)))/(M*TMath::Sqrt(1 - xB)*xB);
        Jem[3][2]=(TMath::Sqrt(2.)*Gm*qpPerp)/TMath::Sqrt(1 - xB);
        Jem[3][2]=0;
        Jem[3][2]=(TMath::Sqrt(2.)*Gm*qpPerp)/TMath::Sqrt(1 - xB);
    	
	
/*------------------------------- BH cross sections ------------------------------------*/

        SigmaBHPol0[0]=(-32*TMath::Power(Ge,2)*TMath::Power(M,2)*(3*TMath::Power(M,8)*t + 
        3*TMath::Power(M,6)*(TMath::Power(Q,4) + TMath::Power(Q,2)*t - 4*s*t) + 
        TMath::Power(M,2)*(3*TMath::Power(Q,4)*TMath::Power(TMath::Power(Q,2) + s,2) - 
           (TMath::Power(Q,2) + s)*(TMath::Power(Q,4) + 3*TMath::Power(Q,2)*s + 12*TMath::Power(s,2))*t + 
           2*(TMath::Power(Q,4) + TMath::Power(Q,2)*s - 3*TMath::Power(s,2))*TMath::Power(t,2)) + 
        TMath::Power(TMath::Power(Q,2) + s,2)*t*(3*s*(s + t) + TMath::Power(Q,2)*(3*s + t)) + 
        TMath::Power(M,4)*(4*TMath::Power(Q,6) + TMath::Power(Q,2)*t*(3*s + t) - 
           TMath::Power(Q,4)*(6*s + t) + 3*s*t*(6*s + t))) + 
     4*TMath::Power(Gm,2)*t*(6*TMath::Power(M,8)*t + 
        t*(3*TMath::Power(TMath::Power(Q,2) + s,2)*
            (TMath::Power(Q,4) + 2*TMath::Power(Q,2)*s + 2*TMath::Power(s,2)) + 
           2*s*(TMath::Power(Q,2) + s)*(2*TMath::Power(Q,2) + 3*s)*t + 
           (3*TMath::Power(Q,4) + 4*TMath::Power(Q,2)*s + 3*TMath::Power(s,2))*TMath::Power(t,2)) - 
        2*TMath::Power(M,6)*(3*TMath::Power(Q,4) - 11*TMath::Power(Q,2)*t + 6*t*(2*s + t)) + 
        TMath::Power(M,4)*(-8*TMath::Power(Q,6) - 26*TMath::Power(Q,2)*t*(s + t) + 
           TMath::Power(Q,4)*(12*s + 25*t) + 
           3*t*(12*TMath::Power(s,2) + 10*s*t + TMath::Power(t,2))) - 
        2*TMath::Power(M,2)*(3*TMath::Power(Q,8) + TMath::Power(Q,6)*(6*s - 5*t) + 
           3*s*t*TMath::Power(2*s + t,2) + 
           TMath::Power(Q,2)*t*(7*TMath::Power(s,2) + 2*s*t - 3*TMath::Power(t,2)) + 
           TMath::Power(Q,4)*(3*TMath::Power(s,2) - 5*s*t + 7*TMath::Power(t,2)))))/
   ((TMath::Power(M,4) + 2*TMath::Power(M,2)*(TMath::Power(Q,2) - s) + TMath::Power(TMath::Power(Q,2) + s,2))*
     (4*TMath::Power(M,2) - t));
        SigmaBHPol0[1]=(-32*TMath::Power(Ge,2)*TMath::Power(M,2)*(TMath::Power(M,8)*t + 
        TMath::Power(M,6)*(TMath::Power(Q,4) + TMath::Power(Q,2)*t - 4*s*t) + 
        TMath::Power(M,2)*(TMath::Power(Q,4)*TMath::Power(TMath::Power(Q,2) + s,2) + 
           (TMath::Power(Q,2) - s)*(TMath::Power(Q,2) + s)*(5*TMath::Power(Q,2) + 4*s)*t - 
           2*(TMath::Power(Q,4) + TMath::Power(Q,2)*s + TMath::Power(s,2))*TMath::Power(t,2)) + 
        TMath::Power(TMath::Power(Q,2) + s,2)*t*(TMath::Power(Q,2)*(s - t) + s*(s + t)) + 
        TMath::Power(M,4)*(-4*TMath::Power(Q,6) + TMath::Power(Q,2)*(s - t)*t + s*t*(6*s + t) + 
           TMath::Power(Q,4)*(-2*s + 5*t))) + 
     4*TMath::Power(Gm,2)*t*(2*TMath::Power(M,8)*t + 
        t*(TMath::Power(TMath::Power(Q,2) + s,2)*
            (TMath::Power(Q,4) + 2*TMath::Power(Q,2)*s + 2*TMath::Power(s,2)) + 
           2*s*(-2*TMath::Power(Q,4) - TMath::Power(Q,2)*s + TMath::Power(s,2))*t + 
           (TMath::Power(Q,4) - 4*TMath::Power(Q,2)*s + TMath::Power(s,2))*TMath::Power(t,2)) - 
        2*TMath::Power(M,6)*(TMath::Power(Q,4) - 9*TMath::Power(Q,2)*t + 2*t*(2*s + t)) + 
        TMath::Power(M,4)*(8*TMath::Power(Q,6) - 2*TMath::Power(Q,2)*t*(15*s + 7*t) + 
           TMath::Power(Q,4)*(4*s + 19*t) + t*(12*TMath::Power(s,2) + 10*s*t + TMath::Power(t,2)))\
         - 2*TMath::Power(M,2)*(TMath::Power(Q,8) + TMath::Power(Q,6)*(2*s + t) + 
           s*t*TMath::Power(2*s + t,2) - 
           TMath::Power(Q,2)*t*(3*TMath::Power(s,2) + 10*s*t + TMath::Power(t,2)) + 
           TMath::Power(Q,4)*(TMath::Power(s,2) - 7*s*t + 5*TMath::Power(t,2)))))/
   ((TMath::Power(M,4) + 2*TMath::Power(M,2)*(TMath::Power(Q,2) - s) + TMath::Power(TMath::Power(Q,2) + s,2))*
     (4*TMath::Power(M,2) - t));
        SigmaBHPol0[2]=(-64*TMath::Power(Ge,2)*TMath::Power(M,2)*Q*qpPerp*(TMath::Power(M,2) - TMath::Power(Q,2) - s)*
      (2*TMath::Power(M,2)*TMath::Power(Q,2) - (TMath::Power(M,2) + TMath::Power(Q,2) + s)*t) + 
     16*TMath::Power(Gm,2)*Q*qpPerp*t*(TMath::Power(M,4)*(-2*TMath::Power(Q,2) + 3*t) + 
        t*(TMath::Power(Q,2)*(s - t) + s*(s + t)) + 
        TMath::Power(M,2)*(2*TMath::Power(Q,4) - t*(4*s + t) + TMath::Power(Q,2)*(2*s + 5*t))))/
   (TMath::Sqrt(TMath::Power(M,4) + 2*TMath::Power(M,2)*(TMath::Power(Q,2) - s) + 
       TMath::Power(TMath::Power(Q,2) + s,2))*(4*TMath::Power(M,2) - t));
        SigmaBHPol0[3]=(-64*TMath::Power(Ge,2)*TMath::Power(M,4)*TMath::Power(Q,2)*
      (TMath::Power(M,4)*t + s*t*(TMath::Power(Q,2) + s + t) + 
        TMath::Power(M,2)*(TMath::Power(Q,4) + TMath::Power(Q,2)*t - 2*s*t)) - 
     8*TMath::Power(Gm,2)*TMath::Power(Q,2)*(2*TMath::Power(M,2) - t)*t*
      (TMath::Power(M,4)*t + s*t*(TMath::Power(Q,2) + s + t) + 
        TMath::Power(M,2)*(TMath::Power(Q,4) + TMath::Power(Q,2)*t - 2*s*t)))/
   ((TMath::Power(M,4) + 2*TMath::Power(M,2)*(TMath::Power(Q,2) - s) + TMath::Power(TMath::Power(Q,2) + s,2))*
     (4*TMath::Power(M,2) - t));
        SigmaBHPolX[0]=(64*Ge*Gm*M*qpPerp*(-(TMath::Power(Q,2)*(TMath::Power(Q,2) + s)) + 
        TMath::Power(M,2)*(TMath::Power(Q,2) - t) + s*t)*
      (2*TMath::Power(M,2)*TMath::Power(Q,2) - (TMath::Power(M,2) + TMath::Power(Q,2) + s)*t) + 
     32*TMath::Power(Gm,2)*M*qpPerp*t*(-TMath::Power(Q,6) - 3*TMath::Power(Q,4)*s - 
        2*TMath::Power(M,4)*(TMath::Power(Q,2) - t) + s*t*(2*s + t) - 
        TMath::Power(Q,2)*(2*TMath::Power(s,2) + TMath::Power(t,2)) - 
        TMath::Power(M,2)*(TMath::Power(Q,4) - 4*TMath::Power(Q,2)*(s + t) + t*(4*s + t))))/
   (TMath::Sqrt(TMath::Power(M,4) + 2*TMath::Power(M,2)*(TMath::Power(Q,2) - s) + 
       TMath::Power(TMath::Power(Q,2) + s,2))*(4*TMath::Power(M,2) - t));
        SigmaBHPolX[1]=(32*Ge*Gm*M*Q*(-2*TMath::Power(M,2)*TMath::Power(Q,2) + (TMath::Power(M,2) + TMath::Power(Q,2) + s)*t)*
      (4*TMath::Power(M,2)*TMath::Power(Q,4) + 
        (2*TMath::Power(M,2) + TMath::Power(Q,2) - 2*s)*(TMath::Power(M,2) - TMath::Power(Q,2) - s)*t + 
        (TMath::Power(M,2) + TMath::Power(Q,2) + 3*s)*TMath::Power(t,2)) + 
     64*TMath::Power(Gm,2)*M*Q*(TMath::Power(M,2) - s - t)*t*
      (TMath::Power(M,4)*t + s*t*(TMath::Power(Q,2) + s + t) + 
        TMath::Power(M,2)*(TMath::Power(Q,4) + TMath::Power(Q,2)*t - 2*s*t)))/
   ((TMath::Power(M,4) + 2*TMath::Power(M,2)*(TMath::Power(Q,2) - s) + TMath::Power(TMath::Power(Q,2) + s,2))*
     (4*TMath::Power(M,2) - t));
        SigmaBHPolY=32*Ge*Gm*M*Q*(TMath::Power(Q,2) - t)*t;
        SigmaBHPolZ[0]=(-128*Ge*Gm*TMath::Power(M,2)*(-(TMath::Power(Q,2)*(TMath::Power(Q,2) + s)) + 
        TMath::Power(M,2)*(TMath::Power(Q,2) - t) + s*t)*
      (TMath::Power(M,4)*t + s*t*(TMath::Power(Q,2) + s + t) + 
        TMath::Power(M,2)*(TMath::Power(Q,4) + TMath::Power(Q,2)*t - 2*s*t)) + 
     16*TMath::Power(Gm,2)*t*(2*TMath::Power(M,2)*TMath::Power(Q,2) - 
        (TMath::Power(M,2) + TMath::Power(Q,2) + s)*t)*
      (TMath::Power(Q,6) + 3*TMath::Power(Q,4)*s + 2*TMath::Power(M,4)*(TMath::Power(Q,2) - t) - 
        s*t*(2*s + t) + TMath::Power(Q,2)*(2*TMath::Power(s,2) + TMath::Power(t,2)) + 
        TMath::Power(M,2)*(TMath::Power(Q,4) - 4*TMath::Power(Q,2)*(s + t) + t*(4*s + t))))/
   ((TMath::Power(M,4) + 2*TMath::Power(M,2)*(TMath::Power(Q,2) - s) + TMath::Power(TMath::Power(Q,2) + s,2))*
     (4*TMath::Power(M,2) - t));
        SigmaBHPolZ[1]=(64*Ge*Gm*TMath::Power(M,2)*Q*qpPerp*(-4*TMath::Power(M,2)*TMath::Power(Q,4) - 
        (2*TMath::Power(M,2) + TMath::Power(Q,2) - 2*s)*(TMath::Power(M,2) - TMath::Power(Q,2) - s)*t - 
        (TMath::Power(M,2) + TMath::Power(Q,2) + 3*s)*TMath::Power(t,2)) + 
     32*TMath::Power(Gm,2)*Q*qpPerp*t*(-TMath::Power(M,2) + s + t)*
      ((TMath::Power(Q,2) + s)*t + TMath::Power(M,2)*(-2*TMath::Power(Q,2) + t)))/
   (TMath::Sqrt(TMath::Power(M,4) + 2*TMath::Power(M,2)*(TMath::Power(Q,2) - s) + 
       TMath::Power(TMath::Power(Q,2) + s,2))*(4*TMath::Power(M,2) - t));
    
    
    	// Flag
    	
    	InitLeadingBHCrossSections = kTRUE;    	
    	
    	
    	// Validation mode only
    
   	 	if ( Validation == kTRUE )
   	 	{	
    		ofstream outfile;    		
    		
    		outfile.open("CheckFormFactors.dat");
			outfile << F1 << endl;
			outfile <<  F2 << endl;
			outfile <<  Ge << endl;
			outfile << Gm << endl;
			outfile.close();
		
		
    		outfile.open("CheckJemLeading.dat");
    		for(Int_t i=0;i<12;i++)
			{
				outfile << Jem[div(i,4).rem][div(i,4).quot] << endl;	
			} // end for i		
			outfile.close();
		
		
    		outfile.open("CheckSigmaBHLeading.dat");
    		for(Int_t i=0;i<9;i++)
			{
				if (0 <= i && i < 4)
				{
					outfile << SigmaBHPol0[i] << endl;
				} // end if i
				if (4 <= i && i < 6)
				{
					outfile << SigmaBHPolX[i-4] << endl;
				} // end if i
				if (i == 6)
				{
					outfile << SigmaBHPolY << endl;
				} // end if i
				if (7<= i && i < 9)
				{
					outfile << SigmaBHPolZ[i-7] << endl;
				} // end if i
			} // end for i
			outfile.close();		
				
		} // end if Validation  
		
	} // end if InitLeadingBHCrossSections
	 
} // end MakeLeadingBHCrossSections


	
/*------------------- Function MakeLeadingVCSAndInterfCrossSections() ------------------*
 | Computes all the stuff to evaluate the cross section assuming the hadronic helicity  |
 | amplitudes are given, i.e. initializes :                                             |
 |   - the helicity amplitudes Mvcs,                                                    |
 |   - the expansion coefficients Ur,                                                   |
 |   - the SigmaVCSPol's and SigmaIPol's.                                               |
 | All quantities are evaluated at leading order in the 1/Q expansion (see Belitsky,    |
 | Mueller and Kirchner, Nucl.Phys. B629 (2002) 323-392, ArXiv:hep-ph/0112108v2)        |
 *--------------------------------------------------------------------------------------*/
 
void TGVKelly::MakeLeadingVCSAndInterfCrossSections()
{
//	if ( InitLeadingVCSAndInterfCrossSections == kFALSE )
//	{

/*-------------- Harmonic expansion coefficients of the VCS cross section --------------*/

        Ur[0]=(-13*TMath::Power(qpPerp,2)*TMath::Sqrt((TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/
       (1 - xB)))/16.;
        Ur[1]=TMath::Sqrt((TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/(1 - xB))/16.;
        Ur[2]=(-3*TMath::Power(qpPerp,2)*TMath::Sqrt((TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/
       (1 - xB)))/16.;
        Ur[3]=(7*TMath::Sqrt((TMath::Power(qpPerp,4) + TMath::Power(M,2)*TMath::Power(qpPerp,2)*TMath::Power(xB,2))/(1 - xB)))/
   16.;
        Ur[4]=(3*TMath::Sqrt(-((-1 + xB)*(TMath::Power(qpPerp,4) + 
           TMath::Power(M,2)*TMath::Power(qpPerp,2)*TMath::Power(xB,2)))))/(16.*xB);
        Ur[5]=(3*TMath::Sqrt((TMath::Power(qpPerp,4) + TMath::Power(M,2)*TMath::Power(qpPerp,2)*TMath::Power(xB,2))/(1 - xB)))/
   16.;
        Ur[6]=-TMath::Sqrt((TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/(1 - xB))/8.;
        Ur[7]=-TMath::Sqrt((TMath::Power(qpPerp,4) + TMath::Power(M,2)*TMath::Power(qpPerp,2)*TMath::Power(xB,2))/(1 - xB))/2.;
        Ur[8]=(TMath::Power(qpPerp,3)*(13 - 6*xB) + 7*TMath::Power(M,2)*qpPerp*TMath::Power(xB,2))/
   (8.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[9]=qpPerp/(4.*TMath::Sqrt(2.));
        Ur[10]=(TMath::Power(qpPerp,3)*(3 - 2*xB) + TMath::Power(M,2)*qpPerp*TMath::Power(xB,2))/
   (8.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[11]=(-5*TMath::Power(M,2)*TMath::Power(xB,2) + TMath::Power(qpPerp,2)*(-11 + 6*xB))/
   (16.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[12]=-3/(16.*TMath::Sqrt(2.));
        Ur[13]=(-(TMath::Power(M,2)*TMath::Power(xB,2)) + TMath::Power(qpPerp,2)*(-7 + 6*xB))/
   (16.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[14]=-qpPerp/(4.*TMath::Sqrt(2.));
        Ur[15]=(TMath::Power(qpPerp,2)*(4 - 3*xB) + TMath::Power(M,2)*TMath::Power(xB,2))/(4.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[16]=-qpPerp/(4.*TMath::Sqrt(2.));
        Ur[17]=(TMath::Power(qpPerp,2)*(1 - 6*xB) - 5*TMath::Power(M,2)*TMath::Power(xB,2))/
   (16.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[18]=3/(16.*TMath::Sqrt(2.));
        Ur[19]=(TMath::Power(qpPerp,2)*(5 - 6*xB) - TMath::Power(M,2)*TMath::Power(xB,2))/(16.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[20]=-(qpPerp*(2*TMath::Power(M,2)*TMath::Power(xB,2) + TMath::Power(qpPerp,2)*(1 + xB)))/
   (4.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[21]=qpPerp/(4.*TMath::Sqrt(2.));
        Ur[22]=(TMath::Power(qpPerp,3)*(1 - 3*xB) - 2*TMath::Power(M,2)*qpPerp*TMath::Power(xB,2))/
   (4.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[23]=(TMath::Power(M,2)*TMath::Power(xB,2) + TMath::Power(qpPerp,2)*(-2 + 3*xB))/(4.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[24]=(3*TMath::Sqrt((TMath::Power(qpPerp,4) + TMath::Power(M,2)*TMath::Power(qpPerp,2)*TMath::Power(xB,2))/
       (2 - 2*xB)))/8.;
        Ur[25]=-TMath::Sqrt(-((-1 + xB)*(TMath::Power(qpPerp,4) + 
          TMath::Power(M,2)*TMath::Power(qpPerp,2)*TMath::Power(xB,2))))/(8.*TMath::Sqrt(2.)*xB);
        Ur[26]=(-3*TMath::Sqrt((TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/(2 - 2*xB)))/8.;
        Ur[27]=TMath::Sqrt((TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/(2 - 2*xB))/8.;
        Ur[28]=(TMath::Power(qpPerp,2)*(-1 + 2*xB)*TMath::Sqrt((TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/
       (2 - 2*xB)))/(2.*xB);
        Ur[29]=-(TMath::Power(qpPerp,2)*(7 - 6*xB) + TMath::Power(M,2)*TMath::Power(xB,2))/(16.*(-1 + xB));
        Ur[30]=0.0625;
        Ur[31]=-qpPerp/4.;
        Ur[32]=qpPerp/4.;
        Ur[33]=qpPerp/4.;
        Ur[34]=-0.0625;
        Ur[35]=qpPerp/4.;
        Ur[36]=-qpPerp/4.;
        Ur[37]=(TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/(-8 + 8*xB);
        Ur[38]=-qpPerp/4.;
        Ur[39]=-TMath::Sqrt((TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/(1 - xB))/16.;
        Ur[40]=(qpPerp*(5 - 9*xB + 4*TMath::Power(xB,2))*
     TMath::Sqrt(-((TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/TMath::Power(-1 + xB,3))))/
   (16.*xB);
        Ur[41]=(3*TMath::Power(qpPerp,3)*TMath::Sqrt((TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/(1 - xB)))/
   16.;
        Ur[42]=((1 - 4*xB)*TMath::Sqrt((TMath::Power(qpPerp,4) + TMath::Power(M,2)*TMath::Power(qpPerp,2)*TMath::Power(xB,2))/
       (1 - xB)))/(16.*xB);
        Ur[43]=(TMath::Power(qpPerp,2)*TMath::Sqrt((TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/(1 - xB)))/8.;
        Ur[44]=TMath::Sqrt((TMath::Power(qpPerp,4) + TMath::Power(M,2)*TMath::Power(qpPerp,2)*TMath::Power(xB,2))/(1 - xB))/
   (4.*xB);
        Ur[45]=(3*qpPerp)/(4.*TMath::Sqrt(2.));
        Ur[46]=-5/(16.*TMath::Sqrt(2.));
        Ur[47]=(3*TMath::Power(qpPerp,4)*(-2 + xB) - 3*TMath::Power(M,2)*TMath::Power(qpPerp,2)*TMath::Power(xB,2))/
   (16.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[48]=-1/(16.*TMath::Sqrt(2.));
        Ur[49]=TMath::Power(qpPerp,3)/(4.*TMath::Sqrt(2.));
        Ur[50]=-1/(4.*TMath::Sqrt(2.));
        Ur[51]=(TMath::Power(qpPerp,3)*(1 - 2*xB) - TMath::Power(M,2)*qpPerp*TMath::Power(xB,2))/
   (8.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[52]=5/(16.*TMath::Sqrt(2.));
        Ur[53]=(3*xB*(TMath::Power(qpPerp,4) + TMath::Power(M,2)*TMath::Power(qpPerp,2)*xB))/
   (16.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[54]=1/(16.*TMath::Sqrt(2.));
        Ur[55]=-qpPerp/(4.*TMath::Sqrt(2.));
        Ur[56]=TMath::Power(qpPerp,3)/(4.*TMath::Sqrt(2.));
        Ur[57]=(-3*qpPerp)/(4.*TMath::Sqrt(2.));
        Ur[58]=1/(4.*TMath::Sqrt(2.));
        Ur[59]=(-7*TMath::Sqrt((TMath::Power(qpPerp,4) + TMath::Power(M,2)*TMath::Power(qpPerp,2)*TMath::Power(xB,2))/
       (1 - xB)))/16.;
        Ur[60]=(-3*TMath::Sqrt((TMath::Power(qpPerp,4) + TMath::Power(M,2)*TMath::Power(qpPerp,2)*TMath::Power(xB,2))/
       (1 - xB)))/16.;
        Ur[61]=(-3*TMath::Power(qpPerp,2)*TMath::Sqrt((TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/
       (1 - xB)))/8.;
        Ur[62]=(-5*TMath::Power(qpPerp,2)*TMath::Sqrt((TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/
       (1 - xB)))/8.;
        Ur[63]=TMath::Sqrt((TMath::Power(qpPerp,4) + TMath::Power(M,2)*TMath::Power(qpPerp,2)*TMath::Power(xB,2))/(1 - xB))/2.;
        Ur[64]=(TMath::Power(qpPerp,2)*(11 - 6*xB) + 5*TMath::Power(M,2)*TMath::Power(xB,2))/
   (16.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[65]=(TMath::Power(qpPerp,2)*(7 - 6*xB) + TMath::Power(M,2)*TMath::Power(xB,2))/(16.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[66]=(-(TMath::Power(qpPerp,3)*(-3 + xB)) + 2*TMath::Power(M,2)*qpPerp*TMath::Power(xB,2))/
   (4.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[67]=(TMath::Power(qpPerp,3)*(5 - 3*xB) + 2*TMath::Power(M,2)*qpPerp*TMath::Power(xB,2))/
   (4.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[68]=(-(TMath::Power(M,2)*TMath::Power(xB,2)) + TMath::Power(qpPerp,2)*(-4 + 3*xB))/
   (4.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[69]=-(qpPerp*(7*TMath::Power(M,2)*TMath::Power(xB,2) + TMath::Power(qpPerp,2)*(1 + 6*xB)))/
   (8.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[70]=(5*TMath::Power(M,2)*TMath::Power(xB,2) + TMath::Power(qpPerp,2)*(-1 + 6*xB))/
   (16.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[71]=(TMath::Power(M,2)*TMath::Power(xB,2) + TMath::Power(qpPerp,2)*(-5 + 6*xB))/(16.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[72]=(TMath::Power(qpPerp,2)*(2 - 3*xB) - TMath::Power(M,2)*TMath::Power(xB,2))/(4.*TMath::Sqrt(2.)*(-1 + xB));
        Ur[73]=TMath::Sqrt(-((-1 + xB)*(TMath::Power(qpPerp,4) + TMath::Power(M,2)*TMath::Power(qpPerp,2)*TMath::Power(xB,2))))/
   (8.*TMath::Sqrt(2.)*xB);
        Ur[74]=(3*TMath::Sqrt((TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/(2 - 2*xB)))/8.;
        Ur[75]=-TMath::Sqrt((TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/(2 - 2*xB))/8.;
        Ur[76]=TMath::Sqrt((TMath::Power(qpPerp,4) + TMath::Power(M,2)*TMath::Power(qpPerp,2)*TMath::Power(xB,2))/(2 - 2*xB))/4.;
        Ur[77]=(TMath::Power(qpPerp,2)*TMath::Sqrt((TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/(2 - 2*xB)))/
   (2.*xB);
        Ur[78]=-0.0625;
        Ur[79]=qpPerp/4.;
        Ur[80]=-qpPerp/4.;
        Ur[81]=(TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/(8 - 8*xB);
        Ur[82]=-qpPerp/4.;
        Ur[83]=(TMath::Power(M,2)*TMath::Power(xB,2) + TMath::Power(qpPerp,2)*(-5 + 6*xB))/(16.*(-1 + xB));
        Ur[84]=0.0625;
        Ur[85]=-qpPerp/4.;
        Ur[86]=-(qpPerp*(5 - 11*xB + 6*TMath::Power(xB,2))*
      TMath::Sqrt(-((TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/TMath::Power(-1 + xB,3))))/
   (16.*xB);
        Ur[87]=-((1 + 2*xB)*TMath::Sqrt((TMath::Power(qpPerp,4) + TMath::Power(M,2)*TMath::Power(qpPerp,2)*TMath::Power(xB,2))/
        (1 - xB)))/(16.*xB);
        Ur[88]=TMath::Sqrt((TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/(1 - xB))/8.;
        Ur[89]=-TMath::Sqrt((TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/(1 - xB))/8.;
        Ur[90]=-(qpPerp*(1 - 3*xB + 2*TMath::Power(xB,2))*
      TMath::Sqrt(-((TMath::Power(qpPerp,2) + TMath::Power(M,2)*TMath::Power(xB,2))/TMath::Power(-1 + xB,3))))/
   (4.*xB);
        Ur[91]=5/(16.*TMath::Sqrt(2.));
        Ur[92]=1/(16.*TMath::Sqrt(2.));
        Ur[93]=-qpPerp/(4.*TMath::Sqrt(2.));
        Ur[94]=(-3*qpPerp)/(4.*TMath::Sqrt(2.));
        Ur[95]=1/(4.*TMath::Sqrt(2.));
        Ur[96]=(3*qpPerp)/(4.*TMath::Sqrt(2.));
        Ur[97]=-5/(16.*TMath::Sqrt(2.));
        Ur[98]=-1/(16.*TMath::Sqrt(2.));
        Ur[99]=-1/(4.*TMath::Sqrt(2.));
	
	
/*-------------------------------- VCS cross sections ----------------------------------*/

        SigmaVCSPol0[0]=(3*TMath::Power(IMvcs[0][0],2) - 2*TMath::Power(IMvcs[0][1],2) + 
     3*TMath::Power(IMvcs[0][2],2) + 3*TMath::Power(IMvcs[1][0],2) - 
     2*TMath::Power(IMvcs[1][1],2) + 3*TMath::Power(IMvcs[1][2],2) + 
     3*TMath::Power(IMvcs[2][0],2) - 2*TMath::Power(IMvcs[2][1],2) + 
     3*TMath::Power(IMvcs[2][2],2) + 3*TMath::Power(IMvcs[3][0],2) - 
     2*TMath::Power(IMvcs[3][1],2) + 3*TMath::Power(IMvcs[3][2],2) + 
     3*TMath::Power(RMvcs[0][0],2) - 2*TMath::Power(RMvcs[0][1],2) + 
     3*TMath::Power(RMvcs[0][2],2) + 3*TMath::Power(RMvcs[1][0],2) - 
     2*TMath::Power(RMvcs[1][1],2) + 3*TMath::Power(RMvcs[1][2],2) + 
     3*TMath::Power(RMvcs[2][0],2) - 2*TMath::Power(RMvcs[2][1],2) + 
     3*TMath::Power(RMvcs[2][2],2) + 3*TMath::Power(RMvcs[3][0],2) - 
     2*TMath::Power(RMvcs[3][1],2) + 3*TMath::Power(RMvcs[3][2],2))/4.;
        SigmaVCSPol0[1]=(TMath::Power(IMvcs[0][0],2) + 2*TMath::Power(IMvcs[0][1],2) + 
     TMath::Power(IMvcs[0][2],2) + TMath::Power(IMvcs[1][0],2) + 
     2*TMath::Power(IMvcs[1][1],2) + TMath::Power(IMvcs[1][2],2) + 
     TMath::Power(IMvcs[2][0],2) + 2*TMath::Power(IMvcs[2][1],2) + 
     TMath::Power(IMvcs[2][2],2) + TMath::Power(IMvcs[3][0],2) + 
     2*TMath::Power(IMvcs[3][1],2) + TMath::Power(IMvcs[3][2],2) + 
     TMath::Power(RMvcs[0][0],2) + 2*TMath::Power(RMvcs[0][1],2) + 
     TMath::Power(RMvcs[0][2],2) + TMath::Power(RMvcs[1][0],2) + 
     2*TMath::Power(RMvcs[1][1],2) + TMath::Power(RMvcs[1][2],2) + 
     TMath::Power(RMvcs[2][0],2) + 2*TMath::Power(RMvcs[2][1],2) + 
     TMath::Power(RMvcs[2][2],2) + TMath::Power(RMvcs[3][0],2) + 
     2*TMath::Power(RMvcs[3][1],2) + TMath::Power(RMvcs[3][2],2))/4.;
        SigmaVCSPol0[2]=(IMvcs[0][1]*(IMvcs[0][0] - IMvcs[0][2]) + 
     IMvcs[1][1]*(IMvcs[1][0] - IMvcs[1][2]) + 
     IMvcs[2][1]*(IMvcs[2][0] - IMvcs[2][2]) + 
     IMvcs[3][1]*(IMvcs[3][0] - IMvcs[3][2]) + 
     RMvcs[0][1]*(RMvcs[0][0] - RMvcs[0][2]) + 
     RMvcs[1][1]*(RMvcs[1][0] - RMvcs[1][2]) + 
     RMvcs[2][1]*(RMvcs[2][0] - RMvcs[2][2]) + 
     RMvcs[3][1]*(RMvcs[3][0] - RMvcs[3][2]))/TMath::Sqrt(2.);
        SigmaVCSPol0[3]=(IMvcs[0][0]*IMvcs[0][2] + IMvcs[1][0]*IMvcs[1][2] + 
     IMvcs[2][0]*IMvcs[2][2] + IMvcs[3][0]*IMvcs[3][2] + 
     RMvcs[0][0]*RMvcs[0][2] + RMvcs[1][0]*RMvcs[1][2] + 
     RMvcs[2][0]*RMvcs[2][2] + RMvcs[3][0]*RMvcs[3][2])/2.;
        SigmaVCSPol0[4]=TMath::Sqrt(2.)*((IMvcs[0][0] - IMvcs[0][2])*RMvcs[0][1] + 
     IMvcs[0][1]*(-RMvcs[0][0] + RMvcs[0][2]) + 
     (IMvcs[1][0] - IMvcs[1][2])*RMvcs[1][1] + 
     IMvcs[1][1]*(-RMvcs[1][0] + RMvcs[1][2]) + 
     (IMvcs[2][0] - IMvcs[2][2])*RMvcs[2][1] + 
     IMvcs[2][1]*(-RMvcs[2][0] + RMvcs[2][2]) + 
     (IMvcs[3][0] - IMvcs[3][2])*RMvcs[3][1] + 
     IMvcs[3][1]*(-RMvcs[3][0] + RMvcs[3][2]));
        SigmaVCSPolX[0]=2*(IMvcs[0][0]*IMvcs[1][0] - IMvcs[0][2]*IMvcs[1][2] - 
     IMvcs[2][0]*IMvcs[3][0] + IMvcs[2][2]*IMvcs[3][2] + 
     RMvcs[0][0]*RMvcs[1][0] - RMvcs[0][2]*RMvcs[1][2] - 
     RMvcs[2][0]*RMvcs[3][0] + RMvcs[2][2]*RMvcs[3][2]);
        SigmaVCSPolX[1]=TMath::Sqrt(2.)*((IMvcs[0][0] + IMvcs[0][2])*IMvcs[1][1] + 
     IMvcs[0][1]*(IMvcs[1][0] + IMvcs[1][2]) - 
     (IMvcs[2][0] + IMvcs[2][2])*IMvcs[3][1] - 
     IMvcs[2][1]*(IMvcs[3][0] + IMvcs[3][2]) + 
     (RMvcs[0][0] + RMvcs[0][2])*RMvcs[1][1] + 
     RMvcs[0][1]*(RMvcs[1][0] + RMvcs[1][2]) - 
     (RMvcs[2][0] + RMvcs[2][2])*RMvcs[3][1] - 
     RMvcs[2][1]*(RMvcs[3][0] + RMvcs[3][2]));
        SigmaVCSPolX[2]=((IMvcs[1][0] + IMvcs[1][2])*RMvcs[0][1] - 
     IMvcs[1][1]*(RMvcs[0][0] + RMvcs[0][2]) + 
     (IMvcs[0][0] + IMvcs[0][2])*RMvcs[1][1] - 
     IMvcs[0][1]*(RMvcs[1][0] + RMvcs[1][2]) - 
     (IMvcs[3][0] + IMvcs[3][2])*RMvcs[2][1] + 
     IMvcs[3][1]*(RMvcs[2][0] + RMvcs[2][2]) - 
     (IMvcs[2][0] + IMvcs[2][2])*RMvcs[3][1] + 
     IMvcs[2][1]*(RMvcs[3][0] + RMvcs[3][2]))/TMath::Sqrt(2.);
        SigmaVCSPolX[3]=(-(IMvcs[1][2]*RMvcs[0][0]) + IMvcs[1][0]*RMvcs[0][2] - 
     IMvcs[0][2]*RMvcs[1][0] + IMvcs[0][0]*RMvcs[1][2] + 
     IMvcs[3][2]*RMvcs[2][0] - IMvcs[3][0]*RMvcs[2][2] + 
     IMvcs[2][2]*RMvcs[3][0] - IMvcs[2][0]*RMvcs[3][2])/2.;
        SigmaVCSPolY[0]=TMath::Sqrt(2.)*((IMvcs[0][0] - IMvcs[0][2])*IMvcs[2][1] + 
     IMvcs[0][1]*(-IMvcs[2][0] + IMvcs[2][2]) + 
     (IMvcs[1][0] - IMvcs[1][2])*IMvcs[3][1] + 
     IMvcs[1][1]*(-IMvcs[3][0] + IMvcs[3][2]) + 
     (RMvcs[0][0] - RMvcs[0][2])*RMvcs[2][1] + 
     RMvcs[0][1]*(-RMvcs[2][0] + RMvcs[2][2]) + 
     (RMvcs[1][0] - RMvcs[1][2])*RMvcs[3][1] + 
     RMvcs[1][1]*(-RMvcs[3][0] + RMvcs[3][2]));
        SigmaVCSPolY[1]=(3*IMvcs[2][0]*RMvcs[0][0] - 2*IMvcs[2][1]*RMvcs[0][1] + 
     3*IMvcs[2][2]*RMvcs[0][2] + 3*IMvcs[3][0]*RMvcs[1][0] - 
     2*IMvcs[3][1]*RMvcs[1][1] + 3*IMvcs[3][2]*RMvcs[1][2] - 
     3*IMvcs[0][0]*RMvcs[2][0] + 2*IMvcs[0][1]*RMvcs[2][1] - 
     3*IMvcs[0][2]*RMvcs[2][2] - 3*IMvcs[1][0]*RMvcs[3][0] + 
     2*IMvcs[1][1]*RMvcs[3][1] - 3*IMvcs[1][2]*RMvcs[3][2])/2.;
        SigmaVCSPolY[2]=(IMvcs[2][0]*RMvcs[0][0] + 2*IMvcs[2][1]*RMvcs[0][1] + 
     IMvcs[2][2]*RMvcs[0][2] + IMvcs[3][0]*RMvcs[1][0] + 
     2*IMvcs[3][1]*RMvcs[1][1] + IMvcs[3][2]*RMvcs[1][2] - 
     IMvcs[0][0]*RMvcs[2][0] - 2*IMvcs[0][1]*RMvcs[2][1] - 
     IMvcs[0][2]*RMvcs[2][2] - IMvcs[1][0]*RMvcs[3][0] - 
     2*IMvcs[1][1]*RMvcs[3][1] - IMvcs[1][2]*RMvcs[3][2])/2.;
        SigmaVCSPolY[3]=((IMvcs[2][0] - IMvcs[2][2])*RMvcs[0][1] + 
     IMvcs[2][1]*(RMvcs[0][0] - RMvcs[0][2]) + 
     (IMvcs[3][0] - IMvcs[3][2])*RMvcs[1][1] + 
     IMvcs[3][1]*(RMvcs[1][0] - RMvcs[1][2]) + 
     (-IMvcs[0][0] + IMvcs[0][2])*RMvcs[2][1] + 
     IMvcs[0][1]*(-RMvcs[2][0] + RMvcs[2][2]) + 
     (-IMvcs[1][0] + IMvcs[1][2])*RMvcs[3][1] + 
     IMvcs[1][1]*(-RMvcs[3][0] + RMvcs[3][2]))/TMath::Sqrt(2.);
        SigmaVCSPolY[4]=(IMvcs[2][2]*RMvcs[0][0] + IMvcs[2][0]*RMvcs[0][2] + 
     IMvcs[3][2]*RMvcs[1][0] + IMvcs[3][0]*RMvcs[1][2] - 
     IMvcs[0][2]*RMvcs[2][0] - IMvcs[0][0]*RMvcs[2][2] - 
     IMvcs[1][2]*RMvcs[3][0] - IMvcs[1][0]*RMvcs[3][2])/2.;
        SigmaVCSPolZ[0]=2*(IMvcs[1][0]*IMvcs[2][0] - IMvcs[1][2]*IMvcs[2][2] + 
     IMvcs[0][0]*IMvcs[3][0] - IMvcs[0][2]*IMvcs[3][2] + 
     RMvcs[1][0]*RMvcs[2][0] - RMvcs[1][2]*RMvcs[2][2] + 
     RMvcs[0][0]*RMvcs[3][0] - RMvcs[0][2]*RMvcs[3][2]);
        SigmaVCSPolZ[1]=TMath::Sqrt(2.)*((IMvcs[1][0] + IMvcs[1][2])*IMvcs[2][1] + 
     IMvcs[1][1]*(IMvcs[2][0] + IMvcs[2][2]) + 
     (IMvcs[0][0] + IMvcs[0][2])*IMvcs[3][1] + 
     IMvcs[0][1]*(IMvcs[3][0] + IMvcs[3][2]) + 
     (RMvcs[1][0] + RMvcs[1][2])*RMvcs[2][1] + 
     RMvcs[1][1]*(RMvcs[2][0] + RMvcs[2][2]) + 
     (RMvcs[0][0] + RMvcs[0][2])*RMvcs[3][1] + 
     RMvcs[0][1]*(RMvcs[3][0] + RMvcs[3][2]));
        SigmaVCSPolZ[2]=((IMvcs[3][0] + IMvcs[3][2])*RMvcs[0][1] - 
     IMvcs[3][1]*(RMvcs[0][0] + RMvcs[0][2]) + 
     (IMvcs[2][0] + IMvcs[2][2])*RMvcs[1][1] - 
     IMvcs[2][1]*(RMvcs[1][0] + RMvcs[1][2]) + 
     (IMvcs[1][0] + IMvcs[1][2])*RMvcs[2][1] - 
     IMvcs[1][1]*(RMvcs[2][0] + RMvcs[2][2]) + 
     (IMvcs[0][0] + IMvcs[0][2])*RMvcs[3][1] - 
     IMvcs[0][1]*(RMvcs[3][0] + RMvcs[3][2]))/TMath::Sqrt(2.);
        SigmaVCSPolZ[3]=(-(IMvcs[3][2]*RMvcs[0][0]) + IMvcs[3][0]*RMvcs[0][2] - 
     IMvcs[2][2]*RMvcs[1][0] + IMvcs[2][0]*RMvcs[1][2] - 
     IMvcs[1][2]*RMvcs[2][0] + IMvcs[1][0]*RMvcs[2][2] - 
     IMvcs[0][2]*RMvcs[3][0] + IMvcs[0][0]*RMvcs[3][2])/2.;
		
	
/*--------------------------- Interference cross sections ------------------------------*/

		if ( InitLeadingBHCrossSections == kFALSE )
		{
			TGVKelly::MakeLeadingBHCrossSections();
			
			
			// Flag
			
			InitLeadingBHCrossSections = kTRUE;
		} // end if InitLeadingBHCrossSections
		

        SigmaIPol0[0]=Jem[0][1]*(-2*RMvcs[0][0]*Ur[0]*TMath::Power(Q,2) + 
      2*RMvcs[0][1]*Ur[24]*TMath::Power(Q,3) - 
      2*RMvcs[0][2]*Ur[39]*TMath::Power(Q,4)) + 
   Jem[2][1]*(-2*RMvcs[2][0]*Ur[0]*TMath::Power(Q,2) + 
      2*RMvcs[2][1]*Ur[24]*TMath::Power(Q,3) - 
      2*RMvcs[2][2]*Ur[39]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*RMvcs[0][0]*Ur[8]*TMath::Power(Q,2) + 
      2*RMvcs[0][1]*Ur[29]*TMath::Power(Q,3) - 
      2*RMvcs[0][2]*Ur[45]*TMath::Power(Q,4)) + 
   Jem[2][2]*(-2*RMvcs[2][0]*Ur[8]*TMath::Power(Q,2) + 
      2*RMvcs[2][1]*Ur[29]*TMath::Power(Q,3) - 
      2*RMvcs[2][2]*Ur[45]*TMath::Power(Q,4)) + 
   Jem[1][2]*(-2*RMvcs[1][0]*Ur[69]*TMath::Power(Q,2) + 
      2*RMvcs[1][1]*Ur[83]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[96]*TMath::Power(Q,4)) + 
   Jem[3][2]*(-2*RMvcs[3][0]*Ur[69]*TMath::Power(Q,2) + 
      2*RMvcs[3][1]*Ur[83]*TMath::Power(Q,3) - 
      2*RMvcs[3][2]*Ur[96]*TMath::Power(Q,4));
        SigmaIPol0[1]=Jem[0][1]*(-2*RMvcs[0][0]*Ur[2]*TMath::Power(Q,2) - 
      2*RMvcs[0][1]*Ur[24]*TMath::Power(Q,3) - 
      2*RMvcs[0][2]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[2][1]*(-2*RMvcs[2][0]*Ur[2]*TMath::Power(Q,2) - 
      2*RMvcs[2][1]*Ur[24]*TMath::Power(Q,3) - 
      2*RMvcs[2][2]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[1][2]*(-2*RMvcs[1][0]*Ur[51]*TMath::Power(Q,2) - 
      2*RMvcs[1][1]*Ur[83]*TMath::Power(Q,3) + 
      2*RMvcs[1][2]*Ur[16]*TMath::Power(Q,4)) + 
   Jem[3][2]*(-2*RMvcs[3][0]*Ur[51]*TMath::Power(Q,2) - 
      2*RMvcs[3][1]*Ur[83]*TMath::Power(Q,3) + 
      2*RMvcs[3][2]*Ur[16]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*RMvcs[0][0]*Ur[10]*TMath::Power(Q,2) - 
      2*RMvcs[0][1]*Ur[29]*TMath::Power(Q,3) - 
      2*RMvcs[0][2]*Ur[9]*TMath::Power(Q,4)) + 
   Jem[2][2]*(-2*RMvcs[2][0]*Ur[10]*TMath::Power(Q,2) - 
      2*RMvcs[2][1]*Ur[29]*TMath::Power(Q,3) - 
      2*RMvcs[2][2]*Ur[9]*TMath::Power(Q,4));
        SigmaIPol0[2]=Jem[0][1]*(-2*RMvcs[0][0]*Ur[3]*TMath::Power(Q,3) - 
      2*RMvcs[0][2]*Ur[40]*TMath::Power(Q,3) + 
      2*RMvcs[0][1]*Ur[26]*TMath::Power(Q,4)) + 
   Jem[2][1]*(-2*RMvcs[2][0]*Ur[3]*TMath::Power(Q,3) - 
      2*RMvcs[2][2]*Ur[40]*TMath::Power(Q,3) + 
      2*RMvcs[2][1]*Ur[26]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*RMvcs[0][0]*Ur[11]*TMath::Power(Q,3) + 
      2*RMvcs[0][1]*Ur[31]*TMath::Power(Q,4) - 
      2*RMvcs[0][2]*Ur[46]*TMath::Power(Q,5)) + 
   Jem[2][2]*(-2*RMvcs[2][0]*Ur[11]*TMath::Power(Q,3) + 
      2*RMvcs[2][1]*Ur[31]*TMath::Power(Q,4) - 
      2*RMvcs[2][2]*Ur[46]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*RMvcs[1][0]*Ur[70]*TMath::Power(Q,3) + 
      2*RMvcs[1][1]*Ur[85]*TMath::Power(Q,4) - 
      2*RMvcs[1][2]*Ur[97]*TMath::Power(Q,5)) + 
   Jem[3][2]*(-2*RMvcs[3][0]*Ur[70]*TMath::Power(Q,3) + 
      2*RMvcs[3][1]*Ur[85]*TMath::Power(Q,4) - 
      2*RMvcs[3][2]*Ur[97]*TMath::Power(Q,5));
        SigmaIPol0[3]=Jem[0][1]*(-2*RMvcs[0][2]*Ur[42]*TMath::Power(Q,3) - 
      2*RMvcs[0][0]*Ur[5]*TMath::Power(Q,3) + 
      2*RMvcs[0][1]*Ur[27]*TMath::Power(Q,4)) + 
   Jem[2][1]*(-2*RMvcs[2][2]*Ur[42]*TMath::Power(Q,3) - 
      2*RMvcs[2][0]*Ur[5]*TMath::Power(Q,3) + 
      2*RMvcs[2][1]*Ur[27]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*RMvcs[0][0]*Ur[13]*TMath::Power(Q,3) + 
      2*RMvcs[0][1]*Ur[32]*TMath::Power(Q,4) - 
      2*RMvcs[0][2]*Ur[48]*TMath::Power(Q,5)) + 
   Jem[2][2]*(-2*RMvcs[2][0]*Ur[13]*TMath::Power(Q,3) + 
      2*RMvcs[2][1]*Ur[32]*TMath::Power(Q,4) - 
      2*RMvcs[2][2]*Ur[48]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*RMvcs[1][0]*Ur[71]*TMath::Power(Q,3) - 
      2*RMvcs[1][1]*Ur[85]*TMath::Power(Q,4) - 
      2*RMvcs[1][2]*Ur[98]*TMath::Power(Q,5)) + 
   Jem[3][2]*(-2*RMvcs[3][0]*Ur[71]*TMath::Power(Q,3) - 
      2*RMvcs[3][1]*Ur[85]*TMath::Power(Q,4) - 
      2*RMvcs[3][2]*Ur[98]*TMath::Power(Q,5));
        SigmaIPol0[4]=Jem[0][1]*(-2*RMvcs[0][2]*Ur[2]*TMath::Power(Q,2) + 
      2*RMvcs[0][1]*Ur[25]*TMath::Power(Q,3) - 
      2*RMvcs[0][0]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[2][1]*(-2*RMvcs[2][2]*Ur[2]*TMath::Power(Q,2) + 
      2*RMvcs[2][1]*Ur[25]*TMath::Power(Q,3) - 
      2*RMvcs[2][0]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*RMvcs[0][2]*Ur[10]*TMath::Power(Q,2) - 
      2*RMvcs[0][0]*Ur[9]*TMath::Power(Q,4) + 
      2*RMvcs[0][1]*Ur[30]*TMath::Power(Q,5)) + 
   Jem[2][2]*(-2*RMvcs[2][2]*Ur[10]*TMath::Power(Q,2) - 
      2*RMvcs[2][0]*Ur[9]*TMath::Power(Q,4) + 
      2*RMvcs[2][1]*Ur[30]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*RMvcs[1][2]*Ur[51]*TMath::Power(Q,2) + 
      2*RMvcs[1][0]*Ur[16]*TMath::Power(Q,4) + 
      2*RMvcs[1][1]*Ur[84]*TMath::Power(Q,5)) + 
   Jem[3][2]*(-2*RMvcs[3][2]*Ur[51]*TMath::Power(Q,2) + 
      2*RMvcs[3][0]*Ur[16]*TMath::Power(Q,4) + 
      2*RMvcs[3][1]*Ur[84]*TMath::Power(Q,5));
        SigmaIPol0[5]=Jem[0][1]*(-2*RMvcs[0][2]*Ur[41]*Q - 
      2*RMvcs[0][0]*Ur[4]*TMath::Power(Q,3)) + 
   Jem[2][1]*(-2*RMvcs[2][2]*Ur[41]*Q - 
      2*RMvcs[2][0]*Ur[4]*TMath::Power(Q,3)) + 
   Jem[0][2]*(-2*RMvcs[0][2]*Ur[47]*Q - 
      2*RMvcs[0][0]*Ur[12]*TMath::Power(Q,5)) + 
   Jem[2][2]*(-2*RMvcs[2][2]*Ur[47]*Q - 
      2*RMvcs[2][0]*Ur[12]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*RMvcs[1][2]*Ur[53]*Q + 
      2*RMvcs[1][0]*Ur[18]*TMath::Power(Q,5)) + 
   Jem[3][2]*(-2*RMvcs[3][2]*Ur[53]*Q + 
      2*RMvcs[3][0]*Ur[18]*TMath::Power(Q,5));
        SigmaIPol0[6]=Jem[0][1]*(-2*IMvcs[0][1]*Ur[28]*TMath::Power(Q,2) + 
      2*IMvcs[0][2]*Ur[44]*TMath::Power(Q,3) + 
      2*IMvcs[0][0]*Ur[7]*TMath::Power(Q,3)) + 
   Jem[2][1]*(-2*IMvcs[2][1]*Ur[28]*TMath::Power(Q,2) + 
      2*IMvcs[2][2]*Ur[44]*TMath::Power(Q,3) + 
      2*IMvcs[2][0]*Ur[7]*TMath::Power(Q,3)) + 
   Jem[0][2]*(2*IMvcs[0][0]*Ur[15]*TMath::Power(Q,3) - 
      2*IMvcs[0][1]*Ur[33]*TMath::Power(Q,4) + 
      2*IMvcs[0][2]*Ur[50]*TMath::Power(Q,5)) + 
   Jem[2][2]*(2*IMvcs[2][0]*Ur[15]*TMath::Power(Q,3) - 
      2*IMvcs[2][1]*Ur[33]*TMath::Power(Q,4) + 
      2*IMvcs[2][2]*Ur[50]*TMath::Power(Q,5)) + 
   Jem[1][2]*(2*IMvcs[1][0]*Ur[72]*TMath::Power(Q,3) + 
      2*IMvcs[1][1]*Ur[85]*TMath::Power(Q,4) + 
      2*IMvcs[1][2]*Ur[99]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*IMvcs[3][0]*Ur[72]*TMath::Power(Q,3) + 
      2*IMvcs[3][1]*Ur[85]*TMath::Power(Q,4) + 
      2*IMvcs[3][2]*Ur[99]*TMath::Power(Q,5));
        SigmaIPol0[7]=Jem[0][1]*(2*IMvcs[0][2]*Ur[43]*TMath::Power(Q,2) - 
      4*IMvcs[0][1]*Ur[25]*TMath::Power(Q,3) + 
      2*IMvcs[0][0]*Ur[6]*TMath::Power(Q,4)) + 
   Jem[2][1]*(2*IMvcs[2][2]*Ur[43]*TMath::Power(Q,2) - 
      4*IMvcs[2][1]*Ur[25]*TMath::Power(Q,3) + 
      2*IMvcs[2][0]*Ur[6]*TMath::Power(Q,4)) + 
   Jem[0][2]*(2*IMvcs[0][2]*Ur[49]*TMath::Power(Q,2) + 
      2*IMvcs[0][0]*Ur[14]*TMath::Power(Q,4) - 
      4*IMvcs[0][1]*Ur[30]*TMath::Power(Q,5)) + 
   Jem[2][2]*(2*IMvcs[2][2]*Ur[49]*TMath::Power(Q,2) + 
      2*IMvcs[2][0]*Ur[14]*TMath::Power(Q,4) - 
      4*IMvcs[2][1]*Ur[30]*TMath::Power(Q,5)) + 
   Jem[1][2]*(2*IMvcs[1][2]*Ur[56]*TMath::Power(Q,2) - 
      2*IMvcs[1][0]*Ur[21]*TMath::Power(Q,4) - 
      4*IMvcs[1][1]*Ur[84]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*IMvcs[3][2]*Ur[56]*TMath::Power(Q,2) - 
      2*IMvcs[3][0]*Ur[21]*TMath::Power(Q,4) - 
      4*IMvcs[3][1]*Ur[84]*TMath::Power(Q,5));
        SigmaIPolX[0]=Jem[0][1]*(2*IMvcs[1][0]*Ur[59]*TMath::Power(Q,3) + 
      2*IMvcs[1][2]*Ur[86]*TMath::Power(Q,3) - 
      2*IMvcs[1][1]*Ur[74]*TMath::Power(Q,4)) + 
   Jem[2][1]*(-2*IMvcs[3][0]*Ur[59]*TMath::Power(Q,3) - 
      2*IMvcs[3][2]*Ur[86]*TMath::Power(Q,3) + 
      2*IMvcs[3][1]*Ur[74]*TMath::Power(Q,4)) + 
   Jem[1][2]*(2*IMvcs[0][0]*Ur[17]*TMath::Power(Q,3) - 
      2*IMvcs[0][1]*Ur[35]*TMath::Power(Q,4) + 
      2*IMvcs[0][2]*Ur[52]*TMath::Power(Q,5)) + 
   Jem[3][2]*(-2*IMvcs[2][0]*Ur[17]*TMath::Power(Q,3) + 
      2*IMvcs[2][1]*Ur[35]*TMath::Power(Q,4) - 
      2*IMvcs[2][2]*Ur[52]*TMath::Power(Q,5)) + 
   Jem[0][2]*(2*IMvcs[1][0]*Ur[64]*TMath::Power(Q,3) - 
      2*IMvcs[1][1]*Ur[79]*TMath::Power(Q,4) + 
      2*IMvcs[1][2]*Ur[91]*TMath::Power(Q,5)) + 
   Jem[2][2]*(-2*IMvcs[3][0]*Ur[64]*TMath::Power(Q,3) + 
      2*IMvcs[3][1]*Ur[79]*TMath::Power(Q,4) - 
      2*IMvcs[3][2]*Ur[91]*TMath::Power(Q,5));
        SigmaIPolX[1]=Jem[0][1]*(2*IMvcs[1][0]*Ur[60]*TMath::Power(Q,3) + 
      2*IMvcs[1][2]*Ur[87]*TMath::Power(Q,3) - 
      2*IMvcs[1][1]*Ur[75]*TMath::Power(Q,4)) + 
   Jem[2][1]*(-2*IMvcs[3][0]*Ur[60]*TMath::Power(Q,3) - 
      2*IMvcs[3][2]*Ur[87]*TMath::Power(Q,3) + 
      2*IMvcs[3][1]*Ur[75]*TMath::Power(Q,4)) + 
   Jem[1][2]*(2*IMvcs[0][0]*Ur[19]*TMath::Power(Q,3) - 
      2*IMvcs[0][1]*Ur[36]*TMath::Power(Q,4) + 
      2*IMvcs[0][2]*Ur[54]*TMath::Power(Q,5)) + 
   Jem[3][2]*(-2*IMvcs[2][0]*Ur[19]*TMath::Power(Q,3) + 
      2*IMvcs[2][1]*Ur[36]*TMath::Power(Q,4) - 
      2*IMvcs[2][2]*Ur[54]*TMath::Power(Q,5)) + 
   Jem[0][2]*(2*IMvcs[1][0]*Ur[65]*TMath::Power(Q,3) - 
      2*IMvcs[1][1]*Ur[80]*TMath::Power(Q,4) + 
      2*IMvcs[1][2]*Ur[92]*TMath::Power(Q,5)) + 
   Jem[2][2]*(-2*IMvcs[3][0]*Ur[65]*TMath::Power(Q,3) + 
      2*IMvcs[3][1]*Ur[80]*TMath::Power(Q,4) - 
      2*IMvcs[3][2]*Ur[92]*TMath::Power(Q,5));
        SigmaIPolX[2]=Jem[0][1]*(2*IMvcs[1][2]*Ur[2]*TMath::Power(Q,2) - 
      2*IMvcs[1][1]*Ur[73]*TMath::Power(Q,3) - 
      2*IMvcs[1][0]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[2][1]*(-2*IMvcs[3][2]*Ur[2]*TMath::Power(Q,2) + 
      2*IMvcs[3][1]*Ur[73]*TMath::Power(Q,3) + 
      2*IMvcs[3][0]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[1][2]*(2*IMvcs[0][2]*Ur[51]*TMath::Power(Q,2) + 
      2*IMvcs[0][0]*Ur[16]*TMath::Power(Q,4) - 
      2*IMvcs[0][1]*Ur[34]*TMath::Power(Q,5)) + 
   Jem[3][2]*(-2*IMvcs[2][2]*Ur[51]*TMath::Power(Q,2) - 
      2*IMvcs[2][0]*Ur[16]*TMath::Power(Q,4) + 
      2*IMvcs[2][1]*Ur[34]*TMath::Power(Q,5)) + 
   Jem[0][2]*(2*IMvcs[1][2]*Ur[10]*TMath::Power(Q,2) - 
      2*IMvcs[1][0]*Ur[9]*TMath::Power(Q,4) - 
      2*IMvcs[1][1]*Ur[78]*TMath::Power(Q,5)) + 
   Jem[2][2]*(-2*IMvcs[3][2]*Ur[10]*TMath::Power(Q,2) + 
      2*IMvcs[3][0]*Ur[9]*TMath::Power(Q,4) + 
      2*IMvcs[3][1]*Ur[78]*TMath::Power(Q,5));
        SigmaIPolX[3]=Jem[0][1]*(2*IMvcs[1][2]*Ur[41]*Q - 
      2*IMvcs[1][0]*Ur[4]*TMath::Power(Q,3)) + 
   Jem[2][1]*(-2*IMvcs[3][2]*Ur[41]*Q + 
      2*IMvcs[3][0]*Ur[4]*TMath::Power(Q,3)) + 
   Jem[0][2]*(2*IMvcs[1][2]*Ur[47]*Q - 
      2*IMvcs[1][0]*Ur[12]*TMath::Power(Q,5)) + 
   Jem[2][2]*(-2*IMvcs[3][2]*Ur[47]*Q + 
      2*IMvcs[3][0]*Ur[12]*TMath::Power(Q,5)) + 
   Jem[1][2]*(2*IMvcs[0][2]*Ur[53]*Q + 
      2*IMvcs[0][0]*Ur[18]*TMath::Power(Q,5)) + 
   Jem[3][2]*(-2*IMvcs[2][2]*Ur[53]*Q - 
      2*IMvcs[2][0]*Ur[18]*TMath::Power(Q,5));
        SigmaIPolX[4]=Jem[1][2]*(-2*RMvcs[0][0]*Ur[20]*TMath::Power(Q,2) + 
      2*RMvcs[0][1]*Ur[37]*TMath::Power(Q,3) - 
      2*RMvcs[0][2]*Ur[55]*TMath::Power(Q,4)) + 
   Jem[3][2]*(2*RMvcs[2][0]*Ur[20]*TMath::Power(Q,2) - 
      2*RMvcs[2][1]*Ur[37]*TMath::Power(Q,3) + 
      2*RMvcs[2][2]*Ur[55]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*RMvcs[1][0]*Ur[61]*TMath::Power(Q,2) + 
      2*RMvcs[1][1]*Ur[76]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[88]*TMath::Power(Q,4)) + 
   Jem[2][1]*(2*RMvcs[3][0]*Ur[61]*TMath::Power(Q,2) - 
      2*RMvcs[3][1]*Ur[76]*TMath::Power(Q,3) + 
      2*RMvcs[3][2]*Ur[88]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*RMvcs[1][0]*Ur[66]*TMath::Power(Q,2) + 
      2*RMvcs[1][1]*Ur[81]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[93]*TMath::Power(Q,4)) + 
   Jem[2][2]*(2*RMvcs[3][0]*Ur[66]*TMath::Power(Q,2) - 
      2*RMvcs[3][1]*Ur[81]*TMath::Power(Q,3) + 
      2*RMvcs[3][2]*Ur[93]*TMath::Power(Q,4));
        SigmaIPolX[5]=Jem[1][2]*(-2*RMvcs[0][0]*Ur[22]*TMath::Power(Q,2) - 
      2*RMvcs[0][1]*Ur[37]*TMath::Power(Q,3) - 
      2*RMvcs[0][2]*Ur[57]*TMath::Power(Q,4)) + 
   Jem[3][2]*(2*RMvcs[2][0]*Ur[22]*TMath::Power(Q,2) + 
      2*RMvcs[2][1]*Ur[37]*TMath::Power(Q,3) + 
      2*RMvcs[2][2]*Ur[57]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*RMvcs[1][0]*Ur[62]*TMath::Power(Q,2) - 
      2*RMvcs[1][1]*Ur[76]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[89]*TMath::Power(Q,4)) + 
   Jem[2][1]*(2*RMvcs[3][0]*Ur[62]*TMath::Power(Q,2) + 
      2*RMvcs[3][1]*Ur[76]*TMath::Power(Q,3) + 
      2*RMvcs[3][2]*Ur[89]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*RMvcs[1][0]*Ur[67]*TMath::Power(Q,2) - 
      2*RMvcs[1][1]*Ur[81]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[94]*TMath::Power(Q,4)) + 
   Jem[2][2]*(2*RMvcs[3][0]*Ur[67]*TMath::Power(Q,2) + 
      2*RMvcs[3][1]*Ur[81]*TMath::Power(Q,3) + 
      2*RMvcs[3][2]*Ur[94]*TMath::Power(Q,4));
        SigmaIPolX[6]=Jem[0][1]*(2*RMvcs[1][1]*Ur[77]*TMath::Power(Q,2) - 
      2*RMvcs[1][0]*Ur[63]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[90]*TMath::Power(Q,3)) + 
   Jem[2][1]*(-2*RMvcs[3][1]*Ur[77]*TMath::Power(Q,2) + 
      2*RMvcs[3][0]*Ur[63]*TMath::Power(Q,3) + 
      2*RMvcs[3][2]*Ur[90]*TMath::Power(Q,3)) + 
   Jem[1][2]*(-2*RMvcs[0][0]*Ur[23]*TMath::Power(Q,3) + 
      2*RMvcs[0][1]*Ur[38]*TMath::Power(Q,4) - 
      2*RMvcs[0][2]*Ur[58]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*RMvcs[2][0]*Ur[23]*TMath::Power(Q,3) - 
      2*RMvcs[2][1]*Ur[38]*TMath::Power(Q,4) + 
      2*RMvcs[2][2]*Ur[58]*TMath::Power(Q,5)) + 
   Jem[0][2]*(-2*RMvcs[1][0]*Ur[68]*TMath::Power(Q,3) + 
      2*RMvcs[1][1]*Ur[82]*TMath::Power(Q,4) - 
      2*RMvcs[1][2]*Ur[95]*TMath::Power(Q,5)) + 
   Jem[2][2]*(2*RMvcs[3][0]*Ur[68]*TMath::Power(Q,3) - 
      2*RMvcs[3][1]*Ur[82]*TMath::Power(Q,4) + 
      2*RMvcs[3][2]*Ur[95]*TMath::Power(Q,5));
        SigmaIPolX[7]=Jem[0][1]*(-2*RMvcs[1][2]*Ur[43]*TMath::Power(Q,2) + 
      4*RMvcs[1][1]*Ur[73]*TMath::Power(Q,3) + 
      2*RMvcs[1][0]*Ur[6]*TMath::Power(Q,4)) + 
   Jem[2][1]*(2*RMvcs[3][2]*Ur[43]*TMath::Power(Q,2) - 
      4*RMvcs[3][1]*Ur[73]*TMath::Power(Q,3) - 
      2*RMvcs[3][0]*Ur[6]*TMath::Power(Q,4)) + 
   Jem[1][2]*(-2*RMvcs[0][2]*Ur[56]*TMath::Power(Q,2) - 
      2*RMvcs[0][0]*Ur[21]*TMath::Power(Q,4) + 
      4*RMvcs[0][1]*Ur[34]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*RMvcs[2][2]*Ur[56]*TMath::Power(Q,2) + 
      2*RMvcs[2][0]*Ur[21]*TMath::Power(Q,4) - 
      4*RMvcs[2][1]*Ur[34]*TMath::Power(Q,5)) + 
   Jem[0][2]*(-2*RMvcs[1][2]*Ur[49]*TMath::Power(Q,2) + 
      2*RMvcs[1][0]*Ur[14]*TMath::Power(Q,4) + 
      4*RMvcs[1][1]*Ur[78]*TMath::Power(Q,5)) + 
   Jem[2][2]*(2*RMvcs[3][2]*Ur[49]*TMath::Power(Q,2) - 
      2*RMvcs[3][0]*Ur[14]*TMath::Power(Q,4) - 
      4*RMvcs[3][1]*Ur[78]*TMath::Power(Q,5));
        SigmaIPolY[0]=Jem[2][1]*(2*IMvcs[0][0]*Ur[0]*TMath::Power(Q,2) - 
      2*IMvcs[0][1]*Ur[24]*TMath::Power(Q,3) + 
      2*IMvcs[0][2]*Ur[39]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*IMvcs[2][0]*Ur[0]*TMath::Power(Q,2) + 
      2*IMvcs[2][1]*Ur[24]*TMath::Power(Q,3) - 
      2*IMvcs[2][2]*Ur[39]*TMath::Power(Q,4)) + 
   Jem[2][2]*(2*IMvcs[0][0]*Ur[8]*TMath::Power(Q,2) - 
      2*IMvcs[0][1]*Ur[29]*TMath::Power(Q,3) + 
      2*IMvcs[0][2]*Ur[45]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*IMvcs[2][0]*Ur[8]*TMath::Power(Q,2) + 
      2*IMvcs[2][1]*Ur[29]*TMath::Power(Q,3) - 
      2*IMvcs[2][2]*Ur[45]*TMath::Power(Q,4)) + 
   Jem[3][2]*(2*IMvcs[1][0]*Ur[69]*TMath::Power(Q,2) - 
      2*IMvcs[1][1]*Ur[83]*TMath::Power(Q,3) + 
      2*IMvcs[1][2]*Ur[96]*TMath::Power(Q,4)) + 
   Jem[1][2]*(-2*IMvcs[3][0]*Ur[69]*TMath::Power(Q,2) + 
      2*IMvcs[3][1]*Ur[83]*TMath::Power(Q,3) - 
      2*IMvcs[3][2]*Ur[96]*TMath::Power(Q,4));
        SigmaIPolY[1]=Jem[2][1]*(2*IMvcs[0][0]*Ur[2]*TMath::Power(Q,2) + 
      2*IMvcs[0][1]*Ur[24]*TMath::Power(Q,3) + 
      2*IMvcs[0][2]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*IMvcs[2][0]*Ur[2]*TMath::Power(Q,2) - 
      2*IMvcs[2][1]*Ur[24]*TMath::Power(Q,3) - 
      2*IMvcs[2][2]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[3][2]*(2*IMvcs[1][0]*Ur[51]*TMath::Power(Q,2) + 
      2*IMvcs[1][1]*Ur[83]*TMath::Power(Q,3) - 
      2*IMvcs[1][2]*Ur[16]*TMath::Power(Q,4)) + 
   Jem[1][2]*(-2*IMvcs[3][0]*Ur[51]*TMath::Power(Q,2) - 
      2*IMvcs[3][1]*Ur[83]*TMath::Power(Q,3) + 
      2*IMvcs[3][2]*Ur[16]*TMath::Power(Q,4)) + 
   Jem[2][2]*(2*IMvcs[0][0]*Ur[10]*TMath::Power(Q,2) + 
      2*IMvcs[0][1]*Ur[29]*TMath::Power(Q,3) + 
      2*IMvcs[0][2]*Ur[9]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*IMvcs[2][0]*Ur[10]*TMath::Power(Q,2) - 
      2*IMvcs[2][1]*Ur[29]*TMath::Power(Q,3) - 
      2*IMvcs[2][2]*Ur[9]*TMath::Power(Q,4));
        SigmaIPolY[2]=Jem[2][1]*(2*IMvcs[0][0]*Ur[3]*TMath::Power(Q,3) + 
      2*IMvcs[0][2]*Ur[40]*TMath::Power(Q,3) - 
      2*IMvcs[0][1]*Ur[26]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*IMvcs[2][0]*Ur[3]*TMath::Power(Q,3) - 
      2*IMvcs[2][2]*Ur[40]*TMath::Power(Q,3) + 
      2*IMvcs[2][1]*Ur[26]*TMath::Power(Q,4)) + 
   Jem[2][2]*(2*IMvcs[0][0]*Ur[11]*TMath::Power(Q,3) - 
      2*IMvcs[0][1]*Ur[31]*TMath::Power(Q,4) + 
      2*IMvcs[0][2]*Ur[46]*TMath::Power(Q,5)) + 
   Jem[0][2]*(-2*IMvcs[2][0]*Ur[11]*TMath::Power(Q,3) + 
      2*IMvcs[2][1]*Ur[31]*TMath::Power(Q,4) - 
      2*IMvcs[2][2]*Ur[46]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*IMvcs[1][0]*Ur[70]*TMath::Power(Q,3) - 
      2*IMvcs[1][1]*Ur[85]*TMath::Power(Q,4) + 
      2*IMvcs[1][2]*Ur[97]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*IMvcs[3][0]*Ur[70]*TMath::Power(Q,3) + 
      2*IMvcs[3][1]*Ur[85]*TMath::Power(Q,4) - 
      2*IMvcs[3][2]*Ur[97]*TMath::Power(Q,5));
        SigmaIPolY[3]=Jem[2][1]*(2*IMvcs[0][2]*Ur[42]*TMath::Power(Q,3) + 
      2*IMvcs[0][0]*Ur[5]*TMath::Power(Q,3) - 
      2*IMvcs[0][1]*Ur[27]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*IMvcs[2][2]*Ur[42]*TMath::Power(Q,3) - 
      2*IMvcs[2][0]*Ur[5]*TMath::Power(Q,3) + 
      2*IMvcs[2][1]*Ur[27]*TMath::Power(Q,4)) + 
   Jem[2][2]*(2*IMvcs[0][0]*Ur[13]*TMath::Power(Q,3) - 
      2*IMvcs[0][1]*Ur[32]*TMath::Power(Q,4) + 
      2*IMvcs[0][2]*Ur[48]*TMath::Power(Q,5)) + 
   Jem[0][2]*(-2*IMvcs[2][0]*Ur[13]*TMath::Power(Q,3) + 
      2*IMvcs[2][1]*Ur[32]*TMath::Power(Q,4) - 
      2*IMvcs[2][2]*Ur[48]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*IMvcs[1][0]*Ur[71]*TMath::Power(Q,3) + 
      2*IMvcs[1][1]*Ur[85]*TMath::Power(Q,4) + 
      2*IMvcs[1][2]*Ur[98]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*IMvcs[3][0]*Ur[71]*TMath::Power(Q,3) - 
      2*IMvcs[3][1]*Ur[85]*TMath::Power(Q,4) - 
      2*IMvcs[3][2]*Ur[98]*TMath::Power(Q,5));
        SigmaIPolY[4]=Jem[2][1]*(2*IMvcs[0][2]*Ur[2]*TMath::Power(Q,2) - 
      2*IMvcs[0][1]*Ur[25]*TMath::Power(Q,3) + 
      2*IMvcs[0][0]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*IMvcs[2][2]*Ur[2]*TMath::Power(Q,2) + 
      2*IMvcs[2][1]*Ur[25]*TMath::Power(Q,3) - 
      2*IMvcs[2][0]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[2][2]*(2*IMvcs[0][2]*Ur[10]*TMath::Power(Q,2) + 
      2*IMvcs[0][0]*Ur[9]*TMath::Power(Q,4) - 
      2*IMvcs[0][1]*Ur[30]*TMath::Power(Q,5)) + 
   Jem[0][2]*(-2*IMvcs[2][2]*Ur[10]*TMath::Power(Q,2) - 
      2*IMvcs[2][0]*Ur[9]*TMath::Power(Q,4) + 
      2*IMvcs[2][1]*Ur[30]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*IMvcs[1][2]*Ur[51]*TMath::Power(Q,2) - 
      2*IMvcs[1][0]*Ur[16]*TMath::Power(Q,4) - 
      2*IMvcs[1][1]*Ur[84]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*IMvcs[3][2]*Ur[51]*TMath::Power(Q,2) + 
      2*IMvcs[3][0]*Ur[16]*TMath::Power(Q,4) + 
      2*IMvcs[3][1]*Ur[84]*TMath::Power(Q,5));
        SigmaIPolY[5]=Jem[2][1]*(2*IMvcs[0][2]*Ur[41]*Q + 
      2*IMvcs[0][0]*Ur[4]*TMath::Power(Q,3)) + 
   Jem[0][1]*(-2*IMvcs[2][2]*Ur[41]*Q - 
      2*IMvcs[2][0]*Ur[4]*TMath::Power(Q,3)) + 
   Jem[2][2]*(2*IMvcs[0][2]*Ur[47]*Q + 
      2*IMvcs[0][0]*Ur[12]*TMath::Power(Q,5)) + 
   Jem[0][2]*(-2*IMvcs[2][2]*Ur[47]*Q - 
      2*IMvcs[2][0]*Ur[12]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*IMvcs[1][2]*Ur[53]*Q - 
      2*IMvcs[1][0]*Ur[18]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*IMvcs[3][2]*Ur[53]*Q + 
      2*IMvcs[3][0]*Ur[18]*TMath::Power(Q,5));
        SigmaIPolY[6]=Jem[2][1]*(-2*RMvcs[0][1]*Ur[28]*TMath::Power(Q,2) + 
      2*RMvcs[0][2]*Ur[44]*TMath::Power(Q,3) + 
      2*RMvcs[0][0]*Ur[7]*TMath::Power(Q,3)) + 
   Jem[0][1]*(2*RMvcs[2][1]*Ur[28]*TMath::Power(Q,2) - 
      2*RMvcs[2][2]*Ur[44]*TMath::Power(Q,3) - 
      2*RMvcs[2][0]*Ur[7]*TMath::Power(Q,3)) + 
   Jem[2][2]*(2*RMvcs[0][0]*Ur[15]*TMath::Power(Q,3) - 
      2*RMvcs[0][1]*Ur[33]*TMath::Power(Q,4) + 
      2*RMvcs[0][2]*Ur[50]*TMath::Power(Q,5)) + 
   Jem[0][2]*(-2*RMvcs[2][0]*Ur[15]*TMath::Power(Q,3) + 
      2*RMvcs[2][1]*Ur[33]*TMath::Power(Q,4) - 
      2*RMvcs[2][2]*Ur[50]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*RMvcs[1][0]*Ur[72]*TMath::Power(Q,3) + 
      2*RMvcs[1][1]*Ur[85]*TMath::Power(Q,4) + 
      2*RMvcs[1][2]*Ur[99]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*RMvcs[3][0]*Ur[72]*TMath::Power(Q,3) - 
      2*RMvcs[3][1]*Ur[85]*TMath::Power(Q,4) - 
      2*RMvcs[3][2]*Ur[99]*TMath::Power(Q,5));
        SigmaIPolY[7]=Jem[2][1]*(2*RMvcs[0][2]*Ur[43]*TMath::Power(Q,2) - 
      4*RMvcs[0][1]*Ur[25]*TMath::Power(Q,3) + 
      2*RMvcs[0][0]*Ur[6]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*RMvcs[2][2]*Ur[43]*TMath::Power(Q,2) + 
      4*RMvcs[2][1]*Ur[25]*TMath::Power(Q,3) - 
      2*RMvcs[2][0]*Ur[6]*TMath::Power(Q,4)) + 
   Jem[2][2]*(2*RMvcs[0][2]*Ur[49]*TMath::Power(Q,2) + 
      2*RMvcs[0][0]*Ur[14]*TMath::Power(Q,4) - 
      4*RMvcs[0][1]*Ur[30]*TMath::Power(Q,5)) + 
   Jem[0][2]*(-2*RMvcs[2][2]*Ur[49]*TMath::Power(Q,2) - 
      2*RMvcs[2][0]*Ur[14]*TMath::Power(Q,4) + 
      4*RMvcs[2][1]*Ur[30]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*RMvcs[1][2]*Ur[56]*TMath::Power(Q,2) - 
      2*RMvcs[1][0]*Ur[21]*TMath::Power(Q,4) - 
      4*RMvcs[1][1]*Ur[84]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*RMvcs[3][2]*Ur[56]*TMath::Power(Q,2) + 
      2*RMvcs[3][0]*Ur[21]*TMath::Power(Q,4) + 
      4*RMvcs[3][1]*Ur[84]*TMath::Power(Q,5));
        SigmaIPolZ[0]=Jem[2][1]*(2*IMvcs[1][0]*Ur[59]*TMath::Power(Q,3) + 
      2*IMvcs[1][2]*Ur[86]*TMath::Power(Q,3) - 
      2*IMvcs[1][1]*Ur[74]*TMath::Power(Q,4)) + 
   Jem[0][1]*(2*IMvcs[3][0]*Ur[59]*TMath::Power(Q,3) + 
      2*IMvcs[3][2]*Ur[86]*TMath::Power(Q,3) - 
      2*IMvcs[3][1]*Ur[74]*TMath::Power(Q,4)) + 
   Jem[3][2]*(2*IMvcs[0][0]*Ur[17]*TMath::Power(Q,3) - 
      2*IMvcs[0][1]*Ur[35]*TMath::Power(Q,4) + 
      2*IMvcs[0][2]*Ur[52]*TMath::Power(Q,5)) + 
   Jem[1][2]*(2*IMvcs[2][0]*Ur[17]*TMath::Power(Q,3) - 
      2*IMvcs[2][1]*Ur[35]*TMath::Power(Q,4) + 
      2*IMvcs[2][2]*Ur[52]*TMath::Power(Q,5)) + 
   Jem[2][2]*(2*IMvcs[1][0]*Ur[64]*TMath::Power(Q,3) - 
      2*IMvcs[1][1]*Ur[79]*TMath::Power(Q,4) + 
      2*IMvcs[1][2]*Ur[91]*TMath::Power(Q,5)) + 
   Jem[0][2]*(2*IMvcs[3][0]*Ur[64]*TMath::Power(Q,3) - 
      2*IMvcs[3][1]*Ur[79]*TMath::Power(Q,4) + 
      2*IMvcs[3][2]*Ur[91]*TMath::Power(Q,5));
        SigmaIPolZ[1]=Jem[2][1]*(2*IMvcs[1][0]*Ur[60]*TMath::Power(Q,3) + 
      2*IMvcs[1][2]*Ur[87]*TMath::Power(Q,3) - 
      2*IMvcs[1][1]*Ur[75]*TMath::Power(Q,4)) + 
   Jem[0][1]*(2*IMvcs[3][0]*Ur[60]*TMath::Power(Q,3) + 
      2*IMvcs[3][2]*Ur[87]*TMath::Power(Q,3) - 
      2*IMvcs[3][1]*Ur[75]*TMath::Power(Q,4)) + 
   Jem[3][2]*(2*IMvcs[0][0]*Ur[19]*TMath::Power(Q,3) - 
      2*IMvcs[0][1]*Ur[36]*TMath::Power(Q,4) + 
      2*IMvcs[0][2]*Ur[54]*TMath::Power(Q,5)) + 
   Jem[1][2]*(2*IMvcs[2][0]*Ur[19]*TMath::Power(Q,3) - 
      2*IMvcs[2][1]*Ur[36]*TMath::Power(Q,4) + 
      2*IMvcs[2][2]*Ur[54]*TMath::Power(Q,5)) + 
   Jem[2][2]*(2*IMvcs[1][0]*Ur[65]*TMath::Power(Q,3) - 
      2*IMvcs[1][1]*Ur[80]*TMath::Power(Q,4) + 
      2*IMvcs[1][2]*Ur[92]*TMath::Power(Q,5)) + 
   Jem[0][2]*(2*IMvcs[3][0]*Ur[65]*TMath::Power(Q,3) - 
      2*IMvcs[3][1]*Ur[80]*TMath::Power(Q,4) + 
      2*IMvcs[3][2]*Ur[92]*TMath::Power(Q,5));
        SigmaIPolZ[2]=Jem[2][1]*(2*IMvcs[1][2]*Ur[2]*TMath::Power(Q,2) - 
      2*IMvcs[1][1]*Ur[73]*TMath::Power(Q,3) - 
      2*IMvcs[1][0]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[0][1]*(2*IMvcs[3][2]*Ur[2]*TMath::Power(Q,2) - 
      2*IMvcs[3][1]*Ur[73]*TMath::Power(Q,3) - 
      2*IMvcs[3][0]*Ur[1]*TMath::Power(Q,4)) + 
   Jem[3][2]*(2*IMvcs[0][2]*Ur[51]*TMath::Power(Q,2) + 
      2*IMvcs[0][0]*Ur[16]*TMath::Power(Q,4) - 
      2*IMvcs[0][1]*Ur[34]*TMath::Power(Q,5)) + 
   Jem[1][2]*(2*IMvcs[2][2]*Ur[51]*TMath::Power(Q,2) + 
      2*IMvcs[2][0]*Ur[16]*TMath::Power(Q,4) - 
      2*IMvcs[2][1]*Ur[34]*TMath::Power(Q,5)) + 
   Jem[2][2]*(2*IMvcs[1][2]*Ur[10]*TMath::Power(Q,2) - 
      2*IMvcs[1][0]*Ur[9]*TMath::Power(Q,4) - 
      2*IMvcs[1][1]*Ur[78]*TMath::Power(Q,5)) + 
   Jem[0][2]*(2*IMvcs[3][2]*Ur[10]*TMath::Power(Q,2) - 
      2*IMvcs[3][0]*Ur[9]*TMath::Power(Q,4) - 
      2*IMvcs[3][1]*Ur[78]*TMath::Power(Q,5));
        SigmaIPolZ[3]=Jem[2][1]*(2*IMvcs[1][2]*Ur[41]*Q - 
      2*IMvcs[1][0]*Ur[4]*TMath::Power(Q,3)) + 
   Jem[0][1]*(2*IMvcs[3][2]*Ur[41]*Q - 
      2*IMvcs[3][0]*Ur[4]*TMath::Power(Q,3)) + 
   Jem[2][2]*(2*IMvcs[1][2]*Ur[47]*Q - 
      2*IMvcs[1][0]*Ur[12]*TMath::Power(Q,5)) + 
   Jem[0][2]*(2*IMvcs[3][2]*Ur[47]*Q - 
      2*IMvcs[3][0]*Ur[12]*TMath::Power(Q,5)) + 
   Jem[3][2]*(2*IMvcs[0][2]*Ur[53]*Q + 
      2*IMvcs[0][0]*Ur[18]*TMath::Power(Q,5)) + 
   Jem[1][2]*(2*IMvcs[2][2]*Ur[53]*Q + 
      2*IMvcs[2][0]*Ur[18]*TMath::Power(Q,5));
        SigmaIPolZ[4]=Jem[3][2]*(-2*RMvcs[0][0]*Ur[20]*TMath::Power(Q,2) + 
      2*RMvcs[0][1]*Ur[37]*TMath::Power(Q,3) - 
      2*RMvcs[0][2]*Ur[55]*TMath::Power(Q,4)) + 
   Jem[1][2]*(-2*RMvcs[2][0]*Ur[20]*TMath::Power(Q,2) + 
      2*RMvcs[2][1]*Ur[37]*TMath::Power(Q,3) - 
      2*RMvcs[2][2]*Ur[55]*TMath::Power(Q,4)) + 
   Jem[2][1]*(-2*RMvcs[1][0]*Ur[61]*TMath::Power(Q,2) + 
      2*RMvcs[1][1]*Ur[76]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[88]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*RMvcs[3][0]*Ur[61]*TMath::Power(Q,2) + 
      2*RMvcs[3][1]*Ur[76]*TMath::Power(Q,3) - 
      2*RMvcs[3][2]*Ur[88]*TMath::Power(Q,4)) + 
   Jem[2][2]*(-2*RMvcs[1][0]*Ur[66]*TMath::Power(Q,2) + 
      2*RMvcs[1][1]*Ur[81]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[93]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*RMvcs[3][0]*Ur[66]*TMath::Power(Q,2) + 
      2*RMvcs[3][1]*Ur[81]*TMath::Power(Q,3) - 
      2*RMvcs[3][2]*Ur[93]*TMath::Power(Q,4));
        SigmaIPolZ[5]=Jem[3][2]*(-2*RMvcs[0][0]*Ur[22]*TMath::Power(Q,2) - 
      2*RMvcs[0][1]*Ur[37]*TMath::Power(Q,3) - 
      2*RMvcs[0][2]*Ur[57]*TMath::Power(Q,4)) + 
   Jem[1][2]*(-2*RMvcs[2][0]*Ur[22]*TMath::Power(Q,2) - 
      2*RMvcs[2][1]*Ur[37]*TMath::Power(Q,3) - 
      2*RMvcs[2][2]*Ur[57]*TMath::Power(Q,4)) + 
   Jem[2][1]*(-2*RMvcs[1][0]*Ur[62]*TMath::Power(Q,2) - 
      2*RMvcs[1][1]*Ur[76]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[89]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*RMvcs[3][0]*Ur[62]*TMath::Power(Q,2) - 
      2*RMvcs[3][1]*Ur[76]*TMath::Power(Q,3) - 
      2*RMvcs[3][2]*Ur[89]*TMath::Power(Q,4)) + 
   Jem[2][2]*(-2*RMvcs[1][0]*Ur[67]*TMath::Power(Q,2) - 
      2*RMvcs[1][1]*Ur[81]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[94]*TMath::Power(Q,4)) + 
   Jem[0][2]*(-2*RMvcs[3][0]*Ur[67]*TMath::Power(Q,2) - 
      2*RMvcs[3][1]*Ur[81]*TMath::Power(Q,3) - 
      2*RMvcs[3][2]*Ur[94]*TMath::Power(Q,4));
        SigmaIPolZ[6]=Jem[2][1]*(2*RMvcs[1][1]*Ur[77]*TMath::Power(Q,2) - 
      2*RMvcs[1][0]*Ur[63]*TMath::Power(Q,3) - 
      2*RMvcs[1][2]*Ur[90]*TMath::Power(Q,3)) + 
   Jem[0][1]*(2*RMvcs[3][1]*Ur[77]*TMath::Power(Q,2) - 
      2*RMvcs[3][0]*Ur[63]*TMath::Power(Q,3) - 
      2*RMvcs[3][2]*Ur[90]*TMath::Power(Q,3)) + 
   Jem[3][2]*(-2*RMvcs[0][0]*Ur[23]*TMath::Power(Q,3) + 
      2*RMvcs[0][1]*Ur[38]*TMath::Power(Q,4) - 
      2*RMvcs[0][2]*Ur[58]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*RMvcs[2][0]*Ur[23]*TMath::Power(Q,3) + 
      2*RMvcs[2][1]*Ur[38]*TMath::Power(Q,4) - 
      2*RMvcs[2][2]*Ur[58]*TMath::Power(Q,5)) + 
   Jem[2][2]*(-2*RMvcs[1][0]*Ur[68]*TMath::Power(Q,3) + 
      2*RMvcs[1][1]*Ur[82]*TMath::Power(Q,4) - 
      2*RMvcs[1][2]*Ur[95]*TMath::Power(Q,5)) + 
   Jem[0][2]*(-2*RMvcs[3][0]*Ur[68]*TMath::Power(Q,3) + 
      2*RMvcs[3][1]*Ur[82]*TMath::Power(Q,4) - 
      2*RMvcs[3][2]*Ur[95]*TMath::Power(Q,5));
        SigmaIPolZ[7]=Jem[2][1]*(-2*RMvcs[1][2]*Ur[43]*TMath::Power(Q,2) + 
      4*RMvcs[1][1]*Ur[73]*TMath::Power(Q,3) + 
      2*RMvcs[1][0]*Ur[6]*TMath::Power(Q,4)) + 
   Jem[0][1]*(-2*RMvcs[3][2]*Ur[43]*TMath::Power(Q,2) + 
      4*RMvcs[3][1]*Ur[73]*TMath::Power(Q,3) + 
      2*RMvcs[3][0]*Ur[6]*TMath::Power(Q,4)) + 
   Jem[3][2]*(-2*RMvcs[0][2]*Ur[56]*TMath::Power(Q,2) - 
      2*RMvcs[0][0]*Ur[21]*TMath::Power(Q,4) + 
      4*RMvcs[0][1]*Ur[34]*TMath::Power(Q,5)) + 
   Jem[1][2]*(-2*RMvcs[2][2]*Ur[56]*TMath::Power(Q,2) - 
      2*RMvcs[2][0]*Ur[21]*TMath::Power(Q,4) + 
      4*RMvcs[2][1]*Ur[34]*TMath::Power(Q,5)) + 
   Jem[2][2]*(-2*RMvcs[1][2]*Ur[49]*TMath::Power(Q,2) + 
      2*RMvcs[1][0]*Ur[14]*TMath::Power(Q,4) + 
      4*RMvcs[1][1]*Ur[78]*TMath::Power(Q,5)) + 
   Jem[0][2]*(-2*RMvcs[3][2]*Ur[49]*TMath::Power(Q,2) + 
      2*RMvcs[3][0]*Ur[14]*TMath::Power(Q,4) + 
      4*RMvcs[3][1]*Ur[78]*TMath::Power(Q,5));
    
    
    	// Flag
    	
    	InitLeadingVCSAndInterfCrossSections = kTRUE;    	
    	
    	
    	// Validation mode only
    
   	 	if ( Validation == kTRUE )
   	 	{
   	 		ofstream outfile;   	 		
   	 			
    		outfile.open("CheckUrLeading.dat");
    		for(Int_t i=0;i<100;i++)
			{
				outfile << Ur[i] << endl;
			} // end for i		
			outfile.close();	
		
			
    		outfile.open("CheckSigmaVCSLeading.dat");
			outfile << SigmaVCSPol0[0] << endl;
			outfile << SigmaVCSPolZ[0] << endl;
			outfile << SigmaVCSPol0[1] << endl;
			outfile << SigmaVCSPolX[0] << endl;
			outfile << SigmaVCSPolX[1] << endl;
			outfile << SigmaVCSPolZ[1] << endl;
			outfile << SigmaVCSPol0[2] << endl;
			outfile << SigmaVCSPol0[3] << endl;
			outfile << SigmaVCSPolY[0] << endl;
			outfile << SigmaVCSPolY[1] << endl;
			outfile << SigmaVCSPolY[2] << endl;
			outfile << SigmaVCSPolY[3] << endl;
			outfile << SigmaVCSPolY[4] << endl;
			outfile << SigmaVCSPol0[4] << endl;
			outfile << SigmaVCSPolX[2] << endl;
			outfile << SigmaVCSPolZ[2] << endl;
			outfile << SigmaVCSPolX[3] << endl;
			outfile << SigmaVCSPolZ[3] << endl;		
			outfile.close();	
		
			
    		outfile.open("CheckSigmaILeading.dat");
			for(Int_t i=0;i<32;i++)
			{
				if (i<8)
				{
					outfile << SigmaIPol0[i] << endl;
				} // end if i
				if ( 8 <= i && i < 16 )
				{
					outfile << SigmaIPolX[i-8] << endl;
				} // end if i
				if ( 16 <= i && i < 24 )
				{
					outfile << SigmaIPolY[i-16] << endl;
				} // end if i
				if ( 24 <= i && i < 32)
				{
					outfile << SigmaIPolZ[i-24] << endl;
				} // end if i
			} // end for i    	
		outfile.close();		
				
		} // end if Validation  
		
//	} // end if InitLeadingVCSAndInterfCrossSections
	 
} // end MakeLeadingVCSAndInterfCrossSections 
 
 
 
/*---------------------------- Function DdirectDcrossed(phi) ---------------------------*
 | Computes the denominator of the Bethe Heitler cross section.                         |
 | Only depends on kinematics.                                                          |
 *--------------------------------------------------------------------------------------*/

Double_t TGVKelly::DdirectDcrossed( Double_t phi )
{
	Double_t DDDC;
		
	DDDC=TMath::Power(TMath::Power(Q,2) + t,2)/4. - TMath::Power((TMath::Power(Q,2)*(TMath::Power(Q,2) + t*(-1 + 2*xB))*
        TMath::CosH(Omega))/(2.*TMath::Sqrt(TMath::Power(Q,4) + 4*TMath::Power(M,2)*TMath::Power(Q,2)*TMath::Power(xB,2)))
       - Q*qpPerp*TMath::Cos(phi)*TMath::SinH(Omega),2);

	return DDDC;
} // end DdirectDcrossed



/*--- Function SqrAmplBH(Q2Input,xBInput,tInput,phi,BeamHeli,TargetPolar,BeamCharge) ---*
 | Computes the expansion of the Bethe Heitler cross section (divided by phase space)   |
 | in terms of the "expansion cross sections" computed in the function                  | 
 | MakeExactCrossSections (or MakeLeadingCrossSections).                                |
 *--------------------------------------------------------------------------------------*/
 
 Double_t TGVKelly::SqrAmplBH( Double_t Q2Input, Double_t xBInput, Double_t tInput, Double_t phi, Int_t BeamHeli, Int_t TargetPolar, Bool_t Exact ) 
 { 	
 	TGVKelly::MakeKinematics( Q2Input, xBInput, tInput );
 	if ( Exact == kTRUE )
    {
    	TGVKelly::MakeExactBHCrossSections();
    } // end if MaxPrecision
	else
    { 	
    	TGVKelly::MakeLeadingBHCrossSections();
    } // end else if MaxPrecision 	
    
 	
 	if ( NoPrint == kFALSE )
 	{
 		cout << "TGVKelly : Computation of Bethe Heitler cross section" << endl;
		cout << "    TGVKelly : Kinematics initialization : " << InitKinematics << endl;
    	cout << "    TGVKelly : Exact Bethe Heitler cross sections initialization : " << InitExactBHCrossSections << endl;
    	cout << "    TGVKelly : Leading Bethe Heitler cross sections initialization : " << InitLeadingBHCrossSections << endl;
	} // end if NoPrint
 	
 	
 	if ( InitKinematics == kTRUE &&  ( InitExactBHCrossSections == kTRUE || InitLeadingBHCrossSections == kTRUE ) )
	{
 		Double_t DDDC=DdirectDcrossed( phi );
 		Double_t SigmaBH;
 	
 		Int_t he = BeamHeli;
 
 	
		if( TargetPolar == 0 )
		{
			SigmaBH=-(SigmaBHPol0[0] - SigmaBHPol0[3]*TMath::Cos(2*phi)*(-1 + TMath::CosH(2*Omega)) + 
      SigmaBHPol0[1]*TMath::CosH(2*Omega) + SigmaBHPol0[2]*TMath::Cos(phi)*TMath::SinH(2*Omega)
      )/(4.*DDDC*TMath::Power(t,2));
		} // end if TargetPolar
	
		if( TargetPolar == 1 )
		{
			SigmaBH=-(SigmaBHPolX[0]*he*TMath::CosH(Omega) + SigmaBHPolX[1]*he*TMath::Cos(phi)*TMath::SinH(Omega))/
   (4.*DDDC*TMath::Power(t,2));
 		} // end if TargetPolar
 		
		if( TargetPolar == 2 )
		{
			SigmaBH=-(he*SigmaBHPolY*TMath::Sin(phi)*TMath::SinH(Omega))/(4.*DDDC*TMath::Power(t,2));
		} // end if TargetPolar
	
		if( TargetPolar == 3 )
		{
			SigmaBH=-(SigmaBHPolZ[0]*he*TMath::CosH(Omega) + SigmaBHPolZ[1]*he*TMath::Cos(phi)*TMath::SinH(Omega))/
   (4.*DDDC*TMath::Power(t,2));
		} // end if TargetPolar	
	
		return SigmaBH;
		
	} // if Init (all variables)
	else
	{
		cout << "TGVKelly : Incomplete initialization !" << endl;
		exit(-1);
	} // end if Init (all variables)
	
} // end SqrAmplBH



/*-- Function SqrAmplVCS(Q2Input,xBInput,tInput,phi,BeamHeli,TargetPolar,BeamCharge) ---*
 | Computes the expansion of the Virtual Compton Scattering cross section (divided by   |
 | phase space) in terms of the "expansion cross sections" computed in the function     | 
 | MakeExactCrossSections (or MakeLeadingCrossSections).                                |
 *--------------------------------------------------------------------------------------*/
 
 Double_t TGVKelly::SqrAmplVCS( Double_t Q2Input, Double_t xBInput, Double_t tInput, Double_t phi, Int_t BeamHeli, Int_t TargetPolar, Double_t ReH, Double_t ImH, Double_t ReE, Double_t ImE, Double_t ReHT, Double_t ImHT, Double_t ReET, Double_t ImET, Bool_t Exact ) 
 { 	
 	TGVKelly::MakeKinematics( Q2Input, xBInput, tInput );
 	TGVKelly::MakeVCSHelicityAmplitudes( ReH, ImH, ReE, ImE, ReHT, ImHT, ReET, ImET );
 	if ( Exact == kTRUE )
    {
    	TGVKelly::MakeExactVCSAndInterfCrossSections();
    } // end if MaxPrecision
	else
    { 	
    	TGVKelly::MakeLeadingVCSAndInterfCrossSections();
    } // end else if MaxPrecision 	
    
 	
 	if ( NoPrint == kFALSE )
 	{
 		cout << "TGVKelly : Computation of Virtual Compton Scattering cross section" << endl;
		cout << "    TGVKelly : Kinematics initialization : " << InitKinematics << endl;
    	cout << "    TGVKelly : Helicity amplitude initialization : " << InitMvcs << endl;
    	cout << "    TGVKelly : Exact Virtual Compton Scattering and Interference cross sections initialization : " << InitExactVCSAndInterfCrossSections << endl;
    	cout << "    TGVKelly : Leading Virtual Compton Scattering and Interference cross sections initialization : " << InitLeadingVCSAndInterfCrossSections << endl;
	} // end if NoPrint
 	
 	
 	if ( InitKinematics == kTRUE && InitMvcs == kTRUE && ( InitExactVCSAndInterfCrossSections == kTRUE || InitLeadingVCSAndInterfCrossSections == kTRUE ) )
	{
 		Double_t SigmaVCS;
 	
 		Int_t he = BeamHeli;
 
 	
		if( TargetPolar == 0 )
		{
			SigmaVCS=(SigmaVCSPol0[0] - SigmaVCSPol0[3]*TMath::Cos(2*phi)*(-1 + TMath::CosH(2*Omega)) + 
     SigmaVCSPol0[1]*TMath::CosH(2*Omega) + 
     SigmaVCSPol0[4]*he*TMath::Sin(phi)*TMath::SinH(Omega) + 
     SigmaVCSPol0[2]*TMath::Cos(phi)*TMath::SinH(2*Omega))/(2.*TMath::Power(Q,2));
		} // end if TargetPolar
	
		if( TargetPolar == 1 )
		{
			SigmaVCS=(SigmaVCSPolX[0]*he*TMath::CosH(Omega) - 
     SigmaVCSPolX[3]*(-1 + TMath::CosH(2*Omega))*TMath::Sin(2*phi) + 
     SigmaVCSPolX[1]*he*TMath::Cos(phi)*TMath::SinH(Omega) + 
     SigmaVCSPolX[2]*TMath::Sin(phi)*TMath::SinH(2*Omega))/(2.*TMath::Power(Q,2));
 		} // end if TargetPolar
	
		if( TargetPolar == 2 )
		{
			SigmaVCS=(SigmaVCSPolY[1] - SigmaVCSPolY[4]*TMath::Cos(2*phi)*(-1 + TMath::CosH(2*Omega)) + 
     SigmaVCSPolY[2]*TMath::CosH(2*Omega) + 
     SigmaVCSPolY[0]*he*TMath::Sin(phi)*TMath::SinH(Omega) + 
     SigmaVCSPolY[3]*TMath::Cos(phi)*TMath::SinH(2*Omega))/(2.*TMath::Power(Q,2));
		} // end if TargetPolar
	
		if( TargetPolar == 3 )
		{
			SigmaVCS=(SigmaVCSPolZ[0]*he*TMath::CosH(Omega) - 
     SigmaVCSPolZ[3]*(-1 + TMath::CosH(2*Omega))*TMath::Sin(2*phi) + 
     SigmaVCSPolZ[1]*he*TMath::Cos(phi)*TMath::SinH(Omega) + 
     SigmaVCSPolZ[2]*TMath::Sin(phi)*TMath::SinH(2*Omega))/(2.*TMath::Power(Q,2));
		} // end if TargetPolar	
	
		return SigmaVCS;
		
	} // if Init (all variables)
	else
	{
		cout << "TGVKelly : Incomplete initialization !" << endl;
		exit(-1);
	} // end if Init (all variables)
	
} // end SqrAmplVCS



/*- Fonction SqrAmplInterf(Q2Input,xBInput,tInput,phi,BeamHeli,TargetPolar,BeamCharge)-*
 | Computes the expansion of the Interference cross section (divided by phase space)   |
 | in terms of the "expansion cross sections" computed in the function                 |
 | MakeExactCrossSections (or MakeLeadingCrossSections).                               |
 *-------------------------------------------------------------------------------------*/
  
 Double_t TGVKelly::SqrAmplInterf( Double_t Q2Input, Double_t xBInput, Double_t tInput, Double_t phi, Int_t BeamHeli, Int_t TargetPolar, Int_t BeamCharge, Double_t ReH, Double_t ImH, Double_t ReE, Double_t ImE, Double_t ReHT, Double_t ImHT, Double_t ReET, Double_t ImET, Bool_t Exact ) 
 { 	
 	TGVKelly::MakeKinematics( Q2Input, xBInput, tInput );
 	TGVKelly::MakeVCSHelicityAmplitudes( ReH, ImH, ReE, ImE, ReHT, ImHT, ReET, ImET );
 	if ( Exact == kTRUE )
    {
    	TGVKelly::MakeExactVCSAndInterfCrossSections();
    } // end if MaxPrecision
	else
    { 	
    	TGVKelly::MakeLeadingVCSAndInterfCrossSections();
    } // end else if MaxPrecision 	
    
 	
 	if ( NoPrint == kFALSE )
 	{
 		cout << "TGVKelly : Computation of Virtual Compton Scattering cross section" << endl;
		cout << "    TGVKelly : Kinematics initialization : " << InitKinematics << endl;
    	cout << "    TGVKelly : Helicity amplitude initialization : " << InitMvcs << endl;
    	cout << "    TGVKelly : Exact Virtual Compton Scattering and Interference cross sections initialization : " << InitExactVCSAndInterfCrossSections << endl;
    	cout << "    TGVKelly : Leading Virtual Compton Scattering and Interference cross sections initialization : " << InitLeadingVCSAndInterfCrossSections << endl;
	} // end if NoPrint
 	
 	
 	if ( InitKinematics == kTRUE && InitMvcs == kTRUE && ( InitExactVCSAndInterfCrossSections == kTRUE || InitLeadingVCSAndInterfCrossSections == kTRUE ) )
	{
 		Double_t DDDC = DdirectDcrossed( phi );
 		Double_t SigmaI;
 	
 		Int_t he = BeamHeli;
 
 	
		if( TargetPolar == 0 )
		{
			SigmaI=-((BeamCharge*(SigmaIPol0[0]*TMath::CosH(Omega) + 
         SigmaIPol0[4]*TMath::Cos(2*phi)*(TMath::CosH(Omega) - TMath::CosH(3*Omega)) + 
         SigmaIPol0[1]*TMath::CosH(3*Omega) - 
         SigmaIPol0[7]*he*(-1 + TMath::CosH(2*Omega))*TMath::Sin(2*phi) + 
         SigmaIPol0[6]*he*TMath::Sin(phi)*TMath::SinH(2*Omega) + 
         (SigmaIPol0[5]*TMath::Cos(3*phi)*(3*TMath::SinH(Omega) - TMath::SinH(3*Omega)))/3. + 
         TMath::Cos(phi)*(SigmaIPol0[2]*TMath::SinH(Omega) + 
            SigmaIPol0[3]*TMath::SinH(3*Omega))))/(DDDC*TMath::Power(Q,2)*t));
		} // end if TargetPolar
	
		if( TargetPolar == 1 )
		{
			SigmaI=-((BeamCharge*(-(SigmaIPolX[7]*he*TMath::Cos(2*phi)*(-1 + TMath::CosH(2*Omega))) + 
         he*(SigmaIPolX[4] + SigmaIPolX[5]*TMath::CosH(2*Omega)) + 
         SigmaIPolX[2]*(TMath::CosH(Omega) - TMath::CosH(3*Omega))*TMath::Sin(2*phi) + 
         SigmaIPolX[6]*he*TMath::Cos(phi)*TMath::SinH(2*Omega) + 
         (SigmaIPolX[3]*TMath::Sin(3*phi)*(3*TMath::SinH(Omega) - TMath::SinH(3*Omega)))/3. + 
         TMath::Sin(phi)*(SigmaIPolX[0]*TMath::SinH(Omega) + 
            SigmaIPolX[1]*TMath::SinH(3*Omega))))/(DDDC*TMath::Power(Q,2)*t));
 		} // end if TargetPolar
	
		if( TargetPolar == 2 )
		{
			SigmaI=-((BeamCharge*(SigmaIPolY[0]*TMath::CosH(Omega) + 
         SigmaIPolY[4]*TMath::Cos(2*phi)*(TMath::CosH(Omega) - TMath::CosH(3*Omega)) + 
         SigmaIPolY[1]*TMath::CosH(3*Omega) - 
         SigmaIPolY[7]*he*(-1 + TMath::CosH(2*Omega))*TMath::Sin(2*phi) + 
         SigmaIPolY[6]*he*TMath::Sin(phi)*TMath::SinH(2*Omega) + 
         (SigmaIPolY[5]*TMath::Cos(3*phi)*(3*TMath::SinH(Omega) - TMath::SinH(3*Omega)))/3. + 
         TMath::Cos(phi)*(SigmaIPolY[2]*TMath::SinH(Omega) + 
            SigmaIPolY[3]*TMath::SinH(3*Omega))))/(DDDC*TMath::Power(Q,2)*t));
		} // end if TargetPolar
	
		if( TargetPolar == 3 )
		{
			SigmaI=-((BeamCharge*(-(SigmaIPolZ[7]*he*TMath::Cos(2*phi)*(-1 + TMath::CosH(2*Omega))) + 
         he*(SigmaIPolZ[4] + SigmaIPolZ[5]*TMath::CosH(2*Omega)) + 
         SigmaIPolZ[2]*(TMath::CosH(Omega) - TMath::CosH(3*Omega))*TMath::Sin(2*phi) + 
         SigmaIPolZ[6]*he*TMath::Cos(phi)*TMath::SinH(2*Omega) + 
         (SigmaIPolZ[3]*TMath::Sin(3*phi)*(3*TMath::SinH(Omega) - TMath::SinH(3*Omega)))/3. + 
         TMath::Sin(phi)*(SigmaIPolZ[0]*TMath::SinH(Omega) + 
            SigmaIPolZ[1]*TMath::SinH(3*Omega))))/(DDDC*TMath::Power(Q,2)*t));
		} // end if TargetPolar	
	
		return SigmaI;
		
	} // if Init (all variables)
	else
	{
		cout << "TGVKelly : Incomplete initialization !" << endl;
		exit(-1);
	} // end if Init (all variables)
	
} // end SqrAmplInterf 	
 
 	

/*-Function CrossSectionBH(Q2Input,xBInput,tInput,phi,BeamHeli,TargetPolar,BeamCharge)--*
 | Computes the expansion of the Bethe Heitler cross section in terms of the "expansion |
 | cross sections" computed in the function MakeExactCrossSections                      | 
 | (or MakeLeadingCrossSections).                                                       |
 *--------------------------------------------------------------------------------------*/
 
 Double_t TGVKelly::CrossSectionBH( Double_t Q2Input, Double_t xBInput, Double_t tInput, Double_t phi, Int_t BeamHeli, Int_t TargetPolar, Bool_t Exact ) 
 {
 	Double_t SA, CS;
 	
 	SA = TGVKelly::SqrAmplBH( Q2Input, xBInput, tInput, phi, BeamHeli, TargetPolar, Exact );
 	
 	CS = SA * PhaseSpace;
 	
 	return CS; 	
 } // end def CrossSectionBH
 
 
  	
/*-Function CrossSectionVCS(Q2Input,xBInput,tInput,phi,BeamHeli,TargetPolar,BeamCharge)-*
 | Computes the expansion of the Virtual Compton Scattering cross section in terms of   |
 | the "expansion cross sections" computed in the function MakeExactCrossSections       | 
 | (or MakeLeadingCrossSections).                                                       |
 *--------------------------------------------------------------------------------------*/
 
 Double_t TGVKelly::CrossSectionVCS( Double_t Q2Input, Double_t xBInput, Double_t tInput, Double_t phi, Int_t BeamHeli, Int_t TargetPolar, Double_t ReH, Double_t ImH, Double_t ReE, Double_t ImE, Double_t ReHT, Double_t ImHT, Double_t ReET, Double_t ImET, Bool_t Exact ) 
 {
 	Double_t SA, CS;
 	
 	SA = TGVKelly::SqrAmplVCS( Q2Input, xBInput, tInput, phi, BeamHeli, TargetPolar, ReH, ImH, ReE, ImE, ReHT, ImHT, ReET, ImET, Exact );
 	
 	CS = SA * PhaseSpace;
 	
 	return CS; 	 	
 } // end def CrossSectionVCS
 
 
 
/*Fonction CrossSectionInterf(Q2Input,xBInput,tInput,phi,BeamHeli,TargetPolar,BeamCharge)*
 | Computes the expansion of the Interference cross section in terms of the "expansion |
 | cross sections" computed in the function MakeExactCrossSections                     |
 | (or MakeLeadingCrossSections).                                                      |
 *-------------------------------------------------------------------------------------*/
  
 Double_t TGVKelly::CrossSectionInterf( Double_t Q2Input, Double_t xBInput, Double_t tInput, Double_t phi, Int_t BeamHeli, Int_t TargetPolar, Int_t BeamCharge, Double_t ReH, Double_t ImH, Double_t ReE, Double_t ImE, Double_t ReHT, Double_t ImHT, Double_t ReET, Double_t ImET, Bool_t Exact ) 
 {
 	Double_t SA, CS;
 	
 	SA = TGVKelly::SqrAmplInterf( Q2Input, xBInput, tInput, phi, BeamHeli, TargetPolar, BeamCharge, ReH, ImH, ReE, ImE, ReHT, ImHT, ReET, ImET, Exact );
 	
 	CS = SA * PhaseSpace;
 	
 	return CS; 	 	
 } // end def CrossSectionInterf
