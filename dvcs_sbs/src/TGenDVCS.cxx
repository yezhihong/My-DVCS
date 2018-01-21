//
// TGenDVCS.cxx, v1.0, Tue Aug 4 11:13:57
// Author: C. Munoz Camacho
//

#include <fstream>
#include <iostream>
#include <stdlib.h>
//#include "mydvcs.hh"

#ifndef __TGenDVCS__
#include "TGenDVCS.h"
#endif

using namespace std;

ClassImp(TGenDVCS)

////////////////////////////////////////////////////////////////////////////////
//
// Event generator base class
// 
////////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________
  TGenDVCS::TGenDVCS(Double_t Ebeam, Int_t TargType, UInt_t seed1, UInt_t seed2):TGenBase(Ebeam,TargType, seed1,seed2)
{
  // Default constructor
  // Initial 4-vectors are initialized

   // if(!fgIsInit) Init();

  cout<<"TGenDVCS constructor"<<endl;

  fdmom=0;
  frho=0;
  fTargType=TargType;
  if(fTargType!=1 && fTargType!=0 && fTargType!=2 && fTargType!=3){
    cout<<"Unknown target for g Event generator"<<endl;
    exit(1);
  }
  
  if(fTargType==0) fm=PMass();
  if(fTargType==1) fm=NMass();
  if(fTargType==2) fm=DMass();
  if(fTargType==3) fm=He3Mass();

  fpini=0;

  ftmax=0.;
  ftmin=0.;

  fxb=0;
  fQ2=0;
  fe=0;
  fq=0;
  fg=0;
  fp=0;
  fSigmaP=0;
  fSigmaM=0;

  //Default parameters for the cross-section calculation
  fphasespace=1;
  fprop=1;
  fb=1.;
  ftdep=1;
  ftcoef=0.8;
  fDD=1;
  fJu=0.3;
  fJd=-0.1;
  fdterm=1;
  fpipole=1;

  // Initialization of TGV + read datafile
  V=new Double_t***[8];
  for(Int_t i=0;i<8;i++){
    V[i]=new Double_t**[13];
    for(Int_t j=0;j<13;j++){
      V[i][j]=new Double_t*[80];
      for(Int_t k=0;k<80;k++) V[i][j][k]=new Double_t[90];
    }
  }

  CFF=new Double_t[8];
  ifstream f("./src/CFFoutput_LO.dat");
  Double_t dum1,dum2,dum3;
  for ( register unsigned int iQ2 = 0; iQ2 < 13; iQ2++ ) {
    for ( register unsigned int iXb = 0; iXb < 80; iXb++ ) {
      for ( register unsigned int it = 0; it < 90; it++ ) {
	f>>dum1>>dum2>>dum3;
	for(Int_t j=0;j<8;j++) f>>V[j][iQ2][iXb][it];
      }
    }
  }
  tgv=new TGVKelly(Ebeam,kFALSE,kTRUE);
}

//_____________________________________________________________________________
Double_t TGenDVCS::KellyE(Double_t q2)
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

//_____________________________________________________________________________
Double_t TGenDVCS::KellyM(Double_t q2)
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

//_____________________________________________________________________________
Double_t* TGenDVCS::Interpol_CFF(Double_t Q2, Double_t xb, Double_t t)
{

  //  if(Q2<1||Q2>10||xb<0.2||xb>0.7||t>tmin||t<tmin-1)
  // cout<<Q2<<" "<<xb<<" "<<t<<" ==========="<<endl;
  Double_t eps2=4.*xb*xb*0.938272*0.938272/Q2;
  Double_t tmin = -Q2*(2.*(1.-xb)*(1-TMath::Sqrt(1+eps2))+eps2)/(4.*xb*(1.-xb)+eps2);

  if(Q2<1||Q2>13||xb<0.1||xb>0.9||t<-3.){
    cout<<"Kinematics (Q2,xb,t,tmin) out of range for cross-section evaluation: "<<Q2<<" "<<xb<<" "<<t<<" "<<tmin<<endl;
    return 0;
  }

  Int_t Q2_0=int(Q2-1.);
  Int_t Q2_1=int(Q2-1.)+1;
  Int_t xb_0=int( (xb-0.1)*79/0.8 );
  Int_t xb_1=int( (xb-0.1)*79/0.8 )+1;
  Int_t t_0=int( -(t-0.)*89/3. );
  Int_t t_1=int( -(t-0.)*89/3. )+1;

  Double_t Q2d=(Q2-(1+Q2_0))/((1+Q2_1)-(1+Q2_0));
  Double_t xbd=(xb-(0.1+0.8*xb_0/79.))/((0.1+0.8*xb_1/79.)-(0.1+0.8*xb_0/79.));
  Double_t td=(t-(0.-3.*t_0/89.))/((0.-3.*t_1/89.)-(0.-3.*t_0/89.));

  for(Int_t GPD=0;GPD<8;GPD++){
    
    Double_t c00=V[GPD][Q2_0][xb_0][t_0]*(1-Q2d)+V[GPD][Q2_1][xb_0][t_0]*Q2d;
    Double_t c10=V[GPD][Q2_0][xb_1][t_0]*(1-Q2d)+V[GPD][Q2_1][xb_1][t_0]*Q2d;
    Double_t c01=V[GPD][Q2_0][xb_0][t_1]*(1-Q2d)+V[GPD][Q2_1][xb_0][t_1]*Q2d;
    Double_t c11=V[GPD][Q2_0][xb_1][t_1]*(1-Q2d)+V[GPD][Q2_1][xb_1][t_1]*Q2d;
    
    Double_t c0=c00*(1-xbd)+c10*xbd;
    Double_t c1=c01*(1-xbd)+c11*xbd;
    
    CFF[GPD]=c0*(1-td)+c1*td;
  }
  return CFF;

}

//_____________________________________________________________________________
Double_t TGenDVCS::XSecSum(int opt)
{
  // returns d4sigma/dQ2dxdtdphi in nb/GeV4
  
  Double_t* check=Interpol_CFF(fQ2,fxb,ft);
  if(!check) return 0;
  
  TGVKelly *tgv2=new TGVKelly(fEbeam,kFALSE,kTRUE);

  Double_t ConvGeV2nbarn = 0.389379304e+6; // Changement d'unités !
  Double_t BHp, BHm, VCSp, VCSm, Ip, Im;
  Double_t SigmaTotPlus, SigmaTotMoins;
  BHp = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, 1, 0, kTRUE );
  VCSp = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, 1, 0, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Ip = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, 1, 0, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  BHm = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, -1, 0, kTRUE );
  VCSm = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, -1, 0, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Im = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, -1, 0, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  SigmaTotPlus = BHp + VCSp + Ip;
  SigmaTotMoins = BHm + VCSm + Im;

  // cout<<CFF[0]<<" "<<CFF[1]<<" "<<CFF[2]<<" "<<CFF[3]<<" "<<CFF[4]<<" "<<CFF[5]<<" "<<CFF[6]<<" "<<CFF[7]<<endl;
  // cout<<" here "<<fQ2<<" "<<fxb<<" "<<ft<<" "<<fphi<<" "<<TMath::Pi() * ( SigmaTotPlus + SigmaTotMoins ) * ConvGeV2nbarn<<endl;

  delete tgv2;

  if(opt==1) {
	  fSigmaP = TMath::Pi() * BHp * ConvGeV2nbarn;
	  fSigmaM = TMath::Pi() * BHm * ConvGeV2nbarn;

	  return TMath::Pi()*(BHp+BHm)* ConvGeV2nbarn;
  }

  fSigmaP = TMath::Pi() * SigmaTotPlus * ConvGeV2nbarn;
  fSigmaM = TMath::Pi() * SigmaTotMoins * ConvGeV2nbarn;
  //Added by Z. Ye to avoid nan
  if(isnan(fSigmaP))  fSigmaP = 1e-33;
  if(isnan(fSigmaM))  fSigmaM = 1e-33;


  fSigmaBHp = TMath::Pi() * BHp * ConvGeV2nbarn;
  fSigmaBHm = TMath::Pi() * BHm * ConvGeV2nbarn;
  if(isnan(fSigmaBHp))  fSigmaBHp = 1e-33;
  if(isnan(fSigmaBHm))  fSigmaBHm = 1e-33;

  return (fSigmaP +fSigmaM );
  
}

//_____________________________________________________________________________
Double_t TGenDVCS::XSecDif(void)
{
  // returns d4Sigma/dQ2dxdtdphi in nb/GeV4
  
  Double_t* check=Interpol_CFF(fQ2,fxb,ft);
  if(!check) return 0;

  TGVKelly *tgv2=new TGVKelly(fEbeam,kFALSE,kTRUE);

  Double_t ConvGeV2nbarn = 0.389379304e+6; // Changement d'unités !
  Double_t SigmaTotPlus, SigmaTotMoins;

  SigmaTotPlus = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, 1, 0, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  SigmaTotMoins = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, -1, 0, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );

  delete tgv2;

  //Added by Z. Ye to avoid nan
  if(isnan(SigmaTotPlus))  SigmaTotPlus = 1e-33;
  if(isnan(SigmaTotMoins))  SigmaTotMoins = 1e-33;

  return  TMath::Pi() * ( SigmaTotPlus - SigmaTotMoins ) * ConvGeV2nbarn;
  
}
//_____________________________________________________________________________
Double_t TGenDVCS::XSecDouble(int opt)
{
   // New subroutine added by to calculate XS with beam polarized and/or target polarized.
   // -- Z. Ye 03/25/2015
  // returns d4sigma/dQ2dxdtdphi in nb/GeV4
  Double_t* check=Interpol_CFF(fQ2,fxb,ft);
  if(!check) return 0;
  
  Double_t ConvGeV2nbarn = 0.389379304e+6; // Changement d'unités !
  Int_t BeamPor = 0; //0->unpolarized, 1->longitudinally polarized, -1-> longitudinally anti-polarized
  Int_t TargetPor = 0; //0->unpolarized, 1-> transversely polarized on x, 2->transversely polarized on y, 3->longitudinally polarized, 

  TGVKelly *tgv2=new TGVKelly(fEbeam,kFALSE,kTRUE);
  
  /*SSA{{{*/
  Double_t BH_LpU, BH_LmU, VCS_LpU, VCS_LmU, I_LpU, I_LmU;
  Double_t Sigma_LpU, Sigma_LmU;
  
  Double_t BH_ULp, BH_ULm, VCS_ULp, VCS_ULm, I_ULp, I_ULm;
  Double_t Sigma_ULp, Sigma_ULm;
 
  Double_t BH_UTXp, BH_UTXm, VCS_UTXp, VCS_UTXm, I_UTXp, I_UTXm;
  Double_t Sigma_UTXp, Sigma_UTXm;

  Double_t BH_UTYp, BH_UTYm, VCS_UTYp, VCS_UTYm, I_UTYp, I_UTYm;
  Double_t Sigma_UTYp, Sigma_UTYm;

  BeamPor = 1; TargetPor = 0;
  BH_LpU = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, kTRUE );
  VCS_LpU = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  I_LpU = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Sigma_LpU = BH_LpU + VCS_LpU + I_LpU;
  BeamPor = -1; TargetPor = 0;
  BH_LmU = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, kTRUE );
  VCS_LmU = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  I_LmU = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Sigma_LmU = BH_LmU + VCS_LmU + I_LmU;
 
  BeamPor = 0; TargetPor = 1;
  BH_UTXp = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, kTRUE );
  VCS_UTXp = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  I_UTXp = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Sigma_UTXp = BH_UTXp + VCS_UTXp + I_UTXp;
  BeamPor = 0; TargetPor = -1;
  BH_UTXm = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, kTRUE );
  VCS_UTXm = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  I_UTXm = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Sigma_UTXm = BH_UTXm + VCS_UTXm + I_UTXm;

  BeamPor = 0; TargetPor = 2;
  BH_UTYp = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, kTRUE );
  VCS_UTYp = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  I_UTYp = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Sigma_UTYp = BH_UTYp + VCS_UTYp + I_UTYp;
  BeamPor = 0; TargetPor = -2;
  BH_UTYm = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, kTRUE );
  VCS_UTYm = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  I_UTYm = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Sigma_UTYm = BH_UTYm + VCS_UTYm + I_UTYm;
  
  BeamPor = 0; TargetPor = 3;
  BH_ULp = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, kTRUE );
  VCS_ULp = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  I_ULp = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Sigma_ULp = BH_ULp + VCS_ULp + I_ULp;
  BeamPor = 0; TargetPor = -3;
  BH_ULm = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, kTRUE );
  VCS_ULm = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  I_ULm = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Sigma_ULm = BH_ULm + VCS_ULm + I_ULm;

  fSigma_LpU = TMath::Pi() * Sigma_LpU * ConvGeV2nbarn;
  fSigma_LmU = TMath::Pi() * Sigma_LmU * ConvGeV2nbarn;
  fSigma_ULp = TMath::Pi() * Sigma_ULp * ConvGeV2nbarn;
  fSigma_ULm = TMath::Pi() * Sigma_ULm * ConvGeV2nbarn;
  fSigma_UTXp = TMath::Pi() * Sigma_UTXp * ConvGeV2nbarn;
  fSigma_UTXm = TMath::Pi() * Sigma_UTXm * ConvGeV2nbarn;
  fSigma_UTYp = TMath::Pi() * Sigma_UTYp * ConvGeV2nbarn;
  fSigma_UTYm = TMath::Pi() * Sigma_UTYm * ConvGeV2nbarn;
  /*}}}*/

  /*DSA{{{*/
  Double_t BH_LpLp, BH_LmLp, VCS_LpLp, VCS_LmLp, I_LpLp, I_LmLp;
  Double_t Sigma_LpLp, Sigma_LmLp;
  Double_t BH_LpLm, BH_LmLm, VCS_LpLm, VCS_LmLm, I_LpLm, I_LmLm;
  Double_t Sigma_LpLm, Sigma_LmLm;

  Double_t BH_LpTXp, BH_LmTXp, VCS_LpTXp, VCS_LmTXp, I_LpTXp, I_LmTXp;
  Double_t Sigma_LpTXp, Sigma_LmTXp;
  Double_t BH_LpTXm, BH_LmTXm, VCS_LpTXm, VCS_LmTXm, I_LpTXm, I_LmTXm;
  Double_t Sigma_LpTXm, Sigma_LmTXm;
   
  Double_t BH_LpTYp, BH_LmTYp, VCS_LpTYp, VCS_LmTYp, I_LpTYp, I_LmTYp;
  Double_t Sigma_LpTYp, Sigma_LmTYp;
  Double_t BH_LpTYm, BH_LmTYm, VCS_LpTYm, VCS_LmTYm, I_LpTYm, I_LmTYm;
  Double_t Sigma_LpTYm, Sigma_LmTYm;

  BeamPor = 1; TargetPor = 1;
  BH_LpTXp = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, kTRUE );
  VCS_LpTXp = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  I_LpTXp = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Sigma_LpTXp = BH_LpTXp + VCS_LpTXp + I_LpTXp;
  BeamPor = 1; TargetPor =-1;
  BH_LpTXm = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, kTRUE );
  VCS_LpTXm = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  I_LpTXm = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Sigma_LpTXm = BH_LpTXm + VCS_LpTXm + I_LpTXm;

  BeamPor = -1; TargetPor = 1;
  BH_LmTXp = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, kTRUE );
  VCS_LmTXp = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  I_LmTXp = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Sigma_LmTXp = BH_LmTXp + VCS_LmTXp + I_LmTXp;
  BeamPor = -1; TargetPor = -1;
  BH_LmTXm = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, kTRUE );
  VCS_LmTXm = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  I_LmTXm = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Sigma_LmTXm = BH_LmTXm + VCS_LmTXm + I_LmTXm;


  BeamPor = 1; TargetPor = 2;
  BH_LpTYp = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, kTRUE );
  VCS_LpTYp = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  I_LpTYp = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Sigma_LpTYp = BH_LpTYp + VCS_LpTYp + I_LpTYp;
  BeamPor = 1; TargetPor =-2;
  BH_LpTYm = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, kTRUE );
  VCS_LpTYm = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  I_LpTYm = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Sigma_LpTYm = BH_LpTYm + VCS_LpTYm + I_LpTYm;

  BeamPor = -1; TargetPor = 2;
  BH_LmTYp = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, kTRUE );
  VCS_LmTYp = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  I_LmTYp = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Sigma_LmTYp = BH_LmTYp + VCS_LmTYp + I_LmTYp;
  BeamPor = -1; TargetPor = -2;
  BH_LmTYm = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, kTRUE );
  VCS_LmTYm = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  I_LmTYm = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Sigma_LmTYm = BH_LmTYm + VCS_LmTYm + I_LmTYm;


  BeamPor = 1; TargetPor = 3;
  BH_LpLp = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, kTRUE );
  VCS_LpLp = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  I_LpLp = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Sigma_LpLp = BH_LpLp + VCS_LpLp + I_LpLp;
  BeamPor = 1; TargetPor =-3;
  BH_LpLm = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, kTRUE );
  VCS_LpLm = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  I_LpLm = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Sigma_LpLm = BH_LpLm + VCS_LpLm + I_LpLm;

  BeamPor = -1; TargetPor = 3;
  BH_LmLp = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, kTRUE );
  VCS_LmLp = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  I_LmLp = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Sigma_LmLp = BH_LmLp + VCS_LmLp + I_LmLp;
  BeamPor = -1; TargetPor = -3;
  BH_LmLm = tgv2->CrossSectionBH( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, kTRUE );
  VCS_LmLm = tgv2->CrossSectionVCS( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  I_LmLm = tgv2->CrossSectionInterf( fQ2, fxb, ft, -fphi, BeamPor, TargetPor, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
  Sigma_LmLm = BH_LmLm + VCS_LmLm + I_LmLm;


  fSigma_LpTXp = TMath::Pi() * Sigma_LpTXp * ConvGeV2nbarn;
  fSigma_LmTXp = TMath::Pi() * Sigma_LmTXp * ConvGeV2nbarn;
  fSigma_LpTXm = TMath::Pi() * Sigma_LpTXm * ConvGeV2nbarn;
  fSigma_LmTXm = TMath::Pi() * Sigma_LmTXm * ConvGeV2nbarn;

  fSigma_LpTYp = TMath::Pi() * Sigma_LpTYp * ConvGeV2nbarn;
  fSigma_LmTYp = TMath::Pi() * Sigma_LmTYp * ConvGeV2nbarn;
  fSigma_LpTYm = TMath::Pi() * Sigma_LpTYm * ConvGeV2nbarn;
  fSigma_LmTYm = TMath::Pi() * Sigma_LmTYm * ConvGeV2nbarn;

  fSigma_LpLp = TMath::Pi() * Sigma_LpLp * ConvGeV2nbarn;
  fSigma_LmLp = TMath::Pi() * Sigma_LmLp * ConvGeV2nbarn;
  fSigma_LpLm = TMath::Pi() * Sigma_LpLm * ConvGeV2nbarn;
  fSigma_LmLm = TMath::Pi() * Sigma_LmLm * ConvGeV2nbarn;
  /*}}}*/


  // cout<<CFF[0]<<" "<<CFF[1]<<" "<<CFF[2]<<" "<<CFF[3]<<" "<<CFF[4]<<" "<<CFF[5]<<" "<<CFF[6]<<" "<<CFF[7]<<endl;
  // cout<<" here "<<fQ2<<" "<<fxb<<" "<<ft<<" "<<fphi<<" "<<TMath::Pi() * ( SigmaTotPlus + SigmaTotMoins ) * ConvGeV2nbarn<<endl;

  delete tgv2;

  /*return{{{*/
  //Added by Z. Ye to avoid nan
  if(isnan(fSigma_LpU))  fSigma_LpU = 1e-33;
  if(isnan(fSigma_LmU))  fSigma_LmU = 1e-33;
  if(isnan(fSigma_ULp))  fSigma_ULp = 1e-33;
  if(isnan(fSigma_ULm))  fSigma_ULm = 1e-33;
  if(isnan(fSigma_LpLp))  fSigma_LpLp = 1e-33;
  if(isnan(fSigma_LpLm))  fSigma_LpLm = 1e-33;
  if(isnan(fSigma_LmLp))  fSigma_LmLp = 1e-33;
  if(isnan(fSigma_LmLm))  fSigma_LmLm = 1e-33;
  if(isnan(fSigma_LpTXp))  fSigma_LpTXp = 1e-33;
  if(isnan(fSigma_LpTXm))  fSigma_LpTXm = 1e-33;
  if(isnan(fSigma_LmTXp))  fSigma_LmTXp = 1e-33;
  if(isnan(fSigma_LmTXm))  fSigma_LmTXm = 1e-33;
  if(isnan(fSigma_LpTYp))  fSigma_LpTYp = 1e-33;
  if(isnan(fSigma_LpTYm))  fSigma_LpTYm = 1e-33;
  if(isnan(fSigma_LmTYp))  fSigma_LmTYp = 1e-33;
  if(isnan(fSigma_LmTYm))  fSigma_LmTYm = 1e-33;

  if(opt==1) 
	  return (fSigma_LpU +fSigma_LpU );
  if(opt==2) 
	  return (fSigma_ULp +fSigma_ULm );

  if(opt==11) 
	  return (fSigma_LpTXp +fSigma_LpTXm );
  if(opt==12) 
	  return (fSigma_LpTYp +fSigma_LpTYm );
  if(opt==13) 
	  return (fSigma_LpLp +fSigma_LpLm );

  if(opt==22) 
	  return (fSigma_LmTXp +fSigma_LmTXm );
  if(opt==23) 
	  return (fSigma_LmTYp +fSigma_LmTYm );
  if(opt==23) 
	  return (fSigma_LmLp +fSigma_LmLm );
  /*}}}*/
}


//_____________________________________________________________________________
TGenDVCS::TGenDVCS(const TGenDVCS& gen):TGenBase(gen)
{
  // Copy constructor
  //  ((TGenDVCS&)TCalobase).Copy(*this);
}

//_____________________________________________________________________________
 TGenDVCS::~TGenDVCS()
{
  // Default destructor

  if(fpini) delete fpini;
  if(fe)    delete fe;
  if(fp)    delete fp;
  if(fg)    delete fg;
}

//_____________________________________________________________________________
void TGenDVCS::IntRCBef(void)
{
  // Makes internal radiative corrections _before_ the vertex (to the initial)
  // electron) using the equivalent radiator technique.
  // Uses \Delta E=E_0 * R^(2./nu) with R randomly between 0 and 1.
  // Factor 1/2 is because internal corrections must be applied twice (before
  // and after the vertex) with equivalent radiator thickness half each time.

  fRadCor=kTRUE;
  Double_t nu=2.*0.007297352533*(TMath::Log(fQ2/TMath::Power(0.000510998902,2))-1.)/TMath::Pi(); 
  
  Double_t eel=(feini->E())*(1-TMath::Power(fRan->Rndm(),1./(nu/2.)));
  
  feini->SetPxPyPzE(0.,0.,eel,eel);
}

//_____________________________________________________________________________
void TGenDVCS::IntRCAft(void)
{
  // Makes internal radiative corrections _after_ the vertex (to the scattered
  // electron) using the equivalent radiator technique.
  // Uses \Delta E=E_0 * R^(2./nu) with R randomly between 0 and 1.
  // Factor 1/2 is because internal corrections must be applied twice (before
  // and after the vertex) with equivalent radiator thickness half each time.

  fRadCor=kTRUE;
  Double_t nu=2.*0.007297352533*(TMath::Log(fQ2/TMath::Power(0.000510998902,2))-1.)/TMath::Pi(); 
  
  Double_t deel=1-TMath::Power(fRan->Rndm(),1./(nu/2.));

  *feprerad=(*fe); // keep vertex scattered electron

  *fe=deel*(*fe);
}
//_____________________________________________________________________________
  void TGenDVCS::Settmin(Double_t tmin) { 
    //BEWARE : When generating events in a loop, the TGenDVCS::Clear() method,
    //which needs to be used at the end of the loop, sets tmin=-2. and tmax=0.
    //If you don't want these default values, you need to use TGenDVCS::Settmin()
    //and TGenDVCS::Settmax() at the beginning of the loop to set your values
ftmin=tmin ; 
}
//_____________________________________________________________________________
  void TGenDVCS::Settmax(Double_t tmax) {
    //BEWARE : When generating events in a loop, the TGenDVCS::Clear() method,
    //which needs to be used at the end of the loop, sets tmin=-2. and tmax=0.
    //If you don't want these default values, you need to use TGenDVCS::Settmin()
    //and TGenDVCS::Settmax() at the beginning of the loop to set your values
 ftmax=tmax ; 
}
//_____________________________________________________________________________
  void TGenDVCS::Clear(char* opt) 
    {
      // Sets initial electron and proton and sets tmin=-2 GeV and tmax=0 GeV
      // If you don't want these default values of tmin and tmax use 
      // TGenDVCS::Settmin() and TGenDVCS::Settmax() at the beginning of your
      // loop

      TGenBase::Clear(opt);
      ftmin=-2.;
      ftmax=0.;

    }

//_____________________________________________________________________________
 void TGenDVCS::ComputeDVCS(void) 
{ 
  //Computes the gamma* p -> gamma p' reaction in the center of mass.
  //The angle phi is generated and all vectors are boost back to the 
  //laboratory. Data members fp and fg are set.

  Compute2Body(0.);
}

//_____________________________________________________________________________
 void TGenDVCS::Compute2Body(Double_t m) 
{ 
  //Computes the gamma* p -> X p' reaction in the center of mass, as a
  //function of the mass m of the particle X.
  //The angle phi is generated and all vectors are boost back to the 
  //laboratory. Data members fp and fg are set.
  
  //////////////////////
  //Computation of t and the phase space factor
  //////////////////////

  Double_t nu=fQ2/(2.*fm*fxb);
  // BEWARE:Calculated for non-fermi smeared kinematics

  Double_t q3=TMath::Sqrt(fQ2+TMath::Power(nu,2.));

  Double_t q0primemax=0.5*fQ2*(1.-fxb)/(fxb*(fm+nu-q3));
  Double_t q0primemin=0.5*fQ2*(1.-fxb)/(fxb*(fm+nu+q3));
  Double_t tmax=-fQ2-2.*q0primemax*(nu-q3);
  Double_t tmin=-fQ2-2.*q0primemin*(nu+q3);
  ftmax=TMath::Min(ftmax,tmax);
  ftmin=TMath::Max(ftmin,tmin);
  
  if(fFermi) ftmax=0;

  do {
    ft=ftmax+(ftmin-ftmax)*fRan->Rndm();
  } while (ft/fQ2>1.);

  //  cout<<ft<<" "<<ftmax<<" "<<ftmin<<endl;


  fPSF=(ftmax-ftmin)*(fxbmax-fxbmin)*(fQ2max-fQ2min);
  ///////////////////

  //  cout << "ftmax: " << ftmax << endl;
  //  cout << "ftmin: " << ftmin << endl;

  // fq=new TLorentzVector(feini->Px()-fe->Px(),feini->Py()-fe->Py(),feini->Pz()-fe->Pz(),feini->E()-fe->E());

  if(fpini==0 && fFermi==kFALSE) fpini=new TLorentzVector(0.,0.,0.,fm);
  if(fpini==0 && fFermi==kTRUE){
    cout<<"You must generate fermi recoil particle first"<<endl;
    exit(1);
  }
  
  TLorentzVector p=*fpini;
  TLorentzVector q=*fq;

  // Boost to center of mass p-q

  TLorentzVector cms = q + p;

  fs=cms.M2();

  q.Boost(-cms.BoostVector() );
  p.Boost(-cms.BoostVector() );

  // Rotation around Y
  
  TVector3 oz(0.,0.,1.);
  TVector3 perpvec=(q.Vect()).Cross(oz);
  Double_t angle=(q.Vect()).Angle(oz);

  q.Rotate(angle,perpvec.Unit());
  p.Rotate(angle,perpvec.Unit());

  //Double_t egammacm=(fs-fm*fm)/(2.*TMath::Sqrt(fs));
  //Double_t thetacm = (ft+fQ2)/(2.*egammacm*q.P())+q(3)/q.P();

  //////////////////////////////////////////

  Double_t ecm=(fs-fm*fm+m*m)/(2.*TMath::Sqrt(fs));
  Double_t pcm=TMath::Sqrt(TMath::Power(ecm,2.)-TMath::Power(m,2.));
  Double_t thetacm=(ft+fQ2-m*m+2.*q(3)*ecm)/(2.*pcm*q.P());

  ///////////////////////////////////////// 

  if (thetacm>1.) thetacm=1.;
  if (thetacm<-1.) thetacm=-1.;
  thetacm=TMath::ACos(thetacm);

  if(!fg) {
    fg=new TLorentzVector(pcm*TMath::Sin(thetacm),0.,TMath::Cos(thetacm)*pcm,ecm);
  }else{
    fg->SetPxPyPzE(pcm*TMath::Sin(thetacm),0.,TMath::Cos(thetacm)*pcm,ecm);
  }
  
  if(!fp){
    fp=new TLorentzVector(q.Px()+p.Px()-fg->Px(),q.Py()+p.Py()-fg->Py(),q.Pz()+p.Pz()-fg->Pz(),q.E()+p.E()-fg->E());
  }else{
    fp->SetPxPyPzE(q.Px()+p.Px()-fg->Px(),q.Py()+p.Py()-fg->Py(),q.Pz()+p.Pz()-fg->Pz(),q.E()+p.E()-fg->E());
  }
   
  // Rotation back
  fp->Rotate(-angle,perpvec.Unit());
  fg->Rotate(-angle,perpvec.Unit());
  q.Rotate(-angle,perpvec.Unit());
  
  // Boost to lab
  fp->Boost(cms.BoostVector());
  fg->Boost(cms.BoostVector());
  q.Boost(cms.BoostVector());

  // Now we have k along Oz and k' on the xOz plane
  // q=k-k' is necessarily on plane as well
  // We just have to rotate along Oy to bring q along Oz
  // and then rotate everything by phi except
  // k and k'. Then we rotate back along Oy.

  angle=(q.Vect()).Angle(oz);
  //  cout<<angle<<endl;

  fp->RotateY(angle);
  fg->RotateY(angle);
  
  fphi=2.*TMath::Pi()*fRan->Rndm();   // phi between 0 and 2pi
  
  fp->RotateZ(fphi);
  fg->RotateZ(fphi);

  fp->RotateY(-angle);
  fg->RotateY(-angle);

}

//_____________________________________________________________________________
void TGenDVCS::ApplyAcc(Double_t avv, Double_t AzimAccMin, Double_t AzimAccMax)
{
  // Applies vertical spectrometer acceptance by rotating all 4-vectors
  // around the beam axis. An angle can be specified, otherwise it's generated
  // randomly between spectrometer acceptances

  const Double_t fSpecPhiAcc = (AzimAccMax-AzimAccMin)/2.;//rad
  Double_t av= avv;
  if(fabs(avv+1)<0.01)//aav=-1
	  av=-fSpecPhiAcc+2.*fSpecPhiAcc*fRan->Rndm();
  
  TGenBase::ApplySpecVerAcc(av);
  fg->RotateZ(av);
  fp->RotateZ(av);
  //fe->RotateZ(av);
  fq->RotateZ(av);
  
  fPSF*=2.*fSpecPhiAcc;

  // Remark:
  // The right phasespace factor is
  // dx.dq2.dt.dphi.dphie/2pi (because the x-section code gives 
  // phi_e-integrated x-sections) but dphi=2pi, so it simplifies to
  // dx.dq2.dt.dphie
}

//_____________________________________________________________________________
void TGenDVCS::ApplySpecVerAcc(Double_t aav)
{
  // Applies vertical spectrometer acceptance by rotating all 4-vectors
  // around the beam axis. An angle can be specified, otherwise it's generated
  // randomly between spectrometer acceptances

  const Double_t fSpecPhiAcc = fSpecVerAcc/TMath::Sin(fSpecAngle-fSpecHorAcc);
  Double_t av=aav;
  if(aav==-1.) {
    av=-fSpecPhiAcc+2.*fSpecPhiAcc*fRan->Rndm();
  }
  TGenBase::ApplySpecVerAcc(av);
  fg->RotateZ(av);
  fp->RotateZ(av);
  //fe->RotateZ(av);
  fq->RotateZ(av);
  

  fPSF*=2.*fSpecPhiAcc;

  // Remark:
  // The right phasespace factor is
  // dx.dq2.dt.dphi.dphie/2pi (because the x-section code gives 
  // phi_e-integrated x-sections) but dphi=2pi, so it simplifies to
  // dx.dq2.dt.dphie
}

//_____________________________________________________________________________
 TLorentzVector* TGenDVCS::GetFinalPhoton(void)
{
  // Returns the final photon 4-vector if it exists
  
  if(!fg) cout<<"Warning : Final photon doesn't exist"<<endl;
  return fg;

}

//_____________________________________________________________________________
 TLorentzVector* TGenDVCS::GetFinalProton(void)
{
  // Returns the final photon 4-vector if it exists
  
  if(!fp) cout<<"Warning : Final proton doesn't exist"<<endl;
  return fp;

}

//_____________________________________________________________________________
 Double_t TGenDVCS::GetFastWeight(void)
{
  // Returns a fast cross-section

  Double_t q2c=2.*fSpecMom*fEbeam*(1-TMath::Cos(fSpecAngle));
  Double_t weight=TMath::Power(q2c/fQ2,3.)*TMath::Power(1./(1.-ft/(fm*fm)),4.);

  return weight;
}

//_____________________________________________________________________________
 void TGenDVCS::XSec(void)
{
  // Computes the g cross-section for both electron helicities

  //  SigmaDVCS(&fSigmaP,&fSigmaM,fTargType,fFermi,feini->E(),fQ2,fxb,ft,fphi,fpini->P(),fpini->Theta(),fpini->Phi(),fphasespace,fprop,fb,ftdep,ftcoef,fDD,fJu,fJd,fdterm,fpipole);

  //  if (TMath::IsNaN(fSigmaP) || TMath::IsNaN(fSigmaM)) {fSigmaP=0.;fSigmaM=0.;}

}

//_____________________________________________________________________________
 void TGenDVCS::SetTheoryParam(Int_t phasespace, Int_t prop, Double_t b, Int_t tdep, Double_t tcoef, Int_t DD, Double_t Ju, Double_t Jd, Int_t dterm, Int_t pipole)
{
  //Sets parameters for the library used to compute the g cross-section

  fphasespace=phasespace;
  fprop=prop;
  fb=b;
  ftdep=tdep;
  ftcoef=tcoef;
  fDD=DD;
  fJu=Ju;
  fJd=Jd;
  fdterm=dterm;
  fpipole=pipole;

}

//_____________________________________________________________________________
void TGenDVCS::Write2File(void)
{

  if(!fNwrite) {
    if(!fFermi) {
      fNwrite=17+6;
    }else{
      fNwrite=17+6+3;
    }
  }
  
  *output<<fVertex->Pz()<<" ";
  *output<<feini->Pz()<<" ";
  *output<<feprerad->Px()<<" "<<feprerad->Py()<<" "<<feprerad->Pz()<<" ";
  if(fFermi)
    *output<<fpini->Px()<<" "<<fpini->Py()<<" "<<fpini->Pz()<<" ";
  *output<<fe->Px()<<" "<<fe->Py()<<" "<<fe->Pz()<<" ";
  *output<<fg->Px()<<" "<<fg->Py()<<" "<<fg->Pz()<<" ";
  *output<<fp->Px()<<" "<<fp->Py()<<" "<<fp->Pz()<<" ";
  *output<<fQ2<<" "<<fxb<<" "<<ft<<" "<<fphi<<" ";//add four kin. variables, Z. Ye 03/11/2015
  *output<<fSigmaP<<" "<<fSigmaM<<" "<<fSigmaBHp<<" "<<fSigmaBHm<<" "<<fPSF<<endl;
  //  *output<<fVertez->Pz()<<endl;
  // *output<<endl;
}

//_____________________________________________________________________________
 void TGenDVCS::Print(char* opt)
{
  // Output on screen. If option "all" is specified the complete setup of the
  // event is printed out. By default only the final state 4-vectors and the 
  // virtual photon are displayed.

  TString option=opt;
  if(option.Contains("all")) TGenBase::Print();
  cout<<"======================================="<<endl;
  cout<<"        4-vectors (Px,Py,Pz,E)         "<<endl;
  cout<<"======================================="<<endl;
  if(feini) {
    cout<<"e("<<feini->Px()<<","<<feini->Py()<<","<<feini->Pz()<<","<<feini->E()<<")"<<endl;
  }else{
    cout<<"NO INITIAL ELECTRON DEFINED"<<endl;
  }
  if(fpini){
    cout<<"p("<<fpini->Px()<<","<<fpini->Py()<<","<<fpini->Pz()<<","<<fpini->E()<<")"<<endl;
  }else{
    cout<<"NO INITIAL TARGET PARTICLE DEFINED"<<endl;
  }
  if(fq){
  cout<<"g*("<<fq->Px()<<","<<fq->Py()<<","<<fq->Pz()<<","<<fq->E()<<")"<<endl;
  }else{
    cout<<"NO VIRTUAL PHOTON DEFINED"<<endl;
  }
  if(fe){
    cout<<"e'("<<fe->Px()<<","<<fe->Py()<<","<<fe->Pz()<<","<<fe->E()<<")"<<endl;
  }else{
    cout<<"NO SCATTERED ELECTRON DEFINED"<<endl;
  }
  if(fg){
    cout<<"g("<<fg->Px()<<","<<fg->Py()<<","<<fg->Pz()<<","<<fg->E()<<")"<<endl;
  }else{
    cout<<"NO EMITTED PHOTON DEFINED"<<endl;
  }
  if(fp){
    cout<<"p'("<<fp->Px()<<","<<fp->Py()<<","<<fp->Pz()<<","<<fp->E()<<")"<<endl;
  }else{
    cout<<"NO RECOIL PARTICLE DEFINED"<<endl;
  }
  cout<<"======================================="<<endl;
}
