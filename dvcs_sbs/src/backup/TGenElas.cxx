//
// TGenElas.cxx, v1.0, Tue Aug 4 11:13:57
// Author: C. Munoz Camacho
//

#include <fstream>
#include <iostream>
#include <stdlib.h>

#ifndef __TGenElas__
#include "TGenElas.h"
#endif

using namespace std;

ClassImp(TGenElas)

////////////////////////////////////////////////////////////////////////////////
//
// Event generator base class
// 
////////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________
  TGenElas::TGenElas(Double_t Ebeam, Int_t TargType, UInt_t seed1, UInt_t seed2):TGenBase(Ebeam,seed1,seed2)
{
  // Default constructor
  // Initial 4-vectors are initialized

   // if(!fgIsInit) Init();

  cout<<"TGenElas constructor"<<endl;

  fdmom=0;
  frho=0;
  fTargType=TargType;

  if(fTargType!=1 && fTargType!=0){
    cout<<"Unknown target for DVCS Event generator"<<endl;
    exit(1);
  }

  if(fTargType==0) fm=PMass();
  if(fTargType==1) fm=NMass();

  fpini=0;

  fxb=0;
  fQ2=0;
  fe=0;
  fq=0;
  fp=0;
    
}

//_____________________________________________________________________________
TGenElas::TGenElas(const TGenElas& gen):TGenBase(gen)
{
  // Copy constructor
  //  ((TGenElas&)TCalobase).Copy(*this);
}

//_____________________________________________________________________________
 TGenElas::~TGenElas()
{
  // Default destructor

  if(fpini) delete fpini;
  if(fe)    delete fe;
  if(fp)    delete fp;

}

//_____________________________________________________________________________
void TGenElas::IntRCBef(void)
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
void TGenElas::IntRCAft(void)
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
 void TGenElas::ComputeElas(void) 
{ 
  //Computes the gamma* p -> p' reaction in the center of mass.
  //Initially generated electron, proton and virtual photon are all rotated
  //by 180 deg along the beam axis in order to generate an electron to be 
  //detected in the calorimeter and a proton in the spectrometer.
  //
  //All vectors are boost back to the laboratory.
  
  fPSF=(fQ2max-fQ2min);

  if(fpini==0 && fFermi==kFALSE) fpini=new TLorentzVector(0.,0.,0.,fm);
  if(fpini==0 && fFermi==kTRUE){
    cout<<"You must generate fermi recoil particle first"<<endl;
    exit(1);
  }

  //In the elastic case, we detect the electron on the calorimeter and the
  //proton in the spectrometer. TGenBase::ComputeElectron() generates the 
  //electron to the left. For the elastic case, we rotate everything 180 deg.
  //here. We could have also overriden TGenBase::ComputeElectron(), but also
  //TGenBase::GenFermiIni()...
  fe->RotateZ(TMath::Pi());
  feprerad->RotateZ(TMath::Pi());
  fq->RotateZ(TMath::Pi());
  fpini->RotateZ(TMath::Pi());
  //////////////////////////////

  if(!fp){
    fp=new TLorentzVector(fq->Px()+fpini->Px(),fq->Py()+fpini->Py(),fq->Pz()+fpini->Pz(),fq->E()+fpini->E());
  }else{
    fp->SetPxPyPzE(fq->Px()+fpini->Px(),fq->Py()+fpini->Py(),fq->Pz()+fpini->Pz(),fq->E()+fpini->E());
  }
}

//_____________________________________________________________________________
void TGenElas::GenKin(void)
{
  // We overload the TGenBase::GenKin() to constrain xb=1 (elastic event)
  // instead of generating it randomly and generate Q2 between the right elastic
  // limits

  Double_t thetapmax=fSpecAngle+fSpecHorAcc;
  Double_t thetapmin=fSpecAngle-fSpecHorAcc;
  Double_t ppmax=fSpecMom*(1+fSpecMomAcc);
  Double_t ppmin=fSpecMom*(1-fSpecMomAcc);

  Double_t eel=feini->E();

  Double_t Q2minP=2.*fm*fm*(TMath::Sqrt(1+TMath::Power(ppmin/fm,2))-1.);
  Double_t Q2maxP=2.*fm*fm*(TMath::Sqrt(1+TMath::Power(ppmax/fm,2))-1.);
  
  Double_t Q2maxA=-(4.*fm*fm*TMath::Power(TMath::Cos(thetapmin),2.))/(TMath::Power(TMath::Cos(thetapmin),2.)-TMath::Power(1+fm/eel,2.));
  Double_t Q2minA=-(4.*fm*fm*TMath::Power(TMath::Cos(thetapmax),2.))/(TMath::Power(TMath::Cos(thetapmax),2.)-TMath::Power(1+fm/eel,2.));
  
  fQ2min=TMath::Max(Q2minP,Q2minA);
  fQ2max=TMath::Min(Q2maxP,Q2maxA);
  
  fQ2=fQ2min+(fQ2max-fQ2min)*fRan->Rndm();
  
  fxb=1.;

}

//_____________________________________________________________________________
Double_t TGenElas::XSec(void)
{
  //Computes the elastic diferential cross-section for the event

  TVector3 Oz(0.,0.,1.);
  Double_t thetae=Oz.Angle(fe->Vect());

  Double_t nscs=TMath::Power(Alpha()*TMath::Cos(thetae/2),2)/(4.*TMath::Power(feini->E(),2)*TMath::Power(TMath::Sin(thetae/2.),4)*(1+2.*feini->E()*TMath::Power(TMath::Sin(thetae/2.),2)/fm));//non-structured cross-section
  
  Double_t ge=TMath::Power(1+fQ2/0.71,-2);
  Double_t gm=2.792847386*ge;
  
  Double_t tau=fQ2/(4.*TMath::Power(fm,2));
  
  Double_t cs=nscs*((TMath::Power(ge,2)+tau*TMath::Power(gm,2))/(1+tau)+2.*tau*
  TMath::Power(gm,2)*TMath::Power(TMath::Tan(thetae/2.),2));//cross-section
  cs*=0.389379e-3; //cross-section in b
  cs*=1e9; // in nb
  return cs;
  
}

Double_t TGenElas::ArringE(Double_t q2)
{
  Double_t a[7]={1.,3.226,1.508,-0.3773,0.611,-0.1853,0.01596};
  Double_t GEa=0;

   for (int ja=0; ja<7; ja++){
     GEa+=a[ja]*pow(q2,ja);
   }
   GEa=1./GEa;

   return GEa;
}

Double_t TGenElas::ArringM(Double_t q2)
{
  Double_t a[7]={1.,3.19,1.355,0.151,-0.0114,0.000533,-0.000009};
  Double_t GMa=0;

   for (int ja=0; ja<7; ja++){
     GMa+=a[ja]*pow(q2,ja);
   }
   GMa=1./GMa;

   return 2.792847386*GMa;
}

Double_t TGenElas::KellyE(Double_t q2)
   /* JJ Kelly PRC 70, 068202 (2004)
    */
{
   Int_t ia=1, ib=3;
   Double_t a[1]={-0.24};
   Double_t b[3]={10.98, 12.82, 21.97};
   Double_t Mp = 0.938272;
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

Double_t TGenElas::KellyM(Double_t q2)
   /* JJ Kelly PRC 70, 068202 (2004)
      Magnetic Form Factor fit
      Returned value is ratio to dipole*mu_p
   */
{
   Int_t ia=1, ib=3;
   Double_t a[1]={0.12};
   Double_t b[3]={10.97, 18.86, 6.55};
   Double_t Mp = 0.938272;
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
Double_t TGenElas::XSecKelly(void)
{
  //Computes the elastic diferential cross-section for the event

  TVector3 Oz(0.,0.,1.);
  Double_t thetae=Oz.Angle(fe->Vect());

  Double_t nscs=TMath::Power(Alpha()*TMath::Cos(thetae/2),2)/(4.*TMath::Power(feini->E(),2)*TMath::Power(TMath::Sin(thetae/2.),4)*(1+2.*feini->E()*TMath::Power(TMath::Sin(thetae/2.),2)/fm));//non-structured cross-section
  
  Double_t ge=KellyE(-fQ2);
  Double_t gm=KellyM(-fQ2);
  
  Double_t tau=fQ2/(4.*TMath::Power(fm,2));
  
  Double_t cs=nscs*((TMath::Power(ge,2)+tau*TMath::Power(gm,2))/(1+tau)+2.*tau*TMath::Power(gm,2)*TMath::Power(TMath::Tan(thetae/2.),2));//cross-section
  cs*=0.389379e-3; //cross-section in bn
  cs*=1e9;
//   cout << "Kelly:  "<<cs <<endl;
//   cout << "Dipole: "<<XSec()<<endl;
  return cs;
  
}

//_____________________________________________________________________________
Double_t TGenElas::XSecArrington(void)
{
  //Computes the elastic diferential cross-section for the event

  TVector3 Oz(0.,0.,1.);
  Double_t thetae=Oz.Angle(fe->Vect());

  Double_t nscs=TMath::Power(Alpha()*TMath::Cos(thetae/2),2)/(4.*TMath::Power(feini->E(),2)*TMath::Power(TMath::Sin(thetae/2.),4)*(1+2.*feini->E()*TMath::Power(TMath::Sin(thetae/2.),2)/fm));//non-structured cross-section
  
  Double_t ge=ArringE(fQ2);
  Double_t gm=ArringM(fQ2);
  
  Double_t tau=fQ2/(4.*TMath::Power(fm,2));
  
  Double_t cs=nscs*((TMath::Power(ge,2)+tau*TMath::Power(gm,2))/(1+tau)+2.*tau*TMath::Power(gm,2)*TMath::Power(TMath::Tan(thetae/2.),2));//cross-section
  cs*=0.389379e-3; //cross-section in bn
  cs*=1e9;
//    cout << "Kelly:  "<<XSecKelly() <<endl;
//    cout << "Dipole: "<<XSec()<<endl;
//    cout << "Arring: "<<cs<<endl<<endl;

  return cs;
  
}

//_____________________________________________________________________________
Double_t TGenElas::XSecTotal(Int_t steps)
{
  //Gives the total cross-section for the elastic setting (in barn).
  //It integrates numerically the differential cross-section. The number of 
  //integration steps can be specified. Default value is 1000.
  //It does not take radiative corrections into account.
  //It doesn't check for electron in calorimeter, etc. For a precise 
  //calculation of the total cross-section of the (generated) run and then
  //the running time, do as in the example at test/macroelas.C
 
  Double_t cs=0.;
  
  Double_t thetapmax=fSpecAngle+fSpecHorAcc;
  Double_t thetapmin=fSpecAngle-fSpecHorAcc;
  Double_t ppmax=fSpecMom*(1+fSpecMomAcc);
  Double_t ppmin=fSpecMom*(1-fSpecMomAcc);

  Double_t Q2minP=2.*fm*fm*(TMath::Sqrt(1+TMath::Power(ppmin/fm,2))-1.);
  Double_t Q2maxP=2.*fm*fm*(TMath::Sqrt(1+TMath::Power(ppmax/fm,2))-1.);
  
  Double_t Q2maxA=-(4.*fm*fm*TMath::Power(TMath::Cos(thetapmin),2.))/(TMath::Power(TMath::Cos(thetapmin),2.)-TMath::Power(1+fm/fEbeam,2.));
  Double_t Q2minA=-(4.*fm*fm*TMath::Power(TMath::Cos(thetapmax),2.))/(TMath::Power(TMath::Cos(thetapmax),2.)-TMath::Power(1+fm/fEbeam,2.));
  
  Double_t Q2min=TMath::Max(Q2minP,Q2minA);
  Double_t Q2max=TMath::Min(Q2maxP,Q2maxA);
  
  Double_t domega=2.*fSpecVerAcc*(Q2max-Q2min)/steps;

  for(Int_t i=0;i<steps;i++){
    Double_t q2=Q2min+i*(Q2max-Q2min)/steps;

    Double_t thetae=TMath::ACos((2.*TMath::Power(fEbeam,2.)*fm-q2*(fm+fEbeam))/(2.*TMath::Power(fEbeam,2.)*fm-q2*fEbeam));
    Double_t nscs=TMath::Power(Alpha()*TMath::Cos(thetae/2),2)/(4.*TMath::Power(fEbeam,2)*TMath::Power(TMath::Sin(thetae/2.),4)*(1+2.*fEbeam*TMath::Power(TMath::Sin(thetae/2.),2)/fm));//non-structured cross-section
    
    Double_t ge=TMath::Power(1+q2/0.71,-2);
    Double_t gm=2.792847386*ge;
    
    Double_t tau=q2/(4.*TMath::Power(fm,2));
    
    Double_t cst=nscs*((TMath::Power(ge,2)+tau*TMath::Power(gm,2))/(1+tau)+2.*tau*TMath::Power(gm,2)*TMath::Power(TMath::Tan(thetae/2.),2));//cross-section
    cst*=0.389379e-3; //cross-section in bn
    cs+=cst*domega;
  }
  return cs;///((Q2max-Q2min)*(2.*fSpecVerAcc));
}

//_____________________________________________________________________________
Double_t TGenElas::XSecMax(void)
{
  //Computes the maximal elastic diferential cross-section for the event

  Double_t eel=feini->E();
  Double_t thetae=TMath::ACos((2.*TMath::Power(eel,2.)*fm-fQ2min*(fm+eel))/(2.*TMath::Power(eel,2.)*fm-fQ2min*eel));
  
  Double_t nscs=TMath::Power(Alpha()*TMath::Cos(thetae/2),2)/(4.*TMath::Power(eel,2)*TMath::Power(TMath::Sin(thetae/2.),4)*(1+2.*eel*TMath::Power(TMath::Sin(thetae/2.),2)/fm));//non-structured cross-section
  
  Double_t ge=TMath::Power(1+fQ2min/0.71,-2);
  Double_t gm=2.792847386*ge;
  
  Double_t tau=fQ2min/(4.*TMath::Power(fm,2));
  
  Double_t cs=nscs*((TMath::Power(ge,2)+tau*TMath::Power(gm,2))/(1+tau)+2.*tau*TMath::Power(gm,2)*TMath::Power(TMath::Tan(thetae/2.),2));//cross-section
  cs*=0.389379e-3; //cross-section in bn

  return cs;

}

//_____________________________________________________________________________
void TGenElas::ApplySpecVerAcc(Double_t aav)
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
 TLorentzVector* TGenElas::GetFinalProton(void)
{
  // Returns the final photon 4-vector if it exists
  
  if(!fp) cout<<"Warning : Final proton doesn't exist"<<endl;
  return fp;

}

//_____________________________________________________________________________
void TGenElas::Write2File(void)
{

  if(!fNwrite) {
    if(!fFermi) {
      fNwrite=13;
    }else{
      fNwrite=13+3;
    }
  }

  *output<<fVertex->Pz()<<" "<<feini->Pz()<<" ";
  *output<<feprerad->Px()<<" "<<feprerad->Py()<<" "<<feprerad->Pz()<<" ";
  if(fFermi)
    *output<<fpini->Px()<<" "<<fpini->Py()<<" "<<fpini->Pz()<<" ";
  *output<<fp->Px()<<" "<<fp->Py()<<" "<<fp->Pz()<<" ";
  *output<<fe->Px()<<" "<<fe->Py()<<" "<<fe->Pz()<<" ";
  *output<< XSecArrington() << " " <<fPSF <<endl;

}

//_____________________________________________________________________________
 void TGenElas::Print(char* opt)
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
  if(fp){
    cout<<"p'("<<fp->Px()<<","<<fp->Py()<<","<<fp->Pz()<<","<<fp->E()<<")"<<endl;
  }else{
    cout<<"NO RECOIL PARTICLE DEFINED"<<endl;
  }
  cout<<"======================================="<<endl;
}
