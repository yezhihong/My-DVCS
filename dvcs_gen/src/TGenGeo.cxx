//
// TGenGeo.cxx, v1.0, Tue Aug 4 11:13:57
// Author: C. Munoz Camacho
//

#include <fstream>
#include <iostream>
#include <stdlib.h>

#ifndef __TGenGeo__
#include "TGenGeo.h"
#endif

#ifndef ROOT_TMath
#include "TMath.h"
#endif

using namespace std;

ClassImp(TGenGeo)

////////////////////////////////////////////////////////////////////////////////
//
// Event generator base class
// 
////////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________
 TGenGeo::TGenGeo()
{
  // Default constructor
  // Default detectors acceptances are set.

  SetDefaultAcceptances();
  SetDefaultAcceptancesGen();

  fSpecAngle=0.;
  fSpecMom=0.;
  fCaloAngle=0.;
  fCaloDist=0.;

  fTargDens=0.;

  cout<<"WARNING: Everything is initialized with default values, except"<<endl;
  cout<<"the target parameters and the configuration of detectors"<<endl;
  cout<<"Please use SetTargetParam and SetGeometry to do this before continuing..."<<endl;
  cout<<"------------------------"<<endl;
}

//_____________________________________________________________________________
 TGenGeo::TGenGeo(const TGenGeo& TCalobase)
{
  // Copy constructor
  //  ((TGenGeo&)TCalobase).Copy(*this);
}

//_____________________________________________________________________________
 TGenGeo::~TGenGeo()
{
  // Default destructor
}

//_____________________________________________________________________________
 void TGenGeo::SetGeometry(Double_t SpecAngle, Double_t SpecMom, Double_t CaloAngle, Double_t CaloDist)
{
  // Sets the angles of the spectro and calorimeter (+PA), the central 
  // momentum of the spectrometer and the calorimeter front face distance 
  // from the center of the target

  fSpecAngle=SpecAngle;
  fSpecMom=SpecMom;
  fCaloAngle=CaloAngle;
  fCaloDist=CaloDist;

}

//_____________________________________________________________________________
 void TGenGeo::SetTargetParam(Double_t targlength, Double_t targden)
{
  // Sets target parameters : length and density

  fTargLength=targlength;//cm
  fTargDens=targden;//g/cm^3
  fTargZoff=0.;

}
//_____________________________________________________________________________
 void TGenGeo::SetTargetParam(Double_t targlength,Double_t targzoff, Double_t targden)
{
  // Sets target parameters : length and density

  fTargLength=targlength;
  fTargDens=targden;
  fTargZoff=targzoff;

}

//_____________________________________________________________________________
 void TGenGeo::SetSpectroAcceptance(Double_t HorAcc, Double_t VerAcc, Double_t MomAcc)
{
  // Sets spectrometer horizontal and vertical acceptance (in rad) 
  // and momentum acceptance: (delta p)/p

  fSpecHorAcc=HorAcc;
  fSpecVerAcc=VerAcc;
  fSpecMomAcc=MomAcc;

}
//_____________________________________________________________________________
 void TGenGeo::SetSpectroAcceptanceGen(Double_t HorAcc, Double_t VerAcc, Double_t MomAcc)
{
  // Sets spectrometer horizontal and vertical acceptance (in rad) 
  // and momentum acceptance: (delta p)/p

  fSpecHorAccGen=HorAcc;
  fSpecVerAccGen=VerAcc;
  fSpecMomAccGen=MomAcc;

}

//_____________________________________________________________________________
 void TGenGeo::SetCaloAcceptance(Double_t HorAccR, Double_t HorAccL, Double_t VerAccU, Double_t VerAccD)
{
  // Sets calorimeter horizontal (right and left) and vertical (up and down)
  // acceptances. All distances are in mm.

  fCaloHorAccR=HorAccR;
  fCaloHorAccL=HorAccL;
  fCaloVerAccU=VerAccU;
  fCaloVerAccD=VerAccD;

}

//_____________________________________________________________________________
 void TGenGeo::SetPAAcceptance(Double_t PolAccMax, Double_t PolAccMin, Double_t AzimAccMax, Double_t AzimAccMin)
{
  // Sets proton array polar (max and min) and azimuthal (max and min)
  // acceptances

  fPAPolarAccMax=PolAccMax;
  fPAPolarAccMin=PolAccMin;
  fPAAzimAccMax=AzimAccMax;
  fPAAzimAccMin=AzimAccMin;

}

//_____________________________________________________________________________
 void TGenGeo::SetDefaultAcceptances(void)
{
  // Sets default acceptances for the spectrometer, calorimeter and
  // proton array
  
  SetSpectroAcceptance(0.08,0.15,0.06);
  //  SetCaloAcceptance(150.,180.,180.,180.);  // calorimeter is not centered
  SetCaloAcceptance(250.,250.,300.,300.);  // calorimeter is not centered
                                           // in final version
  SetPAAcceptance(0.66323,0.31416,5.4978,0.78540);

}

//_____________________________________________________________________________
 void TGenGeo::SetDefaultAcceptancesGen(void)
{
  // Sets default acceptances for the spectrometer, calorimeter and
  // proton array
  
  SetSpectroAcceptanceGen(0.1,0.2,0.07);

}

//_____________________________________________________________________________
 Double_t TGenGeo::PX0(void) 
{ 
  // Returns the radiation length of the LH2 target

  if(fTargDens==0.) cout<<"Target not initialized !"<<endl; 
  return 61.28/fTargDens ; 
}

//_____________________________________________________________________________
 Double_t TGenGeo::NX0(void) 
{ 
  // Returns the radiation length of the LD2 target

  if(fTargDens==0.) cout<<"Target not initialized !"<<endl; 
  return 122.4/fTargDens ; 
}

 Double_t TGenGeo::He3X0(void) 
{ 
  // Returns the radiation length of the LD2 target

  if(fTargDens==0.) cout<<"Target not initialized !"<<endl; 
  return 67.42/fTargDens ; 
}

//_____________________________________________________________________________
 Bool_t TGenGeo::HitsSpectro(TLorentzVector* e) 
{ 
  // Checks if a particle falls into the defined spectrometer acceptances.
  // Assumes particle in the horizontal plane, 

  Double_t thetaemax=fSpecAngle+fSpecHorAcc;
  Double_t thetaemin=fSpecAngle-fSpecHorAcc;
  Double_t phiemax=fSpecVerAcc;
  Double_t phiemin=-fSpecVerAcc;
  Double_t pemax=fSpecMom*(1+fSpecMomAcc);
  Double_t pemin=fSpecMom*(1-fSpecMomAcc);
  
  Double_t Pe=e->P();
  Double_t Thetae=TMath::ATan2(e->Px(),e->Pz());
  Double_t Phie=TMath::ATan2(e->Py(),Pe);

  if(Pe>pemin && Pe<pemax && Thetae>thetaemin && Thetae<thetaemax &&
     Phie>-fSpecVerAcc && Phie<fSpecVerAcc) {
    return kTRUE;
  }else{
    return kFALSE;
  }
}
//_____________________________________________________________________________
 Bool_t TGenGeo::HitsCalo(TLorentzVector* g) 
{ 
  // Checks if a particle falls into the defined calorimeter acceptances.

  g->RotateY(fCaloAngle);

  Double_t ahl=-TMath::ATan(fCaloHorAccL/fCaloDist);
  Double_t ahr= TMath::ATan(fCaloHorAccR/fCaloDist);
  Double_t avu= TMath::ATan(fCaloVerAccU/fCaloDist);
  Double_t avd=-TMath::ATan(fCaloVerAccD/fCaloDist);

  if(TMath::ASin(g->Py()/g->E())<avu && TMath::ASin(g->Py()/g->E())>avd &&
     TMath::ASin(g->Px()/g->E())<ahr && TMath::ASin(g->Px()/g->E())>ahl){
    g->RotateY(-fCaloAngle);//We put it back where it was!
    return kTRUE;
  }else{
    g->RotateY(-fCaloAngle);//We put it back where it was!
    return kFALSE;
  }
}

//_____________________________________________________________________________
 Bool_t TGenGeo::HitsPA(TLorentzVector* p) 
{ 
  // Checks if a particle falls into the defined proton array acceptances.

  p->RotateY(fCaloAngle);

  TVector3 vp=p->Vect();
  TVector3 Oz(0.,0.,1.);

  if(vp.Angle(Oz)>fPAPolarAccMin && vp.Angle(Oz)<fPAPolarAccMax && (p->Px()<0.
  || (p->Px()>0. && (p->Py()>(p->Perp())*TMath::Sin(fPAAzimAccMin)
  ||                 p->Py()<(p->Perp())*TMath::Sin(fPAAzimAccMax))))){
    p->RotateY(-fCaloAngle);//We put it back where it was!
    return kTRUE;
  }else{
    p->RotateY(-fCaloAngle);//We put it back where it was!
    return kFALSE;
  }
}

//_____________________________________________________________________________
 void TGenGeo::Print(char* opt)
{
  // Output on screen

  cout<<"================================================================="<<endl;
  cout<<"TGenGeo"<<endl;
  cout<<"================================================================="<<endl;
  cout<<"Target type : ";
  if(fTargType!=0 && fTargType!=1){
    cout<<" UNKOWN !"<<endl;
  }else{
    if(fTargType==0) {
      cout<<"LH2"<<endl;
    }else{
      cout<<"LD2"<<endl;
    }
  }
  cout<<"Target length and density : "<<fTargLength<<" cm and ";
  cout<<fTargDens<<" g/cm^3"<<endl;
  cout<<"Target Z offset : "<<fTargZoff<<" cm." << endl;
  cout<<"-----------------GEOMETRY-----------------"<<endl;
  cout<<"Spectrometer : "<<fSpecAngle*180./TMath::Pi()<<" deg and ";
  cout<<fSpecMom<<" GeV"<<endl;
  cout<<"Calorimeter : "<<fCaloAngle*180./TMath::Pi()<<" deg and at ";
  cout<<fCaloDist<<" cm"<<endl;
  cout<<"----"<<endl;
  cout<<"Spectro Acceptance : "<<fSpecHorAcc*1000.<<" mrad horizontal "<<fSpecVerAcc*1000.<<" mrad vertical "<<endl;
  cout<<"Spectro Momentum Acceptance : "<<fSpecMomAcc*100<<"%"<<endl;
  cout<<"Calo Acceptance : "<<fCaloHorAccR*180./TMath::Pi()<<" deg hor R "<<fCaloHorAccL*180./TMath::Pi()<<" deg hor L "<<fCaloVerAccU*180./TMath::Pi()<<" deg ver U "<<fCaloVerAccD*180./TMath::Pi()<<" deg ver D"<<endl;
  cout<<"PA Acceptance : "<<fPAPolarAccMax*180./TMath::Pi()<<" deg pol max "<<fPAPolarAccMin*180./TMath::Pi()<<" deg pol min "<<fPAAzimAccMax*180./TMath::Pi()<<" deg azi max "<<fPAAzimAccMin*180./TMath::Pi()<<" deg azi min"<<endl;
  cout<<"================================================================="<<endl;
}
