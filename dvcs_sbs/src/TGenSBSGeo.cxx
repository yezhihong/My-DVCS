//
// TGenSBSGeo.cxx, v1.0, Tue Aug 4 11:13:57
// Author: C. Munoz Camacho
//

#include <fstream>
#include <iostream>
#include <stdlib.h>

#ifndef __TGenSBSGeo__
#include "TGenSBSGeo.h"
#endif

#ifndef ROOT_TMath
#include "TMath.h"
#endif

using namespace std;

ClassImp(TGenSBSGeo)

////////////////////////////////////////////////////////////////////////////////
//
// Event generator base class
// 
////////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________
 TGenSBSGeo::TGenSBSGeo()
{
  // Default constructor
  // Default detectors acceptances are set.
  SetSBSDefaultGeometry();
  SetSBSDefaultAcceptances();
  SetSBSDefaultAcceptancesGen();
  fTargDens=0.;

  cout<<"------------------------"<<endl;
  cout<<"For SBS only: The geometery, acceptance for electron and photon are initialized with default values!"<<endl;
  cout<<"------------------------"<<endl;
}

//_____________________________________________________________________________
 TGenSBSGeo::TGenSBSGeo(const TGenSBSGeo& TECalbase)
{
  // Copy constructor
  //  ((TGenSBSGeo&)TECalbase).Copy(*this);
}

//_____________________________________________________________________________
 TGenSBSGeo::~TGenSBSGeo()
{
  // Default destructor
}

//_____________________________________________________________________________
 void TGenSBSGeo::SetSBSDefaultGeometry(void)
{
  fSBSAngle=0;//rad
  fSBSMom=0;//GeV
  fECalAngle=0;//rad
  fECalDist=790.0;//cm
}


//_____________________________________________________________________________
 void TGenSBSGeo::SetSBSGeometry(Double_t SpecAngle, Double_t SpecMom, Double_t ECalAngle, Double_t ECalDist)
{
  // Sets the angles of the spectro and calorimeter, the central 
  // momentum of the spectrometer and the calorimeter front face distance 
  // from the center of the target

  fSBSAngle=SpecAngle;//rad
  fSBSMom=SpecMom;//GeV
  fECalAngle=ECalAngle;//rad
  fECalDist=ECalDist;//cm

}

//_____________________________________________________________________________
 void TGenSBSGeo::SetSBSAcceptance(Double_t PolAccMin, Double_t PolAccMax, Double_t AzimAccMin, Double_t AzimAccMax, Double_t MomAccMin, Double_t MomAccMax)
{
  // Sets SBS polar (max and min) and azimuthal (max and min)
  // electron acceptances
  // -- Z. Ye, 01/20/2018

  fSBSPolarAccMax=PolAccMax;//rad
  fSBSPolarAccMin=PolAccMin;//rad
  fSBSAzimAccMax=AzimAccMax;//rad
  fSBSAzimAccMin=AzimAccMin;//rad
  fSBSMomAccMax=MomAccMax;//GeV
  fSBSMomAccMin=MomAccMin;//GeV
}

//_____________________________________________________________________________
 void TGenSBSGeo::SetTDISAcceptance(Double_t PolAccMin, Double_t PolAccMax, Double_t AzimAccMin, Double_t AzimAccMax, Double_t MomAccMin, Double_t MomAccMax)
{
  // Sets TDIS polar (max and min) and azimuthal (max and min)
  // recoil proton acceptances
  // -- Z. Ye, 01/20/2018

  fTDISPolarAccMax=PolAccMax;//rad
  fTDISPolarAccMin=PolAccMin;//rad
  fTDISAzimAccMax=AzimAccMax;//rad
  fTDISAzimAccMin=AzimAccMin;//rad
  fTDISMomAccMax=MomAccMax;//GeV
  fTDISMomAccMin=MomAccMin;//GeV
}

//_____________________________________________________________________________
 void TGenSBSGeo::SetECalAcceptance(Double_t PolAccMin, Double_t PolAccMax, Double_t AzimAccMin, Double_t AzimAccMax)
{
  // Sets SBS calorimeter polar (max and min) and azimuthal (max and min)
  // photon acceptances
  // -- Z. Ye, 01/20/2018

  fECalPolarAccMax=PolAccMax;//rad
  fECalPolarAccMin=PolAccMin;//rad
  fECalAzimAccMax=AzimAccMax;//rad
  fECalAzimAccMin=AzimAccMin;//rad
}



//_____________________________________________________________________________
  void TGenSBSGeo::SetSBSDefaultAcceptances(void)
{
  //Set SBS default Acceptance:
  // -- Z. Ye, 01/20/2018 ,, FIXHERE
  SetSBSAcceptance(6*TMath::DegToRad(), 70.*TMath::DegToRad(),0.*TMath::DegToRad(), 360*TMath::DegToRad(),0.5,11.);//rad,rad,rad,rad,GeV,GeV
  SetECalAcceptance(0*TMath::DegToRad(), 30.*TMath::DegToRad(),0.*TMath::DegToRad(), 360*TMath::DegToRad());
  SetTDISAcceptance(0*TMath::DegToRad(), 180.*TMath::DegToRad(),0.*TMath::DegToRad(), 360*TMath::DegToRad(), 0.0,1.0);
}

//_____________________________________________________________________________
  void TGenSBSGeo::SetSBSDefaultAcceptancesGen(void)
{
  //Set SBS default Acceptance:
  // -- Z. Ye, 01/20/2018 ,, FIXHERE
  SetSBSAcceptance(0.1*TMath::DegToRad(), 70.*TMath::DegToRad(),0.*TMath::DegToRad(), 360*TMath::DegToRad(),0.0,11.0);//rad,rad,rad,rad,GeV,GeV
  SetECalAcceptance(0.1*TMath::DegToRad(), 30.*TMath::DegToRad(),0.*TMath::DegToRad(), 360*TMath::DegToRad());
  SetTDISAcceptance(0.1*TMath::DegToRad(), 180.*TMath::DegToRad(),0.*TMath::DegToRad(), 360*TMath::DegToRad(), 0.0,1.0);
}


//_____________________________________________________________________________
 Bool_t TGenSBSGeo::HitsSBS(TLorentzVector* e) 
{ 
  // Checks if a particle falls into the defined SBS acceptances.
  // -- Z. Ye, 01/20/2018 ,, FIXHERE

  TVector3 ve=e->Vect();
  TVector3 Oz(0.,0.,1.);

  if(ve.Angle(Oz)>fSBSPolarAccMin && ve.Angle(Oz)<fSBSPolarAccMax 
                 && e->P()>fSBSMomAccMin && e->P()<fSBSMomAccMax){
    return kTRUE;
  }else{
    return kFALSE;
  }
}

//_____________________________________________________________________________
 Bool_t TGenSBSGeo::HitsECal(TLorentzVector* g) 
{ 
  // Checks if a photon falls into the defined SBS ECal acceptances.
  // -- Z. Ye, 01/20/2018 ,, FIXHERE

  TVector3 vg=g->Vect();
  TVector3 Oz(0.,0.,1.);

  if(vg.Angle(Oz)>fECalPolarAccMin && vg.Angle(Oz)<fECalPolarAccMax){ 
    return kTRUE;
  }else{
    return kFALSE;
  }
}

//_____________________________________________________________________________
 Bool_t TGenSBSGeo::HitsTDIS(TLorentzVector* p) 
{ 
  // Checks if a photon falls into the defined SBS TDIS acceptances.
  // -- Z. Ye, 01/20/2018 ,, FIXHERE

  TVector3 vp=p->Vect();
  TVector3 Oz(0.,0.,1.);

  if(vp.Angle(Oz)>fTDISPolarAccMin && vp.Angle(Oz)<fTDISPolarAccMax 
                 && p->P()>fTDISMomAccMin && p->P()<fTDISMomAccMax){
    return kTRUE;
  }else{
    return kFALSE;
  }
}


//_____________________________________________________________________________
 void TGenSBSGeo::SBSPrint(char* opt)
{
  // Output on screen

  cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
  cout<<"TGenSBSGeo--Additional Geometry Settings for SBS"<<endl;
  cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;

  cout<<"----"<<endl;
  cout<<"SBS Acceptance for electron: "<<fSBSPolarAccMax*180./TMath::Pi()<<" deg pol max "<<fSBSPolarAccMin*180./TMath::Pi()<<" deg pol min "<<fSBSAzimAccMax*180./TMath::Pi()<<" deg azi max "<<fSBSAzimAccMin*180./TMath::Pi()<<" deg azi min" <<fSBSMomAccMax <<" GeV mom max "<<fSBSMomAccMin<<" GeV mom min"<<endl;
  cout<<"ECalrimeter Acceptance for photon: "<<fECalPolarAccMax*180./TMath::Pi()<<" deg pol max "<<fECalPolarAccMin*180./TMath::Pi()<<" deg pol min "<<fECalAzimAccMax*180./TMath::Pi()<<" deg azi max "<<fECalAzimAccMin*180./TMath::Pi()<<" deg azi min"<<endl;
  cout<<"TDIS Acceptance for proton: "<<fTDISPolarAccMax*180./TMath::Pi()<<" deg pol max "<<fTDISPolarAccMin*180./TMath::Pi()<<" deg pol min "<<fTDISAzimAccMax*180./TMath::Pi()<<" deg azi max "<<fTDISAzimAccMin*180./TMath::Pi()<<" deg azi min"<<endl;
  cout<<"================================================================="<<endl;
}
