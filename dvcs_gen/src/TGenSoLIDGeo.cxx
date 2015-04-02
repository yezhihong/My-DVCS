//
// TGenSoLIDGeo.cxx, v1.0, Tue Aug 4 11:13:57
// Author: C. Munoz Camacho
//

#include <fstream>
#include <iostream>
#include <stdlib.h>

#ifndef __TGenSoLIDGeo__
#include "TGenSoLIDGeo.h"
#endif

#ifndef ROOT_TMath
#include "TMath.h"
#endif

using namespace std;

ClassImp(TGenSoLIDGeo)

////////////////////////////////////////////////////////////////////////////////
//
// Event generator base class
// 
////////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________
 TGenSoLIDGeo::TGenSoLIDGeo()
{
  // Default constructor
  // Default detectors acceptances are set.
  SetSoLIDDefaultGeometry();
  SetSoLIDDefaultAcceptances();
  SetSoLIDDefaultAcceptancesGen();
  fTargDens=0.;

  cout<<"------------------------"<<endl;
  cout<<"For SoLID only: The geometery, acceptance for electron and photon are initialized with default values!"<<endl;
  cout<<"------------------------"<<endl;
}

//_____________________________________________________________________________
 TGenSoLIDGeo::TGenSoLIDGeo(const TGenSoLIDGeo& TCalobase)
{
  // Copy constructor
  //  ((TGenSoLIDGeo&)TCalobase).Copy(*this);
}

//_____________________________________________________________________________
 TGenSoLIDGeo::~TGenSoLIDGeo()
{
  // Default destructor
}

//_____________________________________________________________________________
 void TGenSoLIDGeo::SetSoLIDDefaultGeometry(void)
{
  fSoLIDAngle=0;//rad
  fSoLIDMom=0;//GeV
  fSoLIDCaloAngle=0;//rad
  fSoLIDCaloDist=790.0;//cm
}


//_____________________________________________________________________________
 void TGenSoLIDGeo::SetSoLIDGeometry(Double_t SpecAngle, Double_t SpecMom, Double_t CaloAngle, Double_t CaloDist)
{
  // Sets the angles of the spectro and calorimeter, the central 
  // momentum of the spectrometer and the calorimeter front face distance 
  // from the center of the target

  fSoLIDAngle=SpecAngle;//rad
  fSoLIDMom=SpecMom;//GeV
  fSoLIDCaloAngle=CaloAngle;//rad
  fSoLIDCaloDist=CaloDist;//cm

}


//_____________________________________________________________________________
 void TGenSoLIDGeo::SetSoLIDAcceptance(Double_t PolAccMin, Double_t PolAccMax, Double_t AzimAccMin, Double_t AzimAccMax, Double_t MomAccMin, Double_t MomAccMax)
{
  // Sets SoLID polar (max and min) and azimuthal (max and min)
  // electron acceptances
  // -- Z. Ye, 03/03/2015

  fSoLIDPolarAccMax=PolAccMax;//rad
  fSoLIDPolarAccMin=PolAccMin;//rad
  fSoLIDAzimAccMax=AzimAccMax;//rad
  fSoLIDAzimAccMin=AzimAccMin;//rad
  fSoLIDMomAccMax=MomAccMax;//GeV
  fSoLIDMomAccMin=MomAccMin;//GeV
}

//_____________________________________________________________________________
 void TGenSoLIDGeo::SetSoLIDCaloAcceptance(Double_t PolAccMin, Double_t PolAccMax, Double_t AzimAccMin, Double_t AzimAccMax)
{
  // Sets SoLID calorimeter polar (max and min) and azimuthal (max and min)
  // photon acceptances
  // -- Z. Ye, 03/03/2015

  fSoLIDCaloPolarAccMax=PolAccMax;//rad
  fSoLIDCaloPolarAccMin=PolAccMin;//rad
  fSoLIDCaloAzimAccMax=AzimAccMax;//rad
  fSoLIDCaloAzimAccMin=AzimAccMin;//rad
}


//_____________________________________________________________________________
  void TGenSoLIDGeo::SetSoLIDDefaultAcceptances(void)
{
  //Set SoLID default Acceptance:
  // -- Z. Ye, 03/06/2015
  SetSoLIDAcceptance(6*TMath::DegToRad(), 26.*TMath::DegToRad(),0.*TMath::DegToRad(), 360*TMath::DegToRad(),0.5,11.);//rad,rad,rad,rad,GeV,GeV
  SetSoLIDCaloAcceptance(6*TMath::DegToRad(), 26.*TMath::DegToRad(),0.*TMath::DegToRad(), 360*TMath::DegToRad());
}

//_____________________________________________________________________________
  void TGenSoLIDGeo::SetSoLIDDefaultAcceptancesGen(void)
{
  //Set SoLID default Acceptance:
  // -- Z. Ye, 03/06/2015
  SetSoLIDAcceptance(0.1*TMath::DegToRad(), 30.*TMath::DegToRad(),0.*TMath::DegToRad(), 360*TMath::DegToRad(),0.0,11.0);//rad,rad,rad,rad,GeV,GeV
  SetSoLIDCaloAcceptance(0.1*TMath::DegToRad(), 18.*TMath::DegToRad(),0.*TMath::DegToRad(), 360*TMath::DegToRad());
}


//_____________________________________________________________________________
 Bool_t TGenSoLIDGeo::HitsSoLID(TLorentzVector* p) 
{ 
  // Checks if a particle falls into the defined SoLID acceptances.
  // -- Z. Ye, 03/06/2015

  TVector3 vp=p->Vect();
  TVector3 Oz(0.,0.,1.);

  if(vp.Angle(Oz)>fSoLIDPolarAccMin && vp.Angle(Oz)<fSoLIDPolarAccMax 
                 && p->P()>fSoLIDMomAccMin && p->P()<fSoLIDMomAccMax){
    return kTRUE;
  }else{
    return kFALSE;
  }
}

//_____________________________________________________________________________
 Bool_t TGenSoLIDGeo::HitsSoLIDCalo(TLorentzVector* p) 
{ 
  // Checks if a photon falls into the defined SoLID Calo acceptances.
  // -- Z. Ye, 03/06/2015

  TVector3 vp=p->Vect();
  TVector3 Oz(0.,0.,1.);

  if(vp.Angle(Oz)>fSoLIDCaloPolarAccMin && vp.Angle(Oz)<fSoLIDCaloPolarAccMax){ 
    return kTRUE;
  }else{
    return kFALSE;
  }
}


//_____________________________________________________________________________
 void TGenSoLIDGeo::SoLIDPrint(char* opt)
{
  // Output on screen

  cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
  cout<<"TGenSoLIDGeo--Additional Geometry Settings for SoLID"<<endl;
  cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;

  cout<<"----"<<endl;
  cout<<"SoLID Acceptance for electron: "<<fSoLIDPolarAccMax*180./TMath::Pi()<<" deg pol max "<<fSoLIDPolarAccMin*180./TMath::Pi()<<" deg pol min "<<fSoLIDAzimAccMax*180./TMath::Pi()<<" deg azi max "<<fSoLIDAzimAccMin*180./TMath::Pi()<<" deg azi min" <<fSoLIDMomAccMax <<" GeV mom max "<<fSoLIDMomAccMin<<" GeV mom min"<<endl;
  cout<<"SoLID Calorimeter Acceptance for photon: "<<fSoLIDCaloPolarAccMax*180./TMath::Pi()<<" deg pol max "<<fSoLIDCaloPolarAccMin*180./TMath::Pi()<<" deg pol min "<<fSoLIDCaloAzimAccMax*180./TMath::Pi()<<" deg azi max "<<fSoLIDCaloAzimAccMin*180./TMath::Pi()<<" deg azi min"<<endl;
  cout<<"================================================================="<<endl;
}
