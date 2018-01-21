//
// TGenSBSGeo.h, v1.0, Tue Aug 4 11:13:57
// Author: C. Munoz Camacho
//

#ifndef __TGenSBSGeo__
#define __TGenSBSGeo__

#include "TGenGeo.h"

#ifndef ROOT_TLorentzVector
#include "TLorentzVector.h"
#endif

////////////////////////////////////////////////////////////////////////////////
//
// TGenSBSGeo.h
//
// Generator geometry class
// 
//For SBS Geometery, added 03/03/2015 by Z. Ye 
////////////////////////////////////////////////////////////////////////////////

class TGenSBSGeo : public TGenGeo
{
  protected :

	  //For SBS Geometery, added 03/03/2015
  Double_t    fSBSPolarAccMax;// SBS electron maximum polar acceptance (rad)
  Double_t    fSBSPolarAccMin;// SBS electron minimum polar acceptance (rad)
  Double_t    fSBSAzimAccMax; // SBS electron maximum azimuth acceptance (rad)
  Double_t    fSBSAzimAccMin; // SBS electron minimum azimuth acceptance (rad)
  Double_t    fSBSMomAccMin;// SBS electron maximum momentum acceptance (GeV)
  Double_t    fSBSMomAccMax;// SBS electron maximum momentum acceptance (GeV)

  Double_t    fECalPolarAccMax;// ECal photon maximum polar acceptance (rad)
  Double_t    fECalPolarAccMin;// ECal photon minimum polar acceptance (rad)
  Double_t    fECalAzimAccMax; // ECal photon maximum azimuth acceptance (rad)
  Double_t    fECalAzimAccMin; // ECal photon minimum azimuth acceptance (rad)

  Double_t    fTDISPolarAccMax;// TDIS proton maximum polar acceptance (rad)
  Double_t    fTDISPolarAccMin;// TDIS proton minimum polar acceptance (rad)
  Double_t    fTDISAzimAccMax; // TDIS proton maximum azimuth acceptance (rad)
  Double_t    fTDISAzimAccMin; // TDIS proton minimum azimuth acceptance (rad)
  Double_t    fTDISMomAccMin;// TDIS proton maximum momentum acceptance (GeV)
  Double_t    fTDISMomAccMax;// TDIS proton maximum momentum acceptance (GeV)


  // The following variables are unused now since SBS geometry is fixed. 
  // However, I give the option of future modification
  Double_t    fSBSAngle;
  Double_t    fSBSMom;
  Double_t    fECalAngle;
  Double_t    fECalDist;

  public :

  TGenSBSGeo();
  TGenSBSGeo(const TGenSBSGeo&);
  virtual ~TGenSBSGeo();
	
  virtual void SBSPrint(char* opt="");

  void SetSBSDefaultAcceptances(void);
  void SetSBSDefaultAcceptancesGen(void);
  void SetSBSDefaultGeometry(void);

  void SetSBSGeometry(Double_t SpecAngle, Double_t SpecMom, Double_t CaloAngle, Double_t CaloDist);
  Bool_t HitsSBS(TLorentzVector* e);
  Bool_t HitsTDIS(TLorentzVector* p);
  Bool_t HitsECal(TLorentzVector* g);

  void SetSBSAcceptance(Double_t PolAccMin, Double_t PolAccMax, Double_t AzimAccMin, Double_t AzimAccMax, Double_t MomAccMin, Double_t MomAccMax);
  void SetTDISAcceptance(Double_t PolAccMin, Double_t PolAccMax, Double_t AzimAccMin, Double_t AzimAccMax, Double_t MomAccMin, Double_t MomAccMax);
  void SetECalAcceptance(Double_t PolAccMin, Double_t PolAccMax, Double_t AzimAccMin, Double_t AzimAccMax);
  
  ClassDef(TGenSBSGeo,1) // Generator Geometry Class

}; // End of TGenSBSGeo class definition

#endif
