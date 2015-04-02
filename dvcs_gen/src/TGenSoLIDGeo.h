//
// TGenSoLIDGeo.h, v1.0, Tue Aug 4 11:13:57
// Author: C. Munoz Camacho
//

#ifndef __TGenSoLIDGeo__
#define __TGenSoLIDGeo__

#include "TGenGeo.h"

#ifndef ROOT_TLorentzVector
#include "TLorentzVector.h"
#endif

////////////////////////////////////////////////////////////////////////////////
//
// TGenSoLIDGeo.h
//
// Generator geometry class
// 
//For SoLID Geometery, added 03/03/2015 by Z. Ye 
////////////////////////////////////////////////////////////////////////////////

class TGenSoLIDGeo : public TGenGeo
{
  protected :

	  //For SoLID Geometery, added 03/03/2015
  Double_t    fSoLIDPolarAccMax;// SoLID electron maximum polar acceptance (rad)
  Double_t    fSoLIDPolarAccMin;// SoLID electron minimum polar acceptance (rad)
  Double_t    fSoLIDAzimAccMax; // SoLID electron maximum azimuth acceptance (rad)
  Double_t    fSoLIDAzimAccMin; // SoLID electron minimum azimuth acceptance (rad)
  Double_t    fSoLIDMomAccMin;// SoLID electron maximum momentum acceptance (GeV)
  Double_t    fSoLIDMomAccMax;// SoLID electron maximum momentum acceptance (GeV)

  Double_t    fSoLIDCaloPolarAccMax;// SoLID-FAEC photon maximum polar acceptance (rad)
  Double_t    fSoLIDCaloPolarAccMin;// SoLID-FAEC photon minimum polar acceptance (rad)
  Double_t    fSoLIDCaloAzimAccMax; // SoLID-FAEC photon maximum azimuth acceptance (rad)
  Double_t    fSoLIDCaloAzimAccMin; // SoLID-FAEC photon minimum azimuth acceptance (rad)

  // The following variables are unused now since SoLID geometry is fixed. 
  // However, I give the option of future modification
  Double_t    fSoLIDAngle;
  Double_t    fSoLIDMom;
  Double_t    fSoLIDCaloAngle;
  Double_t    fSoLIDCaloDist;

  public :

  TGenSoLIDGeo();
  TGenSoLIDGeo(const TGenSoLIDGeo&);
  virtual ~TGenSoLIDGeo();
	
  virtual void SoLIDPrint(char* opt="");

  void SetSoLIDDefaultAcceptances(void);
  void SetSoLIDDefaultAcceptancesGen(void);
  void SetSoLIDDefaultGeometry(void);

  void SetSoLIDGeometry(Double_t SpecAngle, Double_t SpecMom, Double_t CaloAngle, Double_t CaloDist);
  Bool_t HitsSoLIDCalo(TLorentzVector* g);
  Bool_t HitsSoLID(TLorentzVector* e);

  void SetSoLIDAcceptance(Double_t PolAccMin, Double_t PolAccMax, Double_t AzimAccMin, Double_t AzimAccMax, Double_t MomAccMin, Double_t MomAccMax);
  void SetSoLIDCaloAcceptance(Double_t PolAccMin, Double_t PolAccMax, Double_t AzimAccMin, Double_t AzimAccMax);
  
  ClassDef(TGenSoLIDGeo,1) // Generator Geometry Class

}; // End of TGenSoLIDGeo class definition

#endif
