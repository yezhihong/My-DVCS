//
// TGenGeo.h, v1.0, Tue Aug 4 11:13:57
// Author: C. Munoz Camacho
//

#ifndef __TGenGeo__
#define __TGenGeo__

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef ROOT_TLorentzVector
#include "TLorentzVector.h"
#endif

////////////////////////////////////////////////////////////////////////////////
//
// TGenGeo.h
//
// Generator geometry class
// 
////////////////////////////////////////////////////////////////////////////////

class TGenGeo : public TObject
{
  protected :

  Int_t                   fTargType;           // Target type
  Double_t                fTargLength;         // Target length
  Double_t                fTargDens;           // Target density
  Double_t                fTargZoff;           // Target Z offset in cm
  Double_t                fm;                  // Mass of the recoil particle

  Double_t    fSpecAngle;          // Left spectrometer angle (rad)
  Double_t    fSpecMom;            // Central spectrometer momemtum (GeV)
  Double_t    fCaloAngle;          // Calorimeter+PA angle (rad)
  Double_t    fCaloDist;     // Calorimeter front face distance from the center of the target

  Double_t    fSpecHorAcc;   // Spectrometer horizontal acceptance (rad)
  Double_t    fSpecVerAcc;   // Spectrometer vertical acceptance (rad)
  Double_t    fSpecMomAcc;   // Spectrometer momentum acceptance: (delta p)/p
  Double_t    fSpecHorAccGen;   // Spectrometer horizontal acceptance (rad)
  Double_t    fSpecVerAccGen;   // Spectrometer vertical acceptance (rad)
  Double_t    fSpecMomAccGen;   // Spectrometer momentum acceptance: (delta p)/p
  Double_t    fCaloHorAccR;  // Calorimeter horizontal right acceptance (mm)
  Double_t    fCaloHorAccL;  // Calorimeter horizontal left acceptance (mm)
  Double_t    fCaloVerAccU;  // Calorimeter vertical upwards acceptance (mm)
  Double_t    fCaloVerAccD;  // Calorimeter vertical downwards acceptance (mm)
  Double_t    fPAPolarAccMax;// PA maximum polar acceptance (rad)
  Double_t    fPAPolarAccMin;// PA minimum polar acceptance (rad)
  Double_t    fPAAzimAccMax; // PA maximum azimuth acceptance (rad)
  Double_t    fPAAzimAccMin; // PA minimum azimuth acceptance (rad)

  public :

  TGenGeo();
  TGenGeo(const TGenGeo&);
  virtual ~TGenGeo();

  virtual void Print(char* opt="");

  static Double_t PMass(void) { return 0.938271998 ; }
  static Double_t NMass(void) { return 0.93956533 ; }
  static Double_t DMass(void) { return 1.875613 ; }
  static Double_t He3Mass(void) { return 2.809 ; }
  static Double_t b(void) { return 1.35742 ; }
  static Double_t Alpha(void) {return 7.297352533e-3 ; }

  Double_t PX0(void);
  Double_t NX0(void);
  Double_t He3X0(void);

  void SetTargetParam(Double_t targlength, Double_t targzoffset, Double_t targden);
  void SetTargetParam(Double_t targlength, Double_t targden);
  void SetDefaultAcceptances(void);
  void SetDefaultAcceptancesGen(void);
  void SetSpectroAcceptance(Double_t HorAcc, Double_t VerAcc, Double_t MomAcc);
  void SetSpectroAcceptanceGen(Double_t HorAcc, Double_t VerAcc, Double_t MomAcc);
  void SetCaloAcceptance(Double_t HorAccR, Double_t HorAccL, Double_t VerAccU, Double_t VerAccD);
  void SetPAAcceptance(Double_t PolAccMax, Double_t PolAccMin, Double_t AzimAccMax, Double_t AzimAccMin);
  void SetGeometry(Double_t SpecAngle, Double_t SpecMom, Double_t CaloAngle, Double_t CaloDist);

  void SetSpecAngle(Double_t SpecAngle) { fSpecAngle=SpecAngle ; }
  void SetSpecMomentum(Double_t SpecMom) { fSpecMom=SpecMom ; }
  void SetCaloAngle (Double_t CaloAngle) { fCaloAngle=CaloAngle; }
  void SetCaloDistance (Double_t CaloDist) { fCaloDist=CaloDist ; }

  Double_t GetSpecAngle (void) { return fSpecAngle ; }
  Double_t GetSpecMomentum (void) { return fSpecMom ; }
  Double_t GetCaloAngle (void) { return fCaloAngle ; }
  Double_t GetCaloDistance (void) { return fCaloDist ; }

  Double_t GetSpecVerAcc(void) { return fSpecVerAcc ; }

  Bool_t HitsSpectro(TLorentzVector *e);
  Bool_t HitsCalo(TLorentzVector* g);
  Bool_t HitsPA(TLorentzVector* p);

  ClassDef(TGenGeo,1) // Generator Geometry Class

}; // End of TGenGeo class definition

#endif
