//
// TGenPi0.h, v1.0, Tue Aug 4 11:13:57
// Author: C. Munoz Camacho
//

#ifndef __TGenPi0__
#define __TGenPi0__

#ifndef __TGenBase__
#include "TGenBase.h"
#endif

#ifndef ROOT_TLorentzVector
#include "TLorentzVector.h"
#endif

////////////////////////////////////////////////////////////////////////////////
//
// TGenPip.h
//
// Pip event generator class
// 
////////////////////////////////////////////////////////////////////////////////

class TGenPip : public TGenBase
{
  private :

    TLorentzVector* fpip;  // \pi_0
    TLorentzVector* fn;    // Recoil particle

    //    TLorentzVector* fq;    // Virtual photon

    Double_t    fSigmaP;   // Sigma plus cross-section
    Double_t    fSigmaM;   // Sigma minus cross-section
    Double_t    fPSF;      // Phase space factor of the event

    Double_t ft;           // t for the event
    Double_t fs;           // s of the reaction
    Double_t fphi;         // phi of the event

    Double_t ftmin;         // Minimum t for the event
    Double_t ftmax;         // Maximum t for the event

    Double_t ft1;         // Minimum t for kinematics
    Double_t ft0;         // Maximum t for kinematics

    Double_t fpcm1;       // virtual photon
    Double_t fpcm2;       // proton
    Double_t fpcm3;       // pi0
    Double_t fpcm4;       // delta
 
    Double_t fecm1;
    Double_t fecm3;

  public :

  TGenPip(Double_t Ebeam, Int_t TargType, UInt_t seed1=1, UInt_t seed2=2);
  TGenPip(const TGenPip&);
  virtual ~TGenPip();

  virtual void Print(char* opt="");

  void IntRCBef(void);
  void IntRCAft(void);

  void Compute2Body(Double_t m);
  void ComputePip(void);
  void ComputeTminTmax(Double_t m3,Double_t m4);
  void ApplySpecVerAcc(Double_t aav=-1.);
  void Settmin(Double_t tmin) { ftmin=tmin ; }
  void Settmax(Double_t tmax) { ftmax=tmax ; }
  TLorentzVector* GetFinalPip(void);
  TLorentzVector* GetFinalNeutron(void);
  void Write2File(void);

  ClassDef(TGenPip,1) // Pip Event Generator Class

}; // End of TGenPip class definition

#endif
