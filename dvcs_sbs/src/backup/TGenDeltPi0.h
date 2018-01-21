//
// TGenDeltPi0.h, v1.0, Tue Aug 4 11:13:57
// Author: C. Munoz Camacho
//

#ifndef __TGenDeltPi0__
#define __TGenDeltPi0__

#ifndef __TGenBase__
#include "TGenBase.h"
#endif

#ifndef ROOT_TLorentzVector
#include "TLorentzVector.h"
#endif

////////////////////////////////////////////////////////////////////////////////
//
// TGenDeltPi0.h
//
// DeltPi0 event generator class
// 
////////////////////////////////////////////////////////////////////////////////

class TGenDeltPi0 : public TGenBase
{
  private :

    TLorentzVector* fpi0;  // \pi_0
    TLorentzVector* fg1;   // Emitted photon 1    
    TLorentzVector* fg2;   // Emitted photon 2   
    TLorentzVector* fp;    // Recoil particle
    Double_t fdeltam;        // Delta mass

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




  public :

  TGenDeltPi0(Double_t Ebeam, Int_t TargType, UInt_t seed1=1, UInt_t seed2=2);
  TGenDeltPi0(const TGenDeltPi0&);
  virtual ~TGenDeltPi0();

  virtual void Print(char* opt="");

  void IntRCBef(void);
  void IntRCAft(void);

  void Compute2Body(Double_t m);
  void ComputePi0(void);
  void ComputeTminTmax(Double_t m3,Double_t m4);
  void TwoBodyDecay(Double_t M, Double_t m1, Double_t m2);
  void ApplySpecVerAcc(Double_t aav=-1.);
  void Settmin(Double_t tmin) { ftmin=tmin ; }
  void Settmax(Double_t tmax) { ftmax=tmax ; }
  TLorentzVector* GetFinalPhoton1(void);
  TLorentzVector* GetFinalPhoton2(void);
  TLorentzVector* GetFinalProton(void);
  void Write2File(void);

  ClassDef(TGenDeltPi0,1) // DeltPi0 Event Generator Class

}; // End of TGenDeltPi0 class definition

#endif
