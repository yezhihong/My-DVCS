//
// TGenPi0Photo.h, v1.0, Tue Aug 4 11:13:57
// Author: C. Munoz Camacho
//

#ifndef __TGenPi0Photo__
#define __TGenPi0Photo__

#ifndef __TGenBase__
#include "TGenBase.h"
#endif

#ifndef ROOT_TLorentzVector
#include "TLorentzVector.h"
#endif

////////////////////////////////////////////////////////////////////////////////
//
// TGenPi0Photo.h
//
// DVCS Event generator class
// 
////////////////////////////////////////////////////////////////////////////////

class TGenPi0Photo : public TGenBase
{
  private :

    TLorentzVector* fp;    // Recoil particle

    Double_t    fSigmaP;   // Sigma plus cross-section
    Double_t    fSigmaM;   // Sigma minus cross-section
    Double_t    fPSF;      // Phase space factor of the event

  public :

  TGenPi0Photo(Double_t Ebeam, Int_t TargType, UInt_t seed1=1, UInt_t seed2=2);
  TGenPi0Photo(const TGenPi0Photo&);
  virtual ~TGenPi0Photo();

  virtual void Print(char* opt="");
  void IntRCBef(void);
  void IntRCAft(void);
  void ComputeElas(void);
  Double_t GammaE(void);
  Double_t XSec(void);
  void GenKin(void);
  //  void ComputeDVCS(void);
  //  void Compute2Body(Double_t m);
  void ApplySpecVerAcc(Double_t aav=-1.);
  Double_t GetPSF(void) { return fPSF ; }
  TLorentzVector* GetFinalProton(void);
  void Write2File(void);
  ClassDef(TGenPi0Photo,1) // Elastic Event Generator Class

}; // End of TGenPi0Photo class definition

#endif
