//
// TGenElas.h, v1.0, Tue Aug 4 11:13:57
// Author: C. Munoz Camacho
//

#ifndef __TGenElas__
#define __TGenElas__

#ifndef __TGenBase__
#include "TGenBase.h"
#endif

#ifndef ROOT_TLorentzVector
#include "TLorentzVector.h"
#endif

////////////////////////////////////////////////////////////////////////////////
//
// TGenElas.h
//
// DVCS Event generator class
// 
////////////////////////////////////////////////////////////////////////////////

class TGenElas : public TGenBase
{
  private :

    TLorentzVector* fp;    // Recoil particle

    Double_t    fSigmaP;   // Sigma plus cross-section
    Double_t    fSigmaM;   // Sigma minus cross-section
    Double_t    fPSF;      // Phase space factor of the event

  public :

  TGenElas(Double_t Ebeam, Int_t TargType, UInt_t seed1=1, UInt_t seed2=2);
  TGenElas(const TGenElas&);
  virtual ~TGenElas();

  virtual void Print(char* opt="");
  void IntRCBef(void);
  void IntRCAft(void);
  void ComputeElas(void);
  Double_t KellyE(Double_t);
  Double_t KellyM(Double_t);
  Double_t ArringE(Double_t);
  Double_t ArringM(Double_t);
  Double_t XSec(void);
  Double_t XSecKelly(void);
  Double_t XSecArrington(void);
  Double_t XSecTotal(Int_t steps=1000);
  Double_t XSecMax(void);
  void GenKin(void);
  //  void ComputeDVCS(void);
  //  void Compute2Body(Double_t m);
  void ApplySpecVerAcc(Double_t aav=-1.);
  Double_t GetPSF(void) { return fPSF ; }
  TLorentzVector* GetFinalProton(void);
  void Write2File(void);
  ClassDef(TGenElas,1) // Elastic Event Generator Class

}; // End of TGenElas class definition

#endif
