//
// TGenOme.h, v1.0, Tue Aug 4 11:13:57
// Author: C. Munoz Camacho
//

#ifndef __TGenOme__
#define __TGenOme__

#ifndef __TGenBase__
#include "TGenBase.h"
#endif

#ifndef ROOT_TLorentzVector
#include "TLorentzVector.h"
#endif

////////////////////////////////////////////////////////////////////////////////
//
// TGenOme.h
//
// Ome event generator class
// 
////////////////////////////////////////////////////////////////////////////////

class TGenOme : public TGenBase
{
  private :

    TLorentzVector* fome;  // ome
    TLorentzVector* fpi0;   // pi0 (decay from ome)    
    TLorentzVector* fgam;   // gam (decay from ome) 
    TLorentzVector* fg1;   // Emitted photon 1 (decay from pi0)
    TLorentzVector* fg2;   // Emitted photon 2 (decay from pi0)
    TLorentzVector* fp;    // Recoil particle
    Double_t fomem;        // Ome mass

    //    TLorentzVector* fq;    // Virtual photon

    Double_t    fSigmaP;   // Sigma plus cross-section
    Double_t    fSigmaM;   // Sigma minus cross-section
    Double_t    fPSF;      // Phase space factor of the event
                           // * weight (from phase space gen)

    Double_t ft;           // t for the event
    Double_t fs;           // s of the reaction
    Double_t fphi;         // phi of the event

    Double_t ftmin;         // Minimum t for the event
    Double_t ftmax;         // Maximum t for the event


  public :

  TGenOme(Double_t Ebeam, Int_t TargType, UInt_t seed1=1, UInt_t seed2=2);
  TGenOme(const TGenOme&);
  virtual ~TGenOme();

  virtual void Print(char* opt="");

  void IntRCBef(void);
  void IntRCAft(void);

  void Compute2Body();
  void ComputeOme(void);
  void TwoBodyDecay(Double_t M, Double_t m1, Double_t m2);
  void ApplySpecVerAcc(Double_t aav=-1.);
  void Settmin(Double_t tmin) { ftmin=tmin ; }
  void Settmax(Double_t tmax) { ftmax=tmax ; }
  TLorentzVector* GetFinalPhoton1(void);
  TLorentzVector* GetFinalPhoton2(void);
  TLorentzVector* GetFinalProton(void);
  void Write2File(void);
  Double_t GetOmeMass() {return fomem;}

  ClassDef(TGenOme,1) // Ome Event Generator Class

}; // End of TGenOme class definition

#endif
