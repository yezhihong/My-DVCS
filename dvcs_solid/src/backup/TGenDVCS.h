//
// TGenDVCS.h, v1.0, Tue Aug 4 11:13:57
// Author: C. Munoz Camacho
//

#ifndef __TGenDVCS__
#define __TGenDVCS__

#ifndef __TGenBase__
#include "TGenBase.h"
#endif

#ifndef __TGVKelly__
#include "TGVKelly.h"
#endif

#ifndef ROOT_TLorentzVector
#include "TLorentzVector.h"
#endif

////////////////////////////////////////////////////////////////////////////////
//
// TGenDVCS.h
//
// DVCS Event generator class
// 
////////////////////////////////////////////////////////////////////////////////

class TGenDVCS : public TGenBase
{
  private :

    TLorentzVector* fg;    // Emitted photon    
    TLorentzVector* fp;    // Recoil particle

    Double_t    fSigmaP;   // Sigma plus cross-section
    Double_t    fSigmaM;   // Sigma minus cross-section
    Double_t    fPSF;      // Phase space factor of the event

    Double_t ft;           // t for the event
    Double_t fs;           // s of the reaction
    Double_t fphi;         // phi of the event

    Double_t ftmin;         // Minimum t for the event
    Double_t ftmax;         // Maximum t for the event

    //Theory parameters
    Int_t fphasespace;
    Int_t fdterm;
    Int_t fDD;
    Int_t fpipole;
    Int_t ftdep;
    Int_t fprop;
    Double_t fb;
    Double_t ftcoef;
    Double_t fJu;
    Double_t fJd;
    
    //CFFs (GK for TGV)
    Double_t ****V;//!
    TGVKelly *tgv;
    Double_t *CFF;//!

 public :

  TGenDVCS(Double_t Ebeam, Int_t TargType, UInt_t seed1=1, UInt_t seed2=2);
  TGenDVCS(const TGenDVCS&);
  virtual ~TGenDVCS();

  //   virtual void Init(void); 
  virtual void Print(char* opt="");
  void IntRCBef(void);
  void IntRCAft(void);
  void ComputeDVCS(void);
  void Compute2Body(Double_t m);
  void ApplySpecVerAcc(Double_t aav=-1.);
  Double_t GetFastWeight(void);
  void Settmin(Double_t tmin);
  void Settmax(Double_t tmax);
  void SetTheoryParam(Int_t phasespace, Int_t prop, Double_t b, Int_t tdep, Double_t tcoef, Int_t DD, Double_t Ju, Double_t Jd, Int_t dterm, Int_t pipole);
  void XSec(void);
  TLorentzVector* GetFinalPhoton(void);
  TLorentzVector* GetFinalProton(void);
  void Write2File(void);
  void Clear(char* opt="");


  Double_t GetPSF(){ return fPSF; }
  Double_t KellyE(Double_t q2);
  Double_t KellyM(Double_t q2);
  Double_t* Interpol_CFF(Double_t Q2, Double_t xb, Double_t t);
  Double_t XSecSum(int opt=0);
  Double_t XSecDif(void);
  Double_t GetPhi(void){return fphi;}
  Double_t SetPhi(Double_t phi){fphi=phi;}
 Double_t GetT(void){return ft;}
 Double_t SetT(Double_t t){ft=t;}

  ClassDef(TGenDVCS,1) // DVCS Event Generator Class

}; // End of TGenDVCS class definition

#endif
