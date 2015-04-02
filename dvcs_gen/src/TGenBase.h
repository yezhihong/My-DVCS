//
// TGenBase.h, v1.0, Tue Aug 4 11:13:57
// Author: C. Munoz Camacho
//

#ifndef __TGenBase__
#define __TGenBase__

#ifndef __TGenSoLIDGeo__
#include "TGenSoLIDGeo.h"
#endif

#ifndef ROOT_TVector3
#include "TVector3.h"
#endif

#ifndef ROOT_TRandom2
#include "TRandom2.h"
#endif

////////////////////////////////////////////////////////////////////////////////
//
// TGenBase.h
//
// Event generator base class
// 
////////////////////////////////////////////////////////////////////////////////

class TGenBase : public TGenSoLIDGeo
{
  protected :

  static Bool_t           fgIsInit;            // Is TGenBase initialized ?
  static Bool_t           fgWarnings;          // Display warnings ?
  static Bool_t           fgErrors;            // Display errors ?

  TLorentzVector* feini; // Initial electron
  TLorentzVector* fpini; // Initial proton (or neutron)
  TLorentzVector* fe;    // Scattered electron
  TLorentzVector* feprerad;    // Scattered electron at vertex
  TLorentzVector* fq;    // Virtual photon

  Double_t                fEbeam;              // Beam energy

    Double_t fQ2;          // Q2 for the event
    Double_t fxb;          // xb for the event
    Double_t fQ2min;        // Minimum Q2 for the event
    Double_t fQ2max;        // Maximum Q2 for the event
    Double_t fxbmin;        // Minimum x_b for the event
    Double_t fxbmax;        // Maximum x_b for the event

  Bool_t     fFermi;        // Fermi momentum corrections (?)
  Bool_t     fRadCor;       // External bremsstrahlung before the vertex and internal (elastic) radiative corrections
  
  TVector3*  fVertex;       // Vertex of interaction

  TRandom2*  fRan;          // Random number generator

  Double_t*  fdmom;         //! Deuteron momentum distribution;
  Double_t*  frho;          //! Deuteron density (?) distribution;

  ofstream* output; //! Temporary output file
  Int_t fNwrite; //! Number of parameters per event to write

  public :

  TGenBase(Double_t Ebeam, UInt_t seed1=1, UInt_t seed2=2);
  TGenBase(const TGenBase&);
  virtual ~TGenBase();

  virtual void Init(void);

  virtual void Print(char* opt="");

  Double_t GetBeamEnergy (void) { return fEbeam ; }
  void SetBeamEnergy (Double_t ebeam) { fEbeam=ebeam ; }

  void SetRadCor (Bool_t rc = kTRUE) { fRadCor=rc ; }
  void SetFermi (Bool_t fm = kTRUE) { fFermi=fm ; }

  TRandom2* GetRanGen(void) { return fRan ; }

  TLorentzVector* GetScatteredElectron(void);
  TLorentzVector* GetInitialElectron(void);
  void GenKin(void);
  void GenKinGen(void);
  void GenKinCLAS(void);
  void GenKinSoLID(void);
  TVector3* GetVertex(void);
  Int_t ComputeElectron(void);
  void ExtBrem(void);
  void SetWarnings(Bool_t val) { fgWarnings = val; }
  void SetErrors(Bool_t val) { fgErrors = val; }
  void SetVertex(Double_t x, Double_t y, Double_t z);
  void SetVertex(TVector3* vertex);
  TLorentzVector* GenFermiIni(void);
  void SetQ2(Double_t q2) { fQ2=q2 ; }
  void SetXb(Double_t xb) { fxb=xb ; }
  Double_t GetQ2(void) { return fQ2 ; }
  Double_t GetXb(void) { return fxb ; }
  void GenerateVertex(void);
  void ApplySpecVerAcc(Double_t aav);
  void CloseTmpFile();
  void DumpFinalFile(char* finalfile, Int_t Ngen, Int_t Nacc);
  void DumpRootFile(char* finalfile, Int_t Ngen, Int_t Nacc);
  void DumpRootFilePi0(char* finalfile, Int_t Ngen, Int_t Nacc);
  void Clear(char* opt="");

  ClassDef(TGenBase,1) // Event Generator Base Class

}; // End of TGenBase class definition

#endif
