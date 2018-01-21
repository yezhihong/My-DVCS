//
// TGenBase.cxx, v1.0, Tue Aug 4 11:13:57
// Author: C. Munoz Camacho
//

#include <fstream>
#include <iostream>
#include <stdlib.h>

#ifndef __TGenBase__
#include "TGenBase.h"
#endif

#ifndef ROOT_TMath
#include "TMath.h"
#endif

/*ROOT Includes{{{*/
#include "TRandom2.h"
#include <TFile.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <Riostream.h>
#include "TObjString.h"
#include <TNamed.h>
#include <TPRegexp.h>
#include <TObjArray.h>
#include <TChain.h>
#include <TF1.h>
#include <TVirtualFitter.h>
//#include <TMatrix.h>
/*}}}*/

using namespace std;

ClassImp(TGenBase)

////////////////////////////////////////////////////////////////////////////////
//
// Event generatorS base class
// 
////////////////////////////////////////////////////////////////////////////////

Bool_t           TGenBase::fgIsInit = kFALSE;
Bool_t           TGenBase::fgWarnings = kTRUE;
Bool_t           TGenBase::fgErrors = kTRUE;

//_____________________________________________________________________________
 TGenBase::TGenBase(Double_t Ebeam, Int_t TargType, UInt_t seed1, UInt_t seed2)
{
  // Default constructor
  // The default behaviour is : seeds 1 and 2, we ask for a particle in 
  // the calorimeter, don't ask for a proton in the proton array,
  // radiative corrections activated, fermi momentum disable.
  // Default detectors acceptances are set.

  fEbeam=Ebeam;
  fTargType=TargType;
  feini=new TLorentzVector(0.,0.,fEbeam,fEbeam);//We neglect electron mass
  feprerad=new TLorentzVector(0.,0.,0.,0.);

  fRan=new TRandom2();
  fRan->SetSeed(0);

  fRadCor=kFALSE;
  fFermi=kFALSE;

  fVertex=0;
  fNwrite=0;
  
  fOutTempFile = Form("temp_dvcs_E%d_T%d.txt", (int) (fEbeam), fTargType);
  cerr<<Form("--- Saving data into a temperate file = %s", fOutTempFile.Data())<<endl;
  output=new ofstream(fOutTempFile);
}

//_____________________________________________________________________________
 TGenBase::TGenBase(const TGenBase& TCalobase)
{
  // Copy constructor
  //  ((TGenBase&)TCalobase).Copy(*this);
}

//_____________________________________________________________________________
 TGenBase::~TGenBase()
{
  // Default destructor
  if (fVertex) delete fVertex;
  delete feini;
  delete feprerad;
  delete fRan;
}

//_____________________________________________________________________________
 void TGenBase::ExtBrem(void)
{
  // Make external bremsstrahlung corrections before the vertex (straggling)
  // Uses \Delta E=E_0 * R^(1/b/t) with R randomly between 0 and 1.

  Double_t t=0; //Lengh of matter in radiation length
  Double_t eel=0.;
  fRadCor=kTRUE;

   if(fTargType!=0 && fTargType!=1 && fTargType!=2 && fTargType!=3){
    cout<<" UNKOWN TARGET!"<<endl;
    exit(1);
  }
  if(fTargType==0){
    t=(fVertex->Pz()+fTargLength/2.)/PX0();
  }else if( fTargType==1 || fTargType==2 ){
    t=(fVertex->Pz()+fTargLength/2.)/NX0();
  }
  else
    {
      t=(fVertex->Pz()+fTargLength/2.)/He3X0();
    }

    Double_t toto=TMath::Power(fRan->Rndm(),1./(b()*t));

    //cout<<"toto="<<toto<<endl;    

    eel=fEbeam*(1.-toto);
    feini->SetPxPyPzE(0.,0.,eel,eel);
}

//_____________________________________________________________________________
 void TGenBase::SetVertex(Double_t x, Double_t y, Double_t z) 
{ 
  // Set Vertex coordinates

  if(!fVertex){
    fVertex=new TVector3(x,y,z);
  }else{
    fVertex->SetXYZ(x,y,z);
  }
}

//_____________________________________________________________________________
 void TGenBase::SetVertex(TVector3* vertex) 
{ 
  // Sets vertex

  if(!fVertex){
    fVertex=new TVector3(vertex->Px(),vertex->Py(),vertex->Pz());
  }else{
    fVertex->SetXYZ(vertex->Px(),vertex->Py(),vertex->Pz());
  }
}

//_____________________________________________________________________________
 void TGenBase::GenerateVertex(void) 
{ 
  // Generates an uniformly distributed vertex along the target length
  // (no rastering), with offset from TGenSBSGeo
  SetVertex(0.,0.,fTargLength*fRan->Rndm()-fTargLength/2.+fTargZoff);
}

//_____________________________________________________________________________
 void TGenBase::GenKin(void)
{
  // Generates the electron event : Q2, and xb (uniformly and randomly between limits)

  // These 4 following values could be set as non-persistent data members of TGenBase
  // and avoid recalculating them for each event.
  Double_t thetaemax=fSpecAngle+fSpecHorAccGen;
  Double_t thetaemin=fSpecAngle-fSpecHorAccGen;
  Double_t pemax=fSpecMom*(1+fSpecMomAccGen);
  Double_t pemin=fSpecMom*(1-fSpecMomAccGen);
 
  //cout<<"th_min="<<thetaemin<<"  th_max="<<thetaemax<<
  //    "pemin="<<pemin<<"pemax="<<pemax<<endl;


  Double_t eel=feini->E();

  fxbmin=TMath::Max(0.05,pemin*eel*(1.-TMath::Cos(thetaemin))/(fm*(eel-pemin)));
  fxbmax=TMath::Min(0.95,pemax*eel*(1.-TMath::Cos(thetaemax))/(fm*(eel-pemax)));
  if(eel<pemin){
    fxbmin=0.00000001;
    fxbmax=0.00000001;
  }//Again, a great energy loss might make our equations invalid

  if (fFermi) {
    fxbmin=TMath::Min(0.1,fxbmin*0.8);
    fxbmax=TMath::Max(0.9,fxbmax*1.2);
  }

  // This is also true for Fermi smeared kinematics

  fQ2min=2.*pemin*eel*(1-TMath::Cos(thetaemin));
  fQ2max=2.*pemax*eel*(1-TMath::Cos(thetaemax)); 
  
  fxb=fxbmin+(fxbmax-fxbmin)*fRan->Rndm();
  fQ2=fQ2min+(fQ2max-fQ2min)*fRan->Rndm();
 
}
//_____________________________________________________________________________
 void TGenBase::GenKinGen(void)
{
  // Generates the electron event : Q2, and xb (uniformly and randomly between limits)

  // These 4 following values could be set as non-persistent data members of TGenBase
  // and avoid recalculating them for each event.
//   Double_t thetaemax=TMath::Abs(fCaloAngle)+0.03;
//   Double_t thetaemin=TMath::Abs(fCaloAngle)-0.03;
  Double_t thetaemax=70.*TMath::DegToRad();
  Double_t thetaemin=0.;
  Double_t pemax=10.;
  Double_t pemin=0.1;
 
  Double_t eel=feini->E();

  fxbmin=TMath::Max(0.05,pemin*eel*(1.-TMath::Cos(thetaemin))/(fm*(eel-pemin)));
  fxbmax=TMath::Min(0.95,pemax*eel*(1.-TMath::Cos(thetaemax))/(fm*(eel-pemax)));
  if(eel<pemin){
    fxbmin=0.00000001;
    fxbmax=0.00000001;
  }//Again, a great energy loss might make our equations invalid

  if (fFermi) {
    fxbmin=TMath::Min(0.1,fxbmin*0.8);
    fxbmax=TMath::Max(0.9,fxbmax*1.2);
  }

  // This is also true for Fermi smeared kinematics

  fQ2min=2.*pemin*eel*(1-TMath::Cos(thetaemin));
  fQ2max=2.*pemax*eel*(1-TMath::Cos(thetaemax)); 
  
  fxb=fxbmin+(fxbmax-fxbmin)*fRan->Rndm();
  fQ2=fQ2min+(fQ2max-fQ2min)*fRan->Rndm();
 
}

//_____________________________________________________________________________
 void TGenBase::GenKinSBS(void)
{
  // Added by Zhihong Ye for SBS config. on 01/20/2018

  // Generates the electron event : Q2, and xb (uniformly and randomly between limits)
  // These 4 following values could be set as non-persistent data members of TGenBase
  // and avoid recalculating them for each event.

//  Double_t eel=feini->E();

  fxbmin=0.1;
  fxbmax=0.8;

  fQ2min=1.0;
  fQ2max=13.0;

  Double_t nuu=-1;
  Double_t ww=-1;

  Bool_t goodkin=kFALSE;

  while(!goodkin){

    fxb=fxbmin+(fxbmax-fxbmin)*fRan->Rndm();
    fQ2=fQ2min+(fQ2max-fQ2min)*fRan->Rndm();

  // Check if it's coherent

    nuu=fQ2/2./fxb/fm;
    ww=fm*fm+2*fm*nuu-fQ2;

    if(nuu<(fEbeam-0.1)) goodkin=kTRUE;

  }


}


//_____________________________________________________________________________
 void TGenBase::GenKinCLAS(void)
{
  // Generates the electron event : Q2, and xb (uniformly and randomly between limits)

  // These 4 following values could be set as non-persistent data members of TGenBase
  // and avoid recalculating them for each event.
//   Double_t thetaemax=fSpecAngle+fSpecHorAccGen;
//   Double_t thetaemin=fSpecAngle-fSpecHorAccGen;
//   Double_t pemax=fSpecMom*(1+fSpecMomAccGen);
//   Double_t pemin=fSpecMom*(1-fSpecMomAccGen);
//   Double_t thetaemax=50.*TMath::DegToRad();
//  Double_t thetaemin=0.;
//  Double_t pemax=4.;
//  Double_t pemin=0.1;

//  Double_t eel=feini->E();

  fxbmin=0.1;//TGKV fails at x<0.1
  fxbmax=0.7;

  fQ2min=0.1;
  fQ2max=5.;

  Double_t nuu=-1;
  Double_t ww=-1;

  Bool_t goodkin=kFALSE;

  while(!goodkin){

    fxb=fxbmin+(fxbmax-fxbmin)*fRan->Rndm();
    fQ2=fQ2min+(fQ2max-fQ2min)*fRan->Rndm();

  // Check if it's coherent

    nuu=fQ2/2./fxb/fm;
    ww=fm*fm+2*fm*nuu-fQ2;

    if(nuu<(fEbeam-0.1)) goodkin=kTRUE;

  }


}

TVector3* TGenBase::GetVertex(void)
{
  // Returns the pointer of the vertex TVector3
  return fVertex;
}


//_____________________________________________________________________________
 Int_t TGenBase::ComputeElectron(void) 
{ 
  //Sets the scattered electron data member.
  //
  //Warning: This method generates an electron "to the left". In the elastic 
  //setting, electrons are detected in the calorimeter and then must be 
  //generated to the right. All vectors will be rotated by 180 deg. in 
  //TGenElas::ComputeElas(). Be aware of a possible error when printing 
  //electron values before and after TGenElas::ComputeElas().
  Double_t pe,thetae;//,s;
  Double_t eel=feini->E();

  Int_t success=-1;

  if (!fFermi) {
    //s=fQ2*(1.-fxb)/fxb+TMath::Power(fm,2.);
    //pe=(s-fm*fm)/(2.*fm*(fxb-1.))+eel;
    pe=-fQ2/(2.*fm*fxb)+eel;
    // If pe<0, then the electron has lost too much energy through RC
    // for the equations to be true. In reality this event will not be detected
    // so we just set the momentum of the scattered electron well outside the
    // spectrometer acceptance. The event will be however considered as
    // "generated", so that normalization can be computed correctly later.
    if(pe<0.) pe=1000.;
    thetae=TMath::ACos(-fQ2/(2.*pe*eel)+1.);
  }else{
    // This is the solution of a yucky quadratic equation for the general case
    // inelastic scattering off a nucleon
    Double_t pz=fpini->Pz();
    Double_t px=fpini->Px();
    Double_t Ep=fpini->E();
    Double_t pz2=TMath::Power(pz,2);//
    Double_t q4=TMath::Power(fQ2,2.);
    Double_t Ep2=TMath::Power(Ep,2.);//
    Double_t k2=TMath::Power(eel,2.);
    Double_t xb2=TMath::Power(fxb,2.);
    Double_t px2=TMath::Power(px,2.);//
    Double_t k3=TMath::Power(eel,3.);
    Double_t k4=TMath::Power(k2,2.);
    Double_t pxz2=(px2 + pz2);
    pe=(eel*pz*fQ2 + 
        2*Ep2*k2*fxb + 2*k2*pz2*fxb + px2*fQ2*fxb + pz2*fQ2*fxb - 
        Ep*(eel*fQ2 + 4*k2*pz*fxb + pz*fQ2*fxb) - 
        px*TMath::Sqrt(fQ2*fxb)*TMath::Sqrt(2*eel*pz*fQ2 + 4*k2*pz2*fxb + 
	Ep2*(4*k2 - fQ2)*fxb + pxz2*fQ2*fxb - 
	2*Ep*eel*(fQ2 + 4*eel*pz*fxb)))/(2.*eel*TMath::Power(Ep - pz,2.)*fxb);
 
    
    if(pe<0. || isnan(pe)==1) pe=1000.; // See commentary just above

    thetae=(4*k3*(-Ep + pz)*fQ2*fxb + eel*(Ep + pz)*q4*fxb + 
       4*k4*TMath::Power(Ep - pz,2.)*xb2 + Ep*pz*q4*xb2 + k2*fQ2*
       (fQ2 + 2.*(-Ep2 + pz2)*xb2) - px*TMath::Power(fQ2*fxb,1.5)*
       TMath::Sqrt(2*eel*pz*fQ2 + 4*k2*pz2*fxb + 
       Ep2*(4*k2 - fQ2)*fxb + pxz2*fQ2*fxb - 
       2*Ep*eel*(fQ2 + 4*eel*pz*fxb)))/(4*k3*(-Ep + pz)*fQ2*fxb + 
       2*eel*pz*q4*fxb + 4*k4*TMath::Power(Ep - pz,2.)*xb2 + 
       pxz2*q4*xb2 + k2*fQ2*(fQ2 + 4*pz*(-Ep + pz)*xb2));
    thetae=TMath::ACos(thetae);
    if(isnan(thetae)) thetae=3.14;// See commentary just above
  }

  // cout<<"fe ptr = "<<fe<<endl;
  // cout<<"fq ptr = "<<fq<<endl;
    
  if(!fe) {
    fe=new TLorentzVector(pe*TMath::Sin(thetae),0.,pe*TMath::Cos(thetae),pe);
  }else{
    fe->SetPxPyPzE(pe*TMath::Sin(thetae),0.,pe*TMath::Cos(thetae),pe);
  }
  if(!fq) {
  fq=new TLorentzVector(feini->Px()-fe->Px(),feini->Py()-fe->Py(),feini->Pz()-fe->Pz(),feini->E()-fe->E());
  }else{
    fq->SetPxPyPzE(feini->Px()-fe->Px(),feini->Py()-fe->Py(),feini->Pz()-fe->Pz(),feini->E()-fe->E());
  }

  if(feini->E()<fe->E()) {success=0;} 
  else {success=1;}

  return success;

}

//_____________________________________________________________________________
TLorentzVector* TGenBase::GenFermiIni(void)
{   
  // Generates an initial (later recoil) particle according to a fermi distribution
  // previously read from a file via the Init() method.
  

  if(!fgIsInit) Init();

  if(fdmom==0 || frho==0){
    cout<<"You must initialize the event first !"<<endl;
    exit(1);
  }
  
  Double_t costhetfermi=2.*fRan->Rndm()-1.;
  Double_t thetfermi=TMath::ACos(costhetfermi);
  Double_t phifermi=2.*TMath::Pi()*fRan->Rndm();
  Int_t keepfermi = 0;
  Double_t slope,rho_p,pfermi;
  
  while (!keepfermi) {
    pfermi=1000.*fRan->Rndm(); // in MeV
    Int_t nn=0;
    while(pfermi > fdmom[nn]) {
      nn++;
    }
    
    if (nn<200 && nn>0) {
      slope=(frho[nn-1]-frho[nn])/(fdmom[nn-1]-fdmom[nn]);
      rho_p=slope*pfermi+frho[nn]-slope*fdmom[nn];
    } else if (nn==0) {
      rho_p=0.;
    } else {
      rho_p=0.;
    }
    
    Double_t randn2 = fRan->Rndm();
    if (randn2<rho_p) keepfermi=1;
  }
  pfermi*=0.001;  // In GeV
  
  Double_t px=pfermi*TMath::Sin(thetfermi)*TMath::Cos(phifermi);
  Double_t py=pfermi*TMath::Sin(thetfermi)*TMath::Sin(phifermi);
  Double_t pz=pfermi*TMath::Cos(thetfermi);
  Double_t Ep=TMath::Sqrt(px*px+py*py+pz*pz+fm*fm);

  if(!fpini){
    fpini=new TLorentzVector(px,py,pz,Ep);
  }else{
    fpini->SetPxPyPzE(px,py,pz,Ep);
  }
  return fpini;

}

//_____________________________________________________________________________
 TLorentzVector* TGenBase::GetScatteredElectron(void)
{
  // Returns the scattered electron 4-vector if it exists
  
  if(!fe) cout<<"Warning : Scattered electron doesn't exist"<<endl;
  return fe;

}

 TLorentzVector* TGenBase::GetInitialElectron(void)
{
  // Returns the Initial electron 4-vector if it exists
  
  if(!feini) cout<<"Warning : initial electron doesn't exist"<<endl;
  return feini;
}

//_____________________________________________________________________________
 void TGenBase::ApplySpecVerAcc(Double_t aav)
{
  // Applies vertical spectrometer acceptance by rotating the scattered 
  // electron around the beam axis

  fe->RotateZ(aav);
  feprerad->RotateZ(aav);

}

//_____________________________________________________________________________
 void TGenBase::Init(void)
{
  // Initializes the fermi momentum spread. Reads from file "deut_momdist.dat"
  cout<<"this"<<endl;
  if (fFermi) {
    fdmom=new Double_t[200];
    frho=new Double_t[200];
    ifstream deutmomdist("deut_momdist.dat");
    cout << "Using Fermi momentum spread (D2 target)" << endl;
    if(!deutmomdist) {
      cout << "Cannot open file deut_momdist.dat";
      exit(1);
    }
    for(Int_t ind=0;ind<200;ind++) {
      deutmomdist >> fdmom[ind];
      deutmomdist >> frho[ind];
    }
  }
  fgIsInit=kTRUE;
}

//_____________________________________________________________________________
 void TGenBase::Clear(char* opt)
{
  feini->SetPxPyPzE(0.,0.,fEbeam,fEbeam);
  if (fpini) fpini->SetPxPyPzE(0.,0.,0.,fm);
}

//_____________________________________________________________________________
void TGenBase::CloseTmpFile()
{
  output->close();
}

//_____________________________________________________________________________
void TGenBase::DumpFinalFile(char* finalfile, Int_t Ngen, Int_t Nacc)
{

  cout<<Form("Dumping final file from %s...", fOutTempFile.Data())<<endl;
  ifstream input(fOutTempFile);
  ofstream totofile(finalfile);
  totofile<<Ngen<<endl;
  totofile<<Nacc<<endl;
  totofile<<fEbeam<<" "<<fSpecAngle<<" "<<fSpecMom<<endl;
  totofile<<fCaloDist<<" "<<fCaloAngle<<" "<<endl;
  totofile<<fTargType<<" "<<fTargLength<<" "<<fTargDens<<endl;
  totofile<<fFermi<<" "<<fRadCor<<endl;

  Double_t dum[fNwrite];

  //  cout << fNwrite << endl;

  if(Nacc>0){ 
    while(!input.eof()){
      for(Int_t i=0;i<fNwrite;i++){
	input>>dum[i];
      }    
      if(!input.eof()) {
	for(Int_t i=0;i<fNwrite;i++){
	  totofile<<dum[i]<<" ";
	}
	totofile<<endl;
      }
    }
  }
  input.close();
 // remove(fOutTempFile);
  totofile.close();
}

//_____________________________________________________________________________
/*void TGenBase::DumpRootFile(char* finalfile, Int_t Ngen, Int_t Nacc){{{*/
void TGenBase::DumpRootFile(char* finalfile, Int_t Ngen, Int_t Nacc)
{

	if(Nacc>0){ 
		cout<<Form("--> Dumping #%d events from %s into the ROOT file ...", Nacc, fOutTempFile.Data())<<endl;
		ifstream input(fOutTempFile);

		Double_t vertexz, E0; 
		Double_t ePx_ini, ePy_ini, ePz_ini, hPx_ini, hPy_ini,hPz_ini;
		Double_t ePx, ePy, ePz, gPx, gPy, gPz, hPx, hPy, hPz;
		Double_t eP_ini,hP_ini, eP, gP,hP;
		Double_t Q2, x, t, phi, XS_P, XS_M, XS_BHp, XS_BHm, PSF;
		Double_t ePhi_ini,hPhi_ini, ePhi, gPhi,hPhi;
		Double_t eTheta_ini,hTheta_ini, eTheta, gTheta,hTheta;

		TFile *rootfile = new TFile(finalfile,"recreate");
		TTree *T = new TTree("T","A new tree");

		T->Branch("vertexz", &vertexz, "vertexz/D");
		T->Branch("E0", &E0, "E0/D");

		T->Branch("eP_ini", &eP_ini, "eP_ini/D");
		T->Branch("ePx_ini", &ePx_ini, "ePx_ini/D");
		T->Branch("ePy_ini", &ePy_ini, "ePy_ini/D");
		T->Branch("ePz_ini", &ePz_ini, "ePz_ini/D");
		T->Branch("eTheta_ini",&eTheta_ini, "eTheta_ini/D");
		T->Branch("ePhi_ini",&ePhi_ini, "ePhi_ini/D");
		
		if(fFermi){	
			T->Branch("hP_ini", &hP_ini, "hP_ini/D");
			T->Branch("hPx_ini", &hPx_ini, "hPx_ini/D");
			T->Branch("hPy_ini", &hPy_ini, "hPy_ini/D");
			T->Branch("hPz_ini", &hPz_ini, "hPz_ini/D");
			T->Branch("hTheta_ini",&hTheta_ini, "hTheta_ini/D");
			T->Branch("hPhi_ini",&hPhi_ini, "hPhi_ini/D");
		}

		T->Branch("eP", &eP, "eP/D");
		T->Branch("ePx", &ePx, "ePx/D");
		T->Branch("ePy", &ePy, "ePy/D");
		T->Branch("ePz", &ePz, "ePz/D");
		T->Branch("eTheta",&eTheta, "eTheta/D");
		T->Branch("ePhi",&ePhi, "ePhi/D");
	
		T->Branch("gP", &gP, "gP/D");
		T->Branch("gPx", &gPx, "gPx/D");
		T->Branch("gPy", &gPy, "gPy/D");
		T->Branch("gPz", &gPz, "gPz/D");
		T->Branch("gTheta", &gTheta, "gTheta/D");
		T->Branch("gPhi", &gPhi, "gPhi/D");

		T->Branch("hP", &hP, "hP/D");
		T->Branch("hPx", &hPx, "hPx/D");
		T->Branch("hPy", &hPy, "hPy/D");
		T->Branch("hPz", &hPz, "hPz/D");
		T->Branch("hTheta",&hTheta, "hTheta/D");
		T->Branch("hPhi",&hPhi, "hPhi/D");
	
		T->Branch("Q2", &Q2, "Q2/D");
		T->Branch("x", &x, "x/D");
		T->Branch("t", &t, "t/D");
		T->Branch("phi", &phi, "phi/D");
		T->Branch("XS_P", &XS_P, "XS_P/D");
		T->Branch("XS_M", &XS_M, "XS_M/D");
		T->Branch("PSF", &PSF, "PSF/D");
		T->Branch("Ngen", &Ngen, "Ngen/I");
		T->Branch("Nacc", &Nacc, "Nacc/I");
		T->Branch("XS_BHp", &XS_BHp, "XS_BHp/D");
		T->Branch("XS_BHm", &XS_BHm, "XS_BHm/D");


		for(int i=0;i<Nacc;i++){
			input >> vertexz >> E0>> ePx_ini >> ePy_ini >> ePz_ini;
			if(fFermi)
				input>> hPx_ini >> hPy_ini >> hPz_ini;
			input >> ePx >> ePy >> ePz >> gPx >> gPy >> gPz >> hPx >> hPy >> hPz
				  >> Q2 >> x >> t >> phi >> XS_P >> XS_M >> XS_BHp >> XS_BHm >> PSF;

			phi *= 180.0/3.1415926;
			
			eP_ini = sqrt( pow(ePx_ini,2)+pow(ePy_ini,2)+pow(ePz_ini,2));
			hP_ini = sqrt( pow(hPx_ini,2)+pow(hPy_ini,2)+pow(hPz_ini,2));
			eP = sqrt( pow(ePx,2)+pow(ePy,2)+pow(ePz,2));
			gP = sqrt( pow(gPx,2)+pow(gPy,2)+pow(gPz,2));
			hP = sqrt( pow(hPx,2)+pow(hPy,2)+pow(hPz,2));
		
			eTheta_ini = atan(sqrt(ePx_ini*ePx_ini+ePy_ini*ePy_ini)/ePz_ini)/TMath::DegToRad();
			if(fFermi)
				hTheta_ini = atan(sqrt(hPx_ini*hPx_ini+hPy_ini*hPy_ini)/hPz_ini)/TMath::DegToRad();
			eTheta = atan(sqrt(ePx*ePx+ePy*ePy)/ePz)/TMath::DegToRad();
			gTheta = atan(sqrt(gPx*gPx+gPy*gPy)/gPz)/TMath::DegToRad();
			hTheta = atan(sqrt(hPx*hPx+hPy*hPy)/hPz)/TMath::DegToRad();
			
			ePhi_ini = atan2(ePy_ini,ePx_ini)/TMath::DegToRad();
			if(fFermi)
				hPhi_ini = atan2(hPy_ini,hPx_ini)/TMath::DegToRad();
			ePhi = atan2(ePy,ePx)/TMath::DegToRad();
			gPhi = atan2(gPy,gPx)/TMath::DegToRad();
			hPhi = atan2(hPy,hPx)/TMath::DegToRad();

			T->Fill();
		}
		input.close();

		rootfile->cd(); T->Write(); rootfile->Close();

		remove(fOutTempFile);
	}else{
		cout<<Form("*** #%d events: failed to create the ROOT file ...", Nacc)<<endl;
	}
}
/*void TGenBase::DumpRootFile(char* finalfile, Int_t Ngen, Int_t Nacc)}}}*/

//_____________________________________________________________________________
/*void TGenBase::DumpRootFilePi0(char* finalfile, Int_t Ngen, Int_t Nacc){{{*/
void TGenBase::DumpRootFilePi0(char* finalfile, Int_t Ngen, Int_t Nacc)
{//Added in 03/20/2015 by Z. Ye

	if(Nacc>0){ 
		cout<<Form("--> Dumping #%d events into the ROOT file ...", Nacc)<<endl;
		ifstream input(fOutTempFile);

		Double_t vertexz, E0; 
		Double_t ePx_ini, ePy_ini, ePz_ini, hPx_ini, hPy_ini,hPz_ini;
		Double_t ePx, ePy, ePz, hPx, hPy, hPz;
		Double_t g1Px, g1Py, g1Pz, g2Px, g2Py, g2Pz;
		Double_t eP_ini,hP_ini, eP, g1P, g2P,hP;
		Double_t Q2, x, t, phi, XS_P, XS_M,PSF;
		Double_t ePhi_ini,hPhi_ini, ePhi, g1Phi, g2Phi,hPhi;
		Double_t eTheta_ini,hTheta_ini, eTheta, g1Theta, g2Theta,hTheta;

		TFile *rootfile = new TFile(finalfile,"recreate");
		TTree *T = new TTree("T","A new tree");

		T->Branch("vertexz", &vertexz, "vertexz/D");
		T->Branch("E0", &E0, "E0/D");

		T->Branch("eP_ini", &eP_ini, "eP_ini/D");
		T->Branch("ePx_ini", &ePx_ini, "ePx_ini/D");
		T->Branch("ePy_ini", &ePy_ini, "ePy_ini/D");
		T->Branch("ePz_ini", &ePz_ini, "ePz_ini/D");
		T->Branch("eTheta_ini",&eTheta_ini, "eTheta_ini/D");
		T->Branch("ePhi_ini",&ePhi_ini, "ePhi_ini/D");
		
		if(fFermi){	
			T->Branch("hP_ini", &hP_ini, "hP_ini/D");
			T->Branch("hPx_ini", &hPx_ini, "hPx_ini/D");
			T->Branch("hPy_ini", &hPy_ini, "hPy_ini/D");
			T->Branch("hPz_ini", &hPz_ini, "hPz_ini/D");
			T->Branch("hTheta_ini",&hTheta_ini, "hTheta_ini/D");
			T->Branch("hPhi_ini",&hPhi_ini, "hPhi_ini/D");
		}

		T->Branch("eP", &eP, "eP/D");
		T->Branch("ePx", &ePx, "ePx/D");
		T->Branch("ePy", &ePy, "ePy/D");
		T->Branch("ePz", &ePz, "ePz/D");
		T->Branch("eTheta",&eTheta, "eTheta/D");
		T->Branch("ePhi",&ePhi, "ePhi/D");
	
		T->Branch("g1P",    &g1P,    "g1P/D");
		T->Branch("g1Px",   &g1Px,   "g1Px/D");
		T->Branch("g1Py",   &g1Py,   "g1Py/D");
		T->Branch("g1Pz",   &g1Pz,   "g1Pz/D");
		T->Branch("g1Theta",&g1Theta,"g1Theta/D");
		T->Branch("g1Phi",  &g1Phi,  "g1Phi/D");
	
		T->Branch("g2P",    &g2P,    "g2P/D");
		T->Branch("g2Px",   &g2Px,   "g2Px/D");
		T->Branch("g2Py",   &g2Py,   "g2Py/D");
		T->Branch("g2Pz",   &g2Pz,   "g2Pz/D");
		T->Branch("g2Theta",&g2Theta,"g2Theta/D");
		T->Branch("g2Phi",  &g2Phi,  "g2Phi/D");

		T->Branch("hP", &hP, "hP/D");
		T->Branch("hPx", &hPx, "hPx/D");
		T->Branch("hPy", &hPy, "hPy/D");
		T->Branch("hPz", &hPz, "hPz/D");
		T->Branch("hTheta",&hTheta, "hTheta/D");
		T->Branch("hPhi",&hPhi, "hPhi/D");
	
		T->Branch("Q2", &Q2, "Q2/D");
		T->Branch("x", &x, "x/D");
		T->Branch("t", &t, "t/D");
		T->Branch("phi", &phi, "phi/D");
		T->Branch("XS_P", &XS_P, "XS_P/D");
		T->Branch("XS_M", &XS_M, "XS_M/D");
		T->Branch("PSF", &PSF, "PSF/D");
		T->Branch("Ngen", &Ngen, "Ngen/I");
		T->Branch("Nacc", &Nacc, "Nacc/I");

		for(int i=0;i<Nacc;i++){
			input >> vertexz >> E0>> ePx_ini >> ePy_ini >> ePz_ini;
			if(fFermi)
				input>> hPx_ini >> hPy_ini >> hPz_ini;
			input >> ePx >> ePy >> ePz >> g1Px >> g1Py >> g1Pz >> g2Px >> g2Py >> g2Pz >> hPx >> hPy >> hPz
				  >> Q2 >> x >> t >> phi >> XS_P >> XS_M >> PSF;

			phi *= 180.0/3.1415926;

			eP_ini = sqrt( pow(ePx_ini,2)+pow(ePy_ini,2)+pow(ePz_ini,2));
			hP_ini = sqrt( pow(hPx_ini,2)+pow(hPy_ini,2)+pow(hPz_ini,2));
			eP = sqrt( pow(ePx,2)+pow(ePy,2)+pow(ePz,2));
			g1P = sqrt( pow(g1Px,2)+pow(g1Py,2)+pow(g1Pz,2));
			g2P = sqrt( pow(g2Px,2)+pow(g2Py,2)+pow(g2Pz,2));
			hP = sqrt( pow(hPx,2)+pow(hPy,2)+pow(hPz,2));
		
			eTheta_ini = atan(sqrt(ePx_ini*ePx_ini+ePy_ini*ePy_ini)/ePz_ini)/TMath::DegToRad();
			if(fFermi)
				hTheta_ini = atan(sqrt(hPx_ini*hPx_ini+hPy_ini*hPy_ini)/hPz_ini)/TMath::DegToRad();
			eTheta = atan(sqrt(ePx*ePx+ePy*ePy)/ePz)/TMath::DegToRad();
			g1Theta = atan(sqrt(g1Px*g1Px+g1Py*g1Py)/g1Pz)/TMath::DegToRad();
			g2Theta = atan(sqrt(g2Px*g2Px+g2Py*g2Py)/g2Pz)/TMath::DegToRad();
			hTheta = atan(sqrt(hPx*hPx+hPy*hPy)/hPz)/TMath::DegToRad();
			
			ePhi_ini = atan2(ePy_ini,ePx_ini)/TMath::DegToRad();
			if(fFermi)
				hPhi_ini = atan2(hPy_ini,hPx_ini)/TMath::DegToRad();
			ePhi = atan2(ePy,ePx)/TMath::DegToRad();
			g1Phi = atan2(g1Py,g1Px)/TMath::DegToRad();
			g2Phi = atan2(g2Py,g2Px)/TMath::DegToRad();
			hPhi = atan2(hPy,hPx)/TMath::DegToRad();

			T->Fill();
		}
		input.close();

		rootfile->cd(); T->Write(); rootfile->Close();

		remove(fOutTempFile);
	}else{
		cout<<Form("*** #%d events: failed to create the ROOT file ...", Nacc)<<endl;
	}
}
/*void TGenBase::DumpRootFilePi0(char* finalfile, Int_t Ngen, Int_t Nacc)}}}*/

//_____________________________________________________________________________
/*void TGenBase::Print(char* opt){{{*/
void TGenBase::Print(char* opt)
{
	// Output on screen

	cout<<"================================================================="<<endl;
	cout<<"TGenBase"<<endl;
	cout<<"================================================================="<<endl;
	cout<<"Energy beam : "<<fEbeam<<endl;
	cout<<"------"<<endl;
  if(fFermi){
    cout<<"Fermi momentum enable"<<endl;
  }else{
    cout<<"Fermi momentum disable"<<endl;
  }
  if(fRadCor){
    cout<<"Radiative corrections enable"<<endl;
  }else{
    cout<<"Radiative corrections disable"<<endl;
  }
  cout<<"------"<<endl;
  if(fgWarnings){
    cout<<"Warnings enable"<<endl;
  }else{
    cout<<"Warnings disable"<<endl;
  }
  if(fgErrors){
    cout<<"Errors enable"<<endl;
  }else{
    cout<<"Errors disable"<<endl;
  }
  cout<<"================================================================="<<endl;
}
/*void TGenBase::Print(char* opt)}}}*/
