//
// TGenOme.cxx, v1.0, Tue Aug 4 11:13:57
// Author: C. Munoz Camacho
//

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <TGenPhaseSpace.h>

#ifndef __TGenOme__
#include "TGenOme.h"
#endif

using namespace std;

ClassImp(TGenOme)

////////////////////////////////////////////////////////////////////////////////
//
// Event generator base class
// 
////////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________
  TGenOme::TGenOme(Double_t Ebeam, Int_t TargType, UInt_t seed1, UInt_t seed2):TGenBase(Ebeam,seed1,seed2)
{
  // Default constructor
  // Initial 4-vectors are initialized

   // if(!fgIsInit) Init();

  cout<<"TGenOme constructor"<<endl;

  fdmom=0;
  fome=0;
  fTargType=TargType;

  if(fTargType!=1 && fTargType!=0){
    cout<<"Unknown target for DVCS Event generator"<<endl;
    exit(1);
  }

  if(fTargType==0) fm=PMass();
  if(fTargType==1) fm=NMass();

  fpini=0;

  ftmax=0.;
  ftmin=0.;

  fxb=0;
  fQ2=0;
  fe=0;
  //  fq=0;
  fome=0;
  fg1=0;
  fg2=0;
  fp=0;

}

//_____________________________________________________________________________
TGenOme::TGenOme(const TGenOme& gen):TGenBase(gen)
{
  // Copy constructor
  //  ((TGenOme&)TCalobase).Copy(*this);
}

//_____________________________________________________________________________
 TGenOme::~TGenOme()
{
  // Default destructor

  if(fpini) delete fpini;
  if(fe)    delete fe;
  if(fp)    delete fp;
  if(fome)   delete fome;
  if(fg1)    delete fg1;
  if(fg2)    delete fg2;
}

//_____________________________________________________________________________
void TGenOme::IntRCBef(void)
{
  // Makes internal radiative corrections _before_ the vertex (to the initial)
  // electron) using the equivalent radiator technique.
  // Uses \Delta E=E_0 * R^(2./nu) with R randomly between 0 and 1.
  // Factor 1/2 is because internal corrections must be applied twice (before
  // and after the vertex) with equivalent radiator thickness half each time.

  fRadCor=kTRUE;
  Double_t nu=2.*0.007297352533*(TMath::Log(fQ2/TMath::Power(0.000510998902,2))-1.)/TMath::Pi(); 
  
  Double_t eel=(feini->E())*(1-TMath::Power(fRan->Rndm(),1./(nu/2.)));
  
  feini->SetPxPyPzE(0.,0.,eel,eel);
}

//_____________________________________________________________________________
void TGenOme::IntRCAft(void)
{
  // Makes internal radiative corrections _after_ the vertex (to the scattered
  // electron) using the equivalent radiator technique.
  // Uses \Delta E=E_0 * R^(2./nu) with R randomly between 0 and 1.
  // Factor 1/2 is because internal corrections must be applied twice (before
  // and after the vertex) with equivalent radiator thickness half each time.

  fRadCor=kTRUE;
  Double_t nu=2.*0.007297352533*(TMath::Log(fQ2/TMath::Power(0.000510998902,2))-1.)/TMath::Pi(); 
  
  Double_t deel=1-TMath::Power(fRan->Rndm(),1./(nu/2.));

  *feprerad=(*fe); // keep vertex scattered electron
  
  *fe=deel*(*fe);
}

//_____________________________________________________________________________
 void TGenOme::ComputeOme(void) 
{ 
  //Computes the gamma* p -> ome p' reaction in the center of mass.
  //The angle phi is generated and all vectors are boost back to the 
  //laboratory. Data members fp and fg are set.

  // First off, get the ome mass

  // And then compute 2body kinematics with that.


//  FROM PIERRE:
//
//       data 
// e_ome,gam_ome,mpi0,mgam,mn/.770,.1505,.134976,.13957,.93957/   !valeur 
// du delta 1232

//       ei=pom(51)
     
//       ru=(p_i/2.)*(1.-2.*ranf(ur))
//       m_ome=gam_ome*.5*tan(ru)+e_ome
//       if(m_ome.gt.w)goto 41
//       if(m_ome.lt.(mgam+mpi0))goto 41

  Compute2Body();  // Mass will be randommized in this method
}

//_____________________________________________________________________________
 void TGenOme::Compute2Body() 
{ 
  //Computes the gamma* p -> X p' reaction in the center of mass, as a
  //function of the mass m of the particle X.
  //The angle phi is generated and all vectors are boost back to the 
  //laboratory. Data members fp and fg are set.
  
  //////////////////////
  //Computation of t and the phase space factor
  //////////////////////


  Double_t m;

  fomem=0.;

  Double_t nu=fQ2/(2.*fm*fxb);
  // BEWARE:Calculated for non-fermi smeared kinematics

  Double_t q3=TMath::Sqrt(fQ2+TMath::Power(nu,2.));

  Double_t q0primemax=0.5*fQ2*(1.-fxb)/(fxb*(fm+nu-q3));
  Double_t q0primemin=0.5*fQ2*(1.-fxb)/(fxb*(fm+nu+q3));
  Double_t tmax=-fQ2-2.*q0primemax*(nu-q3);
  Double_t tmin=-fQ2-2.*q0primemin*(nu+q3);
  ftmax=TMath::Min(ftmax,tmax);
  ftmin=TMath::Max(ftmin,tmin);
  
  if(fFermi) ftmax=0;

  ft=ftmax+(ftmin-ftmax)*fRan->Rndm();

  fPSF=-(ftmax-ftmin)*(fxbmax-fxbmin)*(fQ2max-fQ2min);
  ///////////////////

  //  fq=new TLorentzVector(feini->Px()-fe->Px(),feini->Py()-fe->Py(),feini->Pz()-fe->Pz(),feini->E()-fe->E());

  if(fpini==0 && fFermi==kFALSE) fpini=new TLorentzVector(0.,0.,0.,fm);
  if(fpini==0 && fFermi==kTRUE){
    cout<<"You must generate fermi recoil particle first"<<endl;
    exit(1);
  }
  
  TLorentzVector p=*fpini;
  TLorentzVector q=*fq;

  // Boost to center of mass p-q

  TLorentzVector cms = q + p;

//   cout << "e(" << feini->Px() << "," << feini->Py() << "," << feini->Pz()<< "," << feini->E()<< ")" << endl;
//   cout << "p(" << fpini->Px() << "," << fpini->Py() << "," << fpini->Pz()<< "," << fpini->E()<< ")" << endl;
//   cout << "e'(" << fe->Px() << "," << fe->Py() << "," << fe->Pz()<< "," << fe->E()<< ")" << endl;
//   cout << "q(" << q.Px() << "," << q.Py() << "," << q.Pz()<< "," << q.E()<< ")" << endl;
//   cout << "cms(" << cms.Px() << "," << cms.Py() << "," << cms.Pz() << "," << cms.E()<< ")" << endl;

  fs=cms.M2();

  q.Boost(-cms.BoostVector() );
  p.Boost(-cms.BoostVector() );

  // Rotation around Y
  
  TVector3 oz(0.,0.,1.);
  TVector3 perpvec=(q.Vect()).Cross(oz);
  Double_t angle=(q.Vect()).Angle(oz);

  q.Rotate(angle,perpvec.Unit());
  p.Rotate(angle,perpvec.Unit());

  //Double_t egammacm=(fs-fm*fm)/(2.*TMath::Sqrt(fs));
  //Double_t thetacm = (ft+fQ2)/(2.*egammacm*q.P())+q(3)/q.P();

  ///////////////////////////////////////// 

  fomem=0.78259;
  m=fomem;

  //////////////////////////////////////////

  Double_t ecm=(fs-fm*fm+m*m)/(2.*TMath::Sqrt(fs));
  Double_t pcm=TMath::Sqrt(TMath::Power(ecm,2.)-TMath::Power(m,2.));
  Double_t thetacm=(ft+fQ2-m*m+2.*q(3)*ecm)/(2.*pcm*q.P());

  if (thetacm>1.) thetacm=1.;
  if (thetacm<-1.) thetacm=-1.;
  thetacm=TMath::ACos(thetacm);

  if (!fome) {
    fome=new TLorentzVector(pcm*TMath::Sin(thetacm),0.,TMath::Cos(thetacm)*pcm,ecm);
  } else {
    fome->SetPxPyPzE(pcm*TMath::Sin(thetacm),0.,TMath::Cos(thetacm)*pcm,ecm);
  }

  if (!fp) {
    fp=new TLorentzVector(q.Px()+p.Px()-fome->Px(),q.Py()+p.Py()-fome->Py(),q.Pz()+p.Pz()-fome->Pz(),q.E()+p.E()-fome->E());
  } else {
    fp->SetPxPyPzE(q.Px()+p.Px()-fome->Px(),q.Py()+p.Py()-fome->Py(),q.Pz()+p.Pz()-fome->Pz(),q.E()+p.E()-fome->E());
  }

  // Rotation back
  fp->Rotate(-angle,perpvec.Unit());
  fome->Rotate(-angle,perpvec.Unit());
  q.Rotate(-angle,perpvec.Unit());
  
  // Boost to lab
  fp->Boost(cms.BoostVector());
  fome->Boost(cms.BoostVector());
  q.Boost(cms.BoostVector());

  // Now we have k along Oz and k' on the xOz plane
  // q=k-k' is necessarily on plane as well
  // We just have to rotate along Oy to bring q along Oz
  // and then rotate everything by phi except
  // k and k'. Then we rotate back along Oy.

  angle=(q.Vect()).Angle(oz);

  fp->RotateY(angle);
  fome->RotateY(angle);
  
  fphi=2.*TMath::Pi()*fRan->Rndm();   // phi between 0 and 2pi
  
  fp->RotateZ(fphi);
  fome->RotateZ(fphi);

  fp->RotateY(-angle);
  fome->RotateY(-angle);

}

//_____________________________________________________________________________
void TGenOme::TwoBodyDecay(Double_t M, Double_t m1, Double_t m2)
{
  // ome decay in gamma pi0

  Double_t masses[2];
  masses[0]=0.1349766; // pi0
  masses[1]=0;  // gam

  Double_t masses2[2];
  masses2[0]=0.; // gamma
  masses2[1]=0.;  // gamma

  TLorentzVector toto1;

  toto1.SetPx(fome->Px());
  toto1.SetPy(fome->Py());
  toto1.SetPz(fome->Pz());
  toto1.SetE(fome->E());

  TGenPhaseSpace ps;
      if ( !(ps.SetDecay(toto1,2,masses,"")) ) {
        cout << "Error in phase space generation" << endl;
      }

  Double_t wei = ps.Generate();

  if(!fpi0)fpi0=new TLorentzVector();
  if(!fgam)fgam=new TLorentzVector();

  TLorentzVector *q0decay = ps.GetDecay(0);
  TLorentzVector *qpdecay = ps.GetDecay(1);
  fpi0->SetPxPyPzE(q0decay->Px(),q0decay->Py(),q0decay->Pz(),q0decay->E());
  fgam->SetPxPyPzE(qpdecay->Px(),qpdecay->Py(),qpdecay->Pz(),qpdecay->E());

  fPSF*=wei;

  TLorentzVector toto2;

  toto2.SetPx(fpi0->Px());
  toto2.SetPy(fpi0->Py());
  toto2.SetPz(fpi0->Pz());
  toto2.SetE(fpi0->E());

  TGenPhaseSpace ps2;
      if ( !(ps2.SetDecay(toto2,2,masses2,"")) ) {
        cout << "Error in phase space generation" << endl;
      }

  Double_t wei2 = ps2.Generate();

  if(!fg1)fg1=new TLorentzVector();
  if(!fg2)fg2=new TLorentzVector();

  TLorentzVector *q1decay = ps2.GetDecay(0);
  TLorentzVector *q2decay = ps2.GetDecay(1);

  fg1->SetPxPyPzE(q1decay->Px(),q1decay->Py(),q1decay->Pz(),q1decay->E());
  fg2->SetPxPyPzE(q2decay->Px(),q2decay->Py(),q2decay->Pz(),q2decay->E());

  fPSF*=wei2;

}

//_____________________________________________________________________________
void TGenOme::ApplySpecVerAcc(Double_t aav)
{
  // Applies vertical spectrometer acceptance by rotating all 4-vectors
  // around the beam axis. An angle can be specified, otherwise it's generated
  // randomly between spectrometer acceptances

  Double_t av=aav;
  if(aav==-1.) {
    av=-fSpecVerAcc+2.*fSpecVerAcc*fRan->Rndm();
    av*=1/TMath::Sin(fSpecAngle);
  }
  TGenBase::ApplySpecVerAcc(av);
  fome->RotateZ(av);
  fg1->RotateZ(av);
  fg2->RotateZ(av);
  fp->RotateZ(av);
  //fe->RotateZ(av);
  fq->RotateZ(av);
  
  fPSF*=2.*fSpecVerAcc/TMath::Sin(fSpecAngle);

  // Remark:
  // The right phasespace factor is
  // dx.dq2.dt.dphi.dphie/2pi (because the x-section code gives 
  // phi_e-integrated x-sections) but dphi=2pi, so it simplifies to
  // dx.dq2.dt.dphie
}

//_____________________________________________________________________________
 TLorentzVector* TGenOme::GetFinalPhoton1(void)
{
  // Returns the final photon1 4-vector if it exists
  
  if(!fg1) cout<<"Warning : Final photon 1 doesn't exist"<<endl;
  return fg1;

}

//_____________________________________________________________________________
 TLorentzVector* TGenOme::GetFinalPhoton2(void)
{
  // Returns the final photon2 4-vector if it exists
  
  if(!fg2) cout<<"Warning : Final photon 2 doesn't exist"<<endl;
  return fg2;

}

//_____________________________________________________________________________
 TLorentzVector* TGenOme::GetFinalProton(void)
{
  // Returns the final photon 4-vector if it exists
  
  if(!fp) cout<<"Warning : Final proton doesn't exist"<<endl;
  return fp;

}

//_____________________________________________________________________________
void TGenOme::Write2File(void)
{

  if(!fNwrite) {
    if(!fFermi) {
      fNwrite=21;
    }else{
      fNwrite=21+3;
    }
  }
  
  *output<<fVertex->Pz()<<" ";
  *output<<feini->Pz()<<" ";
  *output<<feprerad->Px()<<" "<<feprerad->Py()<<" "<<feprerad->Pz()<<" ";
  if(fFermi)
    *output<<fpini->Px()<<" "<<fpini->Py()<<" "<<fpini->Pz()<<" ";
  *output<<fe->Px()<<" "<<fe->Py()<<" "<<fe->Pz()<<" ";
  *output<<fgam->Px()<<" "<<fgam->Py()<<" "<<fgam->Pz()<<" ";
  *output<<fg1->Px()<<" "<<fg1->Py()<<" "<<fg1->Pz()<<" ";
  *output<<fg2->Px()<<" "<<fg2->Py()<<" "<<fg2->Pz()<<" ";
  *output<<fp->Px()<<" "<<fp->Py()<<" "<<fp->Pz()<<" ";
  *output<< fPSF;
  *output<<endl;



}

//_____________________________________________________________________________
 void TGenOme::Print(char* opt)
{
  // Output on screen. If option "all" is specified the complete setup of the
  // event is printed out. By default only the final state 4-vectors and the 
  // virtual photon are displayed.

  TString option=opt;
  if(option.Contains("all")) TGenBase::Print();
  cout<<"======================================="<<endl;
  cout<<"        4-vectors (Px,Py,Pz,E)         "<<endl;
  cout<<"======================================="<<endl;
  if(feini) {
cout<<"e("<<feini->Px()<<","<<feini->Py()<<","<<feini->Pz()<<","<<feini->E()<<")"<<endl;
  }else{
    cout<<"NO INITIAL ELECTRON DEFINED"<<endl;
  }
  if(fpini){
    cout<<"p("<<fpini->Px()<<","<<fpini->Py()<<","<<fpini->Pz()<<","<<fpini->E()<<")"<<endl;
  }else{
    cout<<"NO INITIAL TARGET PARTICLE DEFINED"<<endl;
  }
  if(fq){
  cout<<"g*("<<fq->Px()<<","<<fq->Py()<<","<<fq->Pz()<<","<<fq->E()<<")"<<endl;
  }else{
    cout<<"NO VIRTUAL PHOTON DEFINED"<<endl;
  }
  if(fe){
    cout<<"e'("<<fe->Px()<<","<<fe->Py()<<","<<fe->Pz()<<","<<fe->E()<<")"<<endl;
  }else{
    cout<<"NO SCATTERED ELECTRON DEFINED"<<endl;
  }
  if(fg1){
    cout<<"g1("<<fg1->Px()<<","<<fg1->Py()<<","<<fg1->Pz()<<","<<fg1->E()<<")"<<endl;
  }else{
    cout<<"NO EMITTED PHOTON 1 DEFINED"<<endl;
  }
  if(fg2){
    cout<<"g2("<<fg2->Px()<<","<<fg2->Py()<<","<<fg2->Pz()<<","<<fg2->E()<<")"<<endl;
  }else{
    cout<<"NO EMITTED PHOTON 2 DEFINED"<<endl;
  }
  if(fp){
    cout<<"p'("<<fp->Px()<<","<<fp->Py()<<","<<fp->Pz()<<","<<fp->E()<<")"<<endl;
  }else{
    cout<<"NO RECOIL PARTICLE DEFINED"<<endl;
  }
  cout<<"======================================="<<endl;
}
