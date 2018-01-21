//
// TGenDeltPi0.cxx, v1.0, Tue Aug 4 11:13:57
// Author: C. Munoz Camacho
//

#include <fstream>
#include <iostream>
#include <stdlib.h>

#ifndef __TGenDeltPi0__
#include "TGenDeltPi0.h"
#endif

using namespace std;

ClassImp(TGenDeltPi0)

////////////////////////////////////////////////////////////////////////////////
//
// Event generator base class
// 
////////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________
  TGenDeltPi0::TGenDeltPi0(Double_t Ebeam, Int_t TargType, UInt_t seed1, UInt_t seed2):TGenBase(Ebeam,seed1,seed2)
{
  // Default constructor
  // Initial 4-vectors are initialized

   // if(!fgIsInit) Init();

  cout<<"TGenDeltPi0 constructor"<<endl;

  //
  fdeltam=0;
  fTargType=TargType;

  if(fTargType!=1 && fTargType!=0){
    cout<<"Unknown target for DVCS Event generator"<<endl;
    exit(1);
  }

  if(fTargType==0) fm=PMass();
  if(fTargType==1) fm=NMass();
  if(fTargType==2) fm=DMass();

  fpini=0;

  ftmax=0.;
  ftmin=0.;

  fxb=0;
  fQ2=0;
  fe=0;
  //  fq=0;
  fpi0=0;
  fg1=0;
  fg2=0;
  fp=0;
    
}

//_____________________________________________________________________________
TGenDeltPi0::TGenDeltPi0(const TGenDeltPi0& gen):TGenBase(gen)
{
  // Copy constructor
  //  ((TGenDeltPi0&)TCalobase).Copy(*this);
}

//_____________________________________________________________________________
 TGenDeltPi0::~TGenDeltPi0()
{
  // Default destructor

  if(fpini) delete fpini;
  if(fe)    delete fe;
  if(fp)    delete fp;
  if(fpi0)   delete fpi0;
  if(fg1)    delete fg1;
  if(fg2)    delete fg2;
}

//_____________________________________________________________________________
void TGenDeltPi0::IntRCBef(void)
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
void TGenDeltPi0::IntRCAft(void)
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
 void TGenDeltPi0::ComputePi0(void) 
{ 
  //Computes the gamma* p -> pi0 p' reaction in the center of mass.
  //The angle phi is generated and all vectors are boost back to the 
  //laboratory. Data members fp and fg are set.

  Compute2Body(0.1349766);
}

//_____________________________________________________________________________
 void TGenDeltPi0::ComputeTminTmax(Double_t m3,Double_t m4) 
{ 

  // Computes the (almost) general case tmin and tmax for varying masses
  // of the final state particles. However, virtual photon and proton at
  // rest always used in this method (we use x,Q2 for virtual photon).

  Double_t ecm1,ecm3;
  Double_t m2;
  Double_t m1sq,m2sq,m3sq,m4sq;

  m2=fm;

  m1sq=-fQ2;
  m2sq=m2*m2;
  m3sq=m3*m3;
  m4sq=m4*m4;

  // Update tmin tmax. ATTENTION tmax is still... the usual tmin (=ft0)
  // tmin is the usual tmax (=ft1).. thank Carlos for that convention

  // fs is already calculated when you all this method !!!

  // CoM energies

  ecm1=(fs+m1sq-m2sq)/2./TMath::Sqrt(fs);
  ecm3=(fs+m3sq-m4sq)/2./TMath::Sqrt(fs);

  // CoM momentum

  fpcm1=TMath::Sqrt(ecm1*ecm1-m1sq);
  fpcm3=TMath::Sqrt(ecm3*ecm3-m3sq);

  // tmin/max calculations

  ft0=TMath::Power((m1sq-m3sq-m2sq+m4sq)/(2.*TMath::Sqrt(fs)),2.)-TMath::Power(fpcm1-fpcm3,2.);
  ft1=TMath::Power((m1sq-m3sq-m2sq+m4sq)/(2.*TMath::Sqrt(fs)),2.)-TMath::Power(fpcm1+fpcm3,2.);

  ftmax=TMath::Min(ftmax,ft0);
  ftmax-=0.0001;
  ftmin=TMath::Max(ftmin,ft1);
  ftmin+=0.0001;

//   cout <<"--------------------------------------------"<<endl;
//   cout << "ecm3="<<(fs+m3sq-m4sq)/2./TMath::Sqrt(fs)<<"    "<<"m3="<<TMath::Sqrt(m3sq)<<endl;
//   cout << "ecm4="<<(fs+m4sq-m3sq)/2./TMath::Sqrt(fs)<<"    "<<"m4="<<TMath::Sqrt(m4sq)<<endl;
//   cout << "mass="<<fdeltam<<endl;
//   cout << "t0  ="<<ft0<<"    "<<"t1  ="<<ft1<<endl;
//   cout << "tmax="<<ftmax<<"    "<<"tmin="<<ftmin<<endl;

}

//_____________________________________________________________________________
 void TGenDeltPi0::Compute2Body(Double_t m) 
{ 
  //Computes the gamma* p -> X p' reaction in the center of mass, as a
  //function of the mass m of the particle X.
  //The angle phi is generated and all vectors are boost back to the 
  //laboratory. Data members fp and fg are set.
  
  //////////////////////
  //Computation of t and the phase space factor
  //////////////////////

  Double_t nu=fQ2/(2.*fm*fxb);
  // BEWARE:Calculated for non-fermi smeared kinematics
  Double_t q3=TMath::Sqrt(fQ2+TMath::Power(nu,2.));


  // Tmin Tmax depend on recoil masses... have to calculate it right !

  if(fpini==0 && fFermi==kFALSE) fpini=new TLorentzVector(0.,0.,0.,fm);
  if(fpini==0 && fFermi==kTRUE){
    cout<<"You must generate fermi recoil particle first"<<endl;
    exit(1);
  }

  TLorentzVector p=*fpini;
  TLorentzVector q=*fq;

  // Boost to center of mass p-q

  TLorentzVector cms = q + p;

  fs=cms.M2();

  q.Boost(-cms.BoostVector() );
  p.Boost(-cms.BoostVector() );

  if (fs>1.5) {

    fdeltam=0;

    while (fdeltam<1.0732 || fdeltam>(TMath::Sqrt(fs)-0.1349766)) {
      fdeltam=fRan->BreitWigner(1.232,0.110);
    }

  } else fdeltam=1.232; // this event will not make it further
                        // than spectro cut anyway

  // Now for the "real" values of tmin and tmax
  ComputeTminTmax(0.1349766,fdeltam);
  //
  
  if(fFermi) ftmax=0;

  ft=ftmax+(ftmin-ftmax)*fRan->Rndm();

  fPSF=(ftmax-ftmin)*(fxbmax-fxbmin)*(fQ2max-fQ2min);
  ///////////////////

  //  fq=new TLorentzVector(feini->Px()-fe->Px(),feini->Py()-fe->Py(),feini->Pz()-fe->Pz(),feini->E()-fe->E());  

  // Rotation around Y
  
  TVector3 oz(0.,0.,1.);
  TVector3 perpvec=(q.Vect()).Cross(oz);
  Double_t angle=(q.Vect()).Angle(oz);

  q.Rotate(angle,perpvec.Unit());
  p.Rotate(angle,perpvec.Unit());

  //Double_t egammacm=(fs-fm*fm)/(2.*TMath::Sqrt(fs));
  //Double_t thetacm = (ft+fQ2)/(2.*egammacm*q.P())+q(3)/q.P();

  //

  //////////////////////////////////////////

  // Following concerns gamma in the reaction
  //  e p -> e pi0 Delta+, Delta+ -> N pi (whatever N or pi)
  // only e and pi0 (2 gamma) info go into the final file

  //////////////////////////////////////////

  Double_t ecm=(fs-fdeltam*fdeltam+m*m)/(2.*TMath::Sqrt(fs));
  Double_t pcm=TMath::Sqrt(TMath::Power(ecm,2.)-TMath::Power(m,2.));
  Double_t thetacm=2.*TMath::ASin(TMath::Sqrt((ft0-ft)/4./fpcm1/fpcm3));

  if (!fpi0) {fpi0=new TLorentzVector(pcm*TMath::Sin(thetacm),0.,TMath::Cos(thetacm)*pcm,ecm);}
  else {fpi0->SetPxPyPzE(pcm*TMath::Sin(thetacm),0.,TMath::Cos(thetacm)*pcm,ecm);}
  if (!fp) {fp=new TLorentzVector(q.Px()+p.Px()-fpi0->Px(),q.Py()+p.Py()-fpi0->Py(),q.Pz()+p.Pz()-fpi0->Pz(),q.E()+p.E()-fpi0->E());} else {fp->SetPxPyPzE(q.Px()+p.Px()-fpi0->Px(),q.Py()+p.Py()-fpi0->Py(),q.Pz()+p.Pz()-fpi0->Pz(),q.E()+p.E()-fpi0->E());}
   
  // Rotation back
  fp->Rotate(-angle,perpvec.Unit());
  fpi0->Rotate(-angle,perpvec.Unit());
  q.Rotate(-angle,perpvec.Unit());
  
  // Boost to lab
  fp->Boost(cms.BoostVector());
  fpi0->Boost(cms.BoostVector());
  q.Boost(cms.BoostVector());

  // Now we have k along Oz and k' on the xOz plane
  // q=k-k' is necessarily on plane as well
  // We just have to rotate along Oy to bring q along Oz
  // and then rotate everything by phi except
  // k and k'. Then we rotate back along Oy.

  angle=(q.Vect()).Angle(oz);

  fp->RotateY(angle);
  fpi0->RotateY(angle);
  
  fphi=2.*TMath::Pi()*fRan->Rndm();   // phi between 0 and 2pi
  
  fp->RotateZ(fphi);
  fpi0->RotateZ(fphi);

  fp->RotateY(-angle);
  fpi0->RotateY(-angle);

}

//_____________________________________________________________________________
void TGenDeltPi0::TwoBodyDecay(Double_t M, Double_t m1, Double_t m2)
{
  //Desintegration du pi0

  Double_t cosalpha=2.*fRan->Rndm()-1.; // Cos alpha uniformly between -1 and 1
  Double_t phi2=2.*TMath::Pi()*fRan->Rndm(); // phi between 0 and 2*pi

  if (!fg1) {fg1=new TLorentzVector(0.5*M*TMath::Sin(TMath::ACos(cosalpha))*TMath::Cos(phi2),0.5*M*TMath::Sin(TMath::ACos(cosalpha))*TMath::Sin(phi2),0.5*M*cosalpha,0.5*M);} else 
    {fg1->SetPxPyPzE(0.5*M*TMath::Sin(TMath::ACos(cosalpha))*TMath::Cos(phi2),0.5*M*TMath::Sin(TMath::ACos(cosalpha))*TMath::Sin(phi2),0.5*M*cosalpha,0.5*M);}
  if (!fg2) {fg2=new TLorentzVector(-0.5*M*TMath::Sin(TMath::ACos(cosalpha))*TMath::Cos(phi2),-0.5*M*TMath::Sin(TMath::ACos(cosalpha))*TMath::Sin(phi2),-0.5*M*cosalpha,0.5*M);} else 
    {fg2->SetPxPyPzE(-0.5*M*TMath::Sin(TMath::ACos(cosalpha))*TMath::Cos(phi2),-0.5*M*TMath::Sin(TMath::ACos(cosalpha))*TMath::Sin(phi2),-0.5*M*cosalpha,0.5*M);}

  fg1->Boost(fpi0->BoostVector()); // boost au labo 
  fg2->Boost(fpi0->BoostVector()); // boost au labo

}

//_____________________________________________________________________________
void TGenDeltPi0::ApplySpecVerAcc(Double_t aav)
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
  fpi0->RotateZ(av);
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
 TLorentzVector* TGenDeltPi0::GetFinalPhoton1(void)
{
  // Returns the final photon1 4-vector if it exists
  
  if(!fg1) cout<<"Warning : Final photon 1 doesn't exist"<<endl;
  return fg1;

}

//_____________________________________________________________________________
 TLorentzVector* TGenDeltPi0::GetFinalPhoton2(void)
{
  // Returns the final photon2 4-vector if it exists
  
  if(!fg2) cout<<"Warning : Final photon 2 doesn't exist"<<endl;
  return fg2;

}

//_____________________________________________________________________________
 TLorentzVector* TGenDeltPi0::GetFinalProton(void)
{
  // Returns the final photon 4-vector if it exists
  
  if(!fp) cout<<"Warning : Final proton doesn't exist"<<endl;
  return fp;

}

//_____________________________________________________________________________
void TGenDeltPi0::Write2File(void)
{

  if(!fNwrite) {
    if(!fFermi) {
      fNwrite=18;
    }else{
      fNwrite=18+3;
    }
  }
  
  *output<<fVertex->Pz()<<" ";
  *output<<feini->Pz()<<" ";
  *output<<feprerad->Px()<<" "<<feprerad->Py()<<" "<<feprerad->Pz()<<" ";

  if(fFermi)
    *output<<fpini->Px()<<" "<<fpini->Py()<<" "<<fpini->Pz()<<" ";

  *output<<fe->Px()<<" "<<fe->Py()<<" "<<fe->Pz()<<" ";
  *output<<fg1->Px()<<" "<<fg1->Py()<<" "<<fg1->Pz()<<" ";
  *output<<fg2->Px()<<" "<<fg2->Py()<<" "<<fg2->Pz()<<" ";
  *output<<fp->Px()<<" "<<fp->Py()<<" "<<fp->Pz()<<" ";
  *output<< fPSF;
  *output<<endl;

//   TVector3 qpt(fg1->Px()+fg2->Px(),fg1->Py()+fg2->Py(),fg1->Pz()+fg2->Pz());
//   TVector3 qt(-fe->Px(),-fe->Py(),feini->Pz()-fe->Pz());
//   TVector3 kt(0,0,feini->Pz());
//   TVector3 qk=qt.Cross(kt);//perp at leptonic plan (initial electron,virtual photon)
//   TVector3 qqp=qt.Cross(qpt);//perp at hadronic plan (final pi0, virtual photon)

//   Double_t ang=(180./TMath::Pi())*(qk.Angle(qqp));//angle between the 2 plans

//   if(qt.Dot(qk.Cross(qqp))<0) ang=360.-ang;

//   cout << "Phi_gen="<<fphi*180./3.141592 << endl;

//   cout << "pi0: (" << qpt.X()<<","<<qpt.Y()<<","<<qpt.Z()<<")"<<endl;

//   cout << "pi0 unit: (" << qpt.Unit().X()<<","<<qpt.Unit().Y()<<","<<qpt.Unit().Z()<<")"<<endl;
//   cout << "  k: (" << kt.X()<<","<<kt.Y()<<","<<kt.Z()<<")"<<endl;
//   cout << " kp: (" << fe->Px()<<","<<fe->Py()<<","<<fe->Pz()<<")"<<endl;
//   cout << "  q: (" << qt.X()<<","<<qt.Y()<<","<<qt.Z()<<")"<<endl;
//   cout << "  q unit: (" << qt.Unit().X()<<","<<qt.Unit().Y()<<","<<qt.Unit().Z()<<")"<<endl;
//   cout << " qk: (" << qk.X()<<","<<qk.Y()<<","<<qk.Z()<<")"<<endl;
//   cout << " qqp: (" << qqp.X()<<","<<qqp.Y()<<","<<qqp.Z()<<")"<<endl;
//   cout << "Phi_rec=" << ang << endl;
//   cout << "Q2="<<fQ2<<endl;
//   cout << fphi*180./3.141592 << "   " << ang << endl;
}

//_____________________________________________________________________________
 void TGenDeltPi0::Print(char* opt)
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
