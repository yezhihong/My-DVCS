//
// TGenPip.cxx, v1.0
//

#include <fstream>
#include <iostream>
#include <stdlib.h>

#ifndef __TGenPip__
#include "TGenPip.h"
#endif

using namespace std;

ClassImp(TGenPip)
 
///////////////////////////////////////////////////////////////////////////////
//
// Event generator base class
// 
///////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________
  TGenPip::TGenPip(Double_t Ebeam, Int_t TargType, UInt_t seed1, UInt_t seed2):TGenBase(Ebeam,seed1,seed2)
{
  // Default constructor
  // Initial 4-vectors are initialized

   // if(!fgIsInit) Init();

  cout<<"TGenPip constructor"<<endl;

  fdmom=0;
  frho=0;
  fTargType=TargType;

  if(fTargType!=1 && fTargType!=0 && fTargType!=2){
    cout<<"Unknown target for DVCS Event generator"<<endl;
    exit(1);
  }

  fm=PMass();

  fpini=0;

  ftmax=0.;
  ftmin=0.;

  fxb=0;
  fQ2=0;
  fe=0;
  //  fq=0;
  fpip=0;
  fn=0;
    
}

//_____________________________________________________________________________
TGenPip::TGenPip(const TGenPip& gen):TGenBase(gen)
{
  // Copy constructor
  //  ((TGenPip&)TCalobase).Copy(*this);
}

//_____________________________________________________________________________
 TGenPip::~TGenPip()
{
  // Default destructor

  if(fpini) delete fpini;
  if(fe)    delete fe;
  if(fn)    delete fn;
  if(fpip)   delete fpip;
}

//_____________________________________________________________________________
void TGenPip::IntRCBef(void)
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
void TGenPip::IntRCAft(void)
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
 void TGenPip::ComputePip(void) 
{ 
  //Computes the gamma* p -> pip n reaction in the center of mass.
  //The angle phi is generated and all vectors are boost back to the 
  //laboratory. Data members fn and fpip are set.

  Compute2Body(0.13957018);
}

//_____________________________________________________________________________
 void TGenPip::ComputeTminTmax(Double_t m3,Double_t m4) 
{ 

  // Computes the (almost) general case tmin and tmax for varying masses
  // of the final state particles. However, virtual photon and proton at
  // rest always used in this method (we use x,Q2 for virtual photon).

  Double_t m2;
  Double_t m1sq,m2sq,m3sq,m4sq;

  m2=fm;     // part1=gamma*
             // part2=proton
  m1sq=-fQ2; // part3=pi+
  m2sq=m2*m2;// part4=neutron
  m3sq=m3*m3;
  m4sq=m4*m4;

  // Update tmin tmax. ATTENTION tmax is still... the usual tmin (=ft0)
  // tmin is the usual tmax (=ft1).. thank Carlos for that convention

  // fs is already calculated when you call this method !!!

  // CoM energies

  fecm1=(fs+m1sq-m2sq)/2./TMath::Sqrt(fs);
  fecm3=(fs+m3sq-m4sq)/2./TMath::Sqrt(fs);

  // CoM momentum

  fpcm1=TMath::Sqrt(fecm1*fecm1-m1sq);
  fpcm3=TMath::Sqrt(fecm3*fecm3-m3sq);

  // tmin/max calculations

  ft0=TMath::Power((m1sq-m3sq-m2sq+m4sq)/(2.*TMath::Sqrt(fs)),2.)-TMath::Power(fpcm1-fpcm3,2.);
  ft1=TMath::Power((m1sq-m3sq-m2sq+m4sq)/(2.*TMath::Sqrt(fs)),2.)-TMath::Power(fpcm1+fpcm3,2.);

  ftmax=TMath::Min(ftmax,ft0);
  ftmax-=0.0001;
  ftmin=TMath::Max(ftmin,ft1);
  ftmin+=0.0001;

//   cout <<"--------------------------------------------"<<endl;
//   cout << "fecm3="<<(fs+m3sq-m4sq)/2./TMath::Sqrt(fs)<<"    "<<"m3="<<TMath::Sqrt(m3sq)<<endl;
//   cout << "ecm4="<<(fs+m4sq-m3sq)/2./TMath::Sqrt(fs)<<"    "<<"m4="<<TMath::Sqrt(m4sq)<<endl;
//   cout << "mass="<<fdeltam<<endl;
//   cout << "t0  ="<<ft0<<"    "<<"t1  ="<<ft1<<endl;
//   cout << "tmax="<<ftmax<<"    "<<"tmin="<<ftmin<<endl;

}

//_____________________________________________________________________________
 void TGenPip::Compute2Body(Double_t m) 
{ 
  //Computes the gamma* p -> X n reaction in the center of mass, as a
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

  // Now for the "real" values of tmin and tmax
  ComputeTminTmax(m,NMass()); // (pip,neutron)
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

  //////////////////////////////////////////

//   Double_t ecm=(fs-fm*fm+m*m)/(2.*TMath::Sqrt(fs));
//   Double_t pcm=TMath::Sqrt(TMath::Power(ecm,2.)-TMath::Power(m,2.));
  Double_t thetacm=2.*TMath::ASin(TMath::Sqrt((ft0-ft)/4./fpcm1/fpcm3));

  if (!fpip) {fpip=new TLorentzVector(fpcm3*TMath::Sin(thetacm),0.,TMath::Cos(thetacm)*fpcm3,fecm3);}
  else {fpip->SetPxPyPzE(fpcm3*TMath::Sin(thetacm),0.,TMath::Cos(thetacm)*fpcm3,fecm3);}
  if (!fn) {fn=new TLorentzVector(q.Px()+p.Px()-fpip->Px(),q.Py()+p.Py()-fpip->Py(),q.Pz()+p.Pz()-fpip->Pz(),q.E()+p.E()-fpip->E());} 
      else {fn->SetPxPyPzE(q.Px()+p.Px()-fpip->Px(),q.Py()+p.Py()-fpip->Py(),q.Pz()+p.Pz()-fpip->Pz(),q.E()+p.E()-fpip->E());}
   
  // Rotation back
  fn->Rotate(-angle,perpvec.Unit());
  fpip->Rotate(-angle,perpvec.Unit());
  q.Rotate(-angle,perpvec.Unit());
  
  // Boost to lab
  fn->Boost(cms.BoostVector());
  fpip->Boost(cms.BoostVector());
  q.Boost(cms.BoostVector());

  // Now we have k along Oz and k' on the xOz plane
  // q=k-k' is necessarily on plane as well
  // We just have to rotate along Oy to bring q along Oz
  // and then rotate everything by phi except
  // k and k'. Then we rotate back along Oy.

  angle=(q.Vect()).Angle(oz);

  fn->RotateY(angle);
  fpip->RotateY(angle);
  
  fphi=2.*TMath::Pi()*fRan->Rndm();   // phi between 0 and 2pi
  
  fn->RotateZ(fphi);
  fpip->RotateZ(fphi);

  fn->RotateY(-angle);
  fpip->RotateY(-angle);

  // Pi rotation because the electron is always generated on the left

  fe->RotateZ(TMath::Pi());
  feprerad->RotateZ(TMath::Pi());
  fq->RotateZ(TMath::Pi());
  fpini->RotateZ(TMath::Pi());
  fn->RotateZ(TMath::Pi());
  fpip->RotateZ(TMath::Pi());


}

//_____________________________________________________________________________
void TGenPip::ApplySpecVerAcc(Double_t aav)
{
  // Applies vertical spectrometer acceptance by rotating all 4-vectors
  // around the beam axis. An angle can be specified, otherwise it's generated
  // randomly between spectrometer acceptances

  Double_t av=aav;
  if(aav==-1.) {
    av=-TMath::ATan(fCaloVerAccD/fCaloDist);
    av+=+(TMath::ATan(fCaloVerAccU/fCaloDist)+TMath::ATan(fCaloVerAccD/fCaloDist))*fRan->Rndm();
    av*=1/TMath::Abs(TMath::Sin(fCaloAngle));
  }

  TGenBase::ApplySpecVerAcc(av);
  fpip->RotateZ(av);
  fn->RotateZ(av);
  //fe->RotateZ(av);
  fq->RotateZ(av);
  
  fPSF*=(TMath::ATan(fCaloVerAccU/fCaloDist)+TMath::ATan(fCaloVerAccD/fCaloDist))/TMath::Abs(TMath::Sin(fCaloAngle));

  // Remark:
  // The right phasespace factor is
  // dx.dq2.dt.dphi.dphie/2pi (because the x-section code gives 
  // phi_e-integrated x-sections) but dphi=2pi, so it simplifies to
  // dx.dq2.dt.dphie
}

//_____________________________________________________________________________
 TLorentzVector* TGenPip::GetFinalNeutron(void)
{
  // Returns the final photon 4-vector if it exists
  
  if(!fn) cout<<"Warning : Final proton doesn't exist"<<endl;
  return fn;

}

//_____________________________________________________________________________
 TLorentzVector* TGenPip::GetFinalPip(void)
{
  // Returns the final photon 4-vector if it exists
  
  if(!fpip) cout<<"Warning : Final proton doesn't exist"<<endl;
  return fpip;

}

//_____________________________________________________________________________
void TGenPip::Write2File(void)
{

   if(!fNwrite) {
     if(!fFermi) {
       fNwrite=17;
     }else{
       fNwrite=17+3;
     }
   }
  
  *output<<fVertex->Pz()<<" ";
  *output<<feini->Pz()<<" ";
  *output<<feprerad->Px()<<" "<<feprerad->Py()<<" "<<feprerad->Pz()<<" ";
  if(fFermi)
    *output<<fpini->Px()<<" "<<fpini->Py()<<" "<<fpini->Pz()<<" ";
  *output<<fe->Px()<<" "<<fe->Py()<<" "<<fe->Pz()<<" ";
  *output<<fpip->Px()<<" "<<fpip->Py()<<" "<<fpip->Pz()<<" ";
  *output<<fn->Px()<<" "<<fn->Py()<<" "<<fn->Pz()<<" ";
  *output<<fPSF<<" "<<ft<<" "<<ftmax; 
  *output<<endl;

//  cout << "." << endl;
//   TVector3 qpt(fg1->Px()+fg2->Px(),fg1->Py()+fg2->Py(),fg1->Pz()+fg2->Pz());
//   TVector3 qt(-fe->Px(),-fe->Py(),feini->Pz()-fe->Pz());
//   TVector3 kt(0,0,feini->Pz());

//   TVector3 qk=qt.Cross(kt);//perp at leptonic plan (initial electron,virtual photon)

//   TVector3 qqp=qt.Cross(qpt);//perp at hadronic plan (final pi0, virtual photon)

//   Double_t ang=(180./TMath::Pi())*(qk.Angle(qqp));//angle between the 2 plans

//   if(qt.Dot(qk.Cross(qqp))<0) ang=360.-ang;

//      cout << "Phi_gen="<<fphi*180./3.141592 << endl;

//      cout << "pi0: (" << qpt.X()<<","<<qpt.Y()<<","<<qpt.Z()<<")"<<endl;

//    cout << "pi0 unit: (" << qpt.Unit().X()<<","<<qpt.Unit().Y()<<","<<qpt.Unit().Z()<<")"<<endl;
//    cout << "  k: (" << kt.X()<<","<<kt.Y()<<","<<kt.Z()<<")"<<endl;
//    cout << " kp: (" << fe->Px()<<","<<fe->Py()<<","<<fe->Pz()<<")"<<endl;
//    cout << "  q: (" << qt.X()<<","<<qt.Y()<<","<<qt.Z()<<")"<<endl;
//    cout << "  q unit: (" << qt.Unit().X()<<","<<qt.Unit().Y()<<","<<qt.Unit().Z()<<")"<<endl;
//     cout << " qk: (" << qk.X()<<","<<qk.Y()<<","<<qk.Z()<<")"<<endl;
//     cout << " qqp: (" << qqp.X()<<","<<qqp.Y()<<","<<qqp.Z()<<")"<<endl;
//     cout << "Phi_rec=" << ang << endl;
//      cout << "Q2="<<fQ2<<endl;
//      cout << "xb="<<fxb<<endl;
//      cout << "t="<<ft<<endl;
//      cout << "tmax="<<ftmax<<endl;
//     cout << "------------------" << endl;

    //cout << fphi*180./3.141592 << "   " << ang << endl;

 }
//___________________________________________________________________
 void TGenPip::Print(char* opt)
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
  if(fpip){
    cout<<"pip("<<fpip->Px()<<","<<fpip->Py()<<","<<fpip->Pz()<<","<<fpip->E()<<")"<<endl;
  }else{
    cout<<"NO EMITTED PI+ DEFINED"<<endl;
  }
  if(fn){
    cout<<"n'("<<fn->Px()<<","<<fn->Py()<<","<<fn->Pz()<<","<<fn->E()<<")"<<endl;
  }else{
    cout<<"NO RECOIL PARTICLE DEFINED"<<endl;
  }
  cout<<"======================================="<<endl;
}
