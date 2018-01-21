#include "TLorentzVector.h"
#include <iostream>
#include <stdlib.h>
#include "TGenBase.h"
#include "TGenGeo.h"
#include "TGenPi0.h"
#include "TObject.h"

using namespace std;

int main(int argc, char *argv[])
{
  cout << "Hello !" << endl;
  if(argc!=3) {
     cout << "Usage: gen <seed1> <seed2>" << endl;
     return -1;
  }
  
  Int_t seed1=atoi(argv[1]);
  Int_t seed2=atoi(argv[2]);

    TGenPi0 *gEv=new TGenPi0(5.7572,0,seed1,seed2); // P target
  //  TGenPi0 *gEv=new TGenPi0(5.7572,1,seed1,seed2); // N target
  //TGenPi0 *gEv=new TGenPi0(5.7572,2,seed1,seed2); // D target

  gEv->SetTargetParam(15.,0.78,0.0723); // P target
  //gEv->SetTargetParam(15.,0.78,0.167); // D target

  //gEv->SetGeometry(0.417308,2.35,0.25836,1119.98); // kin3
  //gEv->SetGeometry(0.3371976,2.95,0.318627,1119.88); // kin2
  gEv->SetGeometry(0.271922,3.55,0.388999,1119.11); // kin1

  //gEv->SetFermi(kTRUE); // For Fermi only

  Int_t nev=7500000; // Good for Proton/Neutron
  //Int_t nev=1000000; // good for Deuteron
  Int_t j=0;
  cout<<"Generating "<<nev<<" events..."<<endl;
  
  for(Int_t i=0;i<nev;i++){
    if ((i%100)==0) cout<<i<<"/"<<nev<<endl;
    gEv->Settmin(-1.); //Use this method to change tmin (default -2 GeV)
    gEv->Settmax(0.);  //Use this method to change tmax (default 0 GeV)
    gEv->GenerateVertex();//Generates the vertex
    gEv->ExtBrem();//Make external pre-vertex radiative corrections
    gEv->GenKin();//Computes scattered electron kinematics

    gEv->IntRCBef();//Internal radiative corrections (before vertex)
    /////////////////////
    //gEv->GenFermiIni();//To use if fFermi is set
    /////////////////////
    if(gEv->ComputeElectron()){//Compute the scattered electron 4-vector
      gEv->IntRCAft();//Internal radiative corrections (after vertex)
      if(gEv->HitsSpectro(gEv->GetScatteredElectron())){//Check hor spectro accep.
        gEv->ComputePi0();//Compute the gamma*p->gamma p collision
	gEv->TwoBodyDecay(0.1349766,0.,0.);
        gEv->ApplySpecVerAcc();//Rotates all final vectors around the beam axis
        if( (gEv->HitsCalo(gEv->GetFinalPhoton1())) || (gEv->HitsCalo(gEv->GetFinalPhoton2()))){//Checks calo acceptance
	  //          gEv->XSec();
          j++;
          gEv->Write2File();//This writes the event in a temporary file
        }
      }
    }
  }
  cout<<"j="<<j<<endl;

  gEv->DumpFinalFile("outpi0.txt",nev,j);

  return 0;
}
