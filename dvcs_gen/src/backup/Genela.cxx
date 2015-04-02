#include "TLorentzVector.h"
#include <iostream>
#include <stdlib.h>
#include "TGenBase.h"
#include "TGenGeo.h"
#include "TGenElas.h"
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

  TGenElas *gEv=new TGenElas(5.7572,0,seed1,seed2);
  gEv->SetTargetParam(15.,0.78,0.0723);
  gEv->SetGeometry(0.65781,2.381,0.35774,5465.66);

  //gEv->SetFermi(kTRUE); // For Fermi only

  Int_t nev=100000;

  Int_t j=0;
  cout<<"Generating "<<nev<<" events..."<<endl;
  
  for(Int_t i=0;i<nev;i++){
    if ((i%100)==0) cout<<i<<"/"<<nev<<endl;
    gEv->GenerateVertex();//Generates the vertex
    gEv->ExtBrem();//Make external pre-vertex radiative corrections
    gEv->GenKin();//Computes scattered electron kinematics
    gEv->IntRCBef();//Internal radiative corrections (before vertex)
    if(gEv->ComputeElectron()){//Compute the scattered electron 4-vector
      gEv->IntRCAft();//Internal radiative corrections (after vertex)
      gEv->ComputeElas();//Compute the scattered electron 4-vector
      if(gEv->HitsSpectro(gEv->GetFinalProton())){//Check hor spectro accep.
        gEv->ApplySpecVerAcc();//Rotates all final vectors around the beam axis
        if( (gEv->HitsCalo(gEv->GetScatteredElectron()))){//Checks calo acceptance
	  //          gEv->XSec();
          j++;
          gEv->Write2File();//This writes the event in a temporary file
        }
      }
    }
  }
  cout<<"j="<<j<<endl;

  gEv->DumpFinalFile("outela.txt",nev,j);

  return 0;
}
