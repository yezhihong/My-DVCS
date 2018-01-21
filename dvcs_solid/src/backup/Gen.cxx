
#include <iostream>
#include <stdlib.h>
#include "TGenBase.h"
#include "TGenGeo.h"
#include "TGenDVCS.h"
#include "TObject.h"
/*ROOT Includes{{{*/
#include <TSystem.h>
#include <TString.h>
#include <TStyle.h>
#include <Riostream.h>
#include "TObjString.h"
#include <TNamed.h>
#include <TPRegexp.h>
#include <TObjArray.h>
#include <TChain.h>
#include <TMath.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TROOT.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TDatime.h>
#include <TError.h>
#include <TVirtualFitter.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TCut.h>
#include <TMultiGraph.h>
#include <TCutG.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>
//#include <TMatrix.h>
/*}}}*/

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

	TGenDVCS *gEv=new TGenDVCS(11.0,0,seed1,seed2); //target:0->P,1->N,2->D
	gEv->SetTargetParam(3.,0.78); // P target, DNP-NH3
	//gEv->SetTargetParam(40.,1.345); // N target, He3, note that N will be 1/3 from the luminosity

	//Geometry of spectrometer for electron and the calorimeter for photon detection
	//          (Spec_Angle, SpecMom, CaloAngle, CaloDist)
	//          EC momentum from 0.5~7.5GeV for SoLID which gives the central value of 4.0 with +/-3.5 range
 	gEv->SetGeometry(0.,4.0,0.,790.); 
	gEv->SetSpectroAcceptance(26.0*TMath::DegToRad(), 26.0*TMath::DegToRad(),3.5/4.0);
	//Only detect photons going forward from 6 to 16, where the actuall values is 8 to 14.8
	gEv->SetCaloAcceptance(16.0*TMath::DegToRad(),-16.0*TMath::DegToRad(),16.0*TMath::DegToRad(),-16.0*TMath::DegToRad());

	////////////////////////
	//gEv->SetFermi(kTRUE); // For Fermi only
	////////////////////////

	Int_t nev=10000000; // Good for Proton/Neutron
	//Int_t nev=1500000; // good for Deuteron
	Int_t j=0;
	cout<<"Generating "<<nev<<" events..."<<endl;

	for(Int_t i=0;i<nev;i++){
		if ((i%1000)==0) cout<<i<<"/"<<nev<<endl;
		gEv->Settmin(-2.); //Use this method to change tmin (default -2 GeV)
		gEv->Settmax(0.);  //Use this method to change tmax (default 0 GeV)
		gEv->GenerateVertex();//Generates the vertex
		gEv->ExtBrem();//Make external pre-vertex radiative corrections
		gEv->GenKinSoLID();//Computes scattered electron kinematics, 
		//-------------------for SoLID(wider), 6<Theta<26, 0.5<Ee<7.5, 0.05<xb<0.8, 0.1<Q2<10.0

		gEv->IntRCBef();//Internal radiative corrections (before vertex)
		/////////////////////
		//gEv->GenFermiIni();//To use if fFermi is set
		/////////////////////
		if(gEv->ComputeElectron()){//Compute the scattered electron 4-vector
			gEv->IntRCAft();//Internal radiative corrections (after vertex)
			if(gEv->HitsSpectro(gEv->GetScatteredElectron())){//Check hor spectro accep.
				gEv->ComputeDVCS();//Compute the gamma*p->gamma p collision
				gEv->ApplySpecVerAcc();//Rotates all final vectors around the beam axis
				if(gEv->HitsCalo(gEv->GetFinalPhoton())){//Checks calo acceptance
					//          gEv->XSec();
					//gEv->HitsCalo(gEv->GetFinalPhoton());
					j++;
					gEv->Write2File();//This writes the event in a temporary file
				}
			}
		}
	}
	gEv->CloseTmpFile();
	cout<<"j="<<j<<endl;

	gEv->DumpFinalFile("outdvcs.txt",nev,j);

	return 0;
}
